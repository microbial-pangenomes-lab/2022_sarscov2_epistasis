#!/bin/bash
set -eo pipefail

tree=$1
mmsa=$2
#subset=$3

treedir=$(basename $tree .zip)
mmsadir=$(basename $mmsa .tar.xz)

echo $treedir
rm -rf $treedir || true
rm -rf $mmsadir || true

# Unzip GISAID tree and Multisequence Alignment
unzip $tree
tar -xvf $mmsa
xz -T 18 $mmsadir/*.fa


# generate list of public GISAID IDs
mkdir -p data
curl --silent -L https://github.com/CDCgov/SARS-CoV-2_Sequencing/raw/master/files/epiToPublic.tsv.gz | gzip -d | awk '{print $1}' > data/public_seqs.txt
# Extract the the public sequences from the total tree
python3 src/public_seq_filter.py $treedir/metadata.csv data/public_seqs.txt > data/tree_subset.txt
# Extract public sequences from GISAID multialignment
xzcat $mmsadir/*.fa.xz | python3 src/prefilter.py data/tree_subset.txt > prefiltered.fasta
xz -T 10 prefiltered.fasta

# Make dict and extract 27 subsets 1 per each month (5 months collapsed between 2019-12 and 2020-04) and randomly take 1000 seqs
# for each
xzcat prefiltered.fasta.xz | python3 src/mk_subsets.py
# Bash loop for concatenating incrementally all the files basing on date order
# Create first subset till april
cat 2019-12.fasta 2020-01.fasta 2020-02.fasta 2020-03.fasta 2020-04.fasta > tmp.fasta

# Delete non concatenated files
rm 2019-12.fasta 2020-01.fasta 2020-02.fasta 2020-03.fasta 2020-04.fasta

input='tmp.fasta'

for i in 20*; do
  echo $i
  output="${i}_out";
  echo $output;
  cat $input $i > $output;
  input="$output";
done

mv tmp.fasta 2020-04.fasta_out


#cd ../../../
# only keep those positions with starting from 266 to 29768
python3 src/make_pos.py > filtered_dca.fasta.pos
# Tree file in data folder
tree=${treedir}/global.tree

# Declarate array of months with '31' days
months=("2020-05" "2020-05" "2020-07" "2020-08" "2020-10" "2020-12" "2021-01" "2021-03" "2021-05" "2021-07" "2021-08" "2021-10" "2021-12" "2022-01" "2022-03" "2022-05")

# Iterate he workflow on each sample
for i in *.fasta_out;
do
  echo $i
  name=$(basename $i .fasta_out)
  IFS=- read -r year month <<< $name
  #echo $name
  new_path=$name
  #echo $new_path
  mkdir -p $new_path
  cat $i | python src/filter_single.py filtered_dca.fasta.pos $new_path/filtered.npz > $new_path/filtered_dca.fasta.names
  # Prune tree in order to get only needed leaves for each subset (TreeSwift needed)
  python3 src/extract_tree.py --tree $tree --name $name
  # Use BioPython to calculate weights based on pruned trees (relative sequence distances) for each subset
  python3 src/get_weights.py -s $new_path/filtered_dca.fasta.names --tree $name.tree > $new_path/filtered_dca_tree.fasta.weights
  mv $name.tree $new_path
  # Get weights based on time points
  # Each month has different number of days
  if [ $month == '02' ]; then
    day='28'; elif [[ "${months[*]}" =~ " ${name}" ]]; then day='31'; else day='30';
  fi
  python3 src/get_time_weight.py ${treedir}/metadata.csv $new_path/filtered_dca.fasta.names $year $month $day > $new_path/filtered_dca_time.fasta.weights
  # Combine the tree and time inferred weights using array product
  python $src/combine_weights.py product $new_path/filtered_dca_tree.fasta.weights $new_path/filtered_dca_time.fasta.weights > $new_path/product_weights.txt
  # Calculate filtered mi treshold
  python $src/spydrpick_alt.py -a $new_path/filtered.npz -p filtered_dca.fasta.pos -w $new_path/product_weights.txt > $new_path/filtered.mi_t
  # run spydrpick for real now on all samples and positions
  # each job focuses on a single positional window
  rm -rf $new_path/mi || true
  mkdir -p $new_path/mi
  for s in $(seq 0 1000 $(cat filtered_dca.fasta.pos | tr ',' '\n' | wc -l));
  do
    echo "python src/spydrpick_alt.py -a $new_path/filtered.npz -p filtered_dca.fasta.pos -w $new_path/product_weights.txt --threshold $(cat $new_path/filtered.mi_t) --start $s --cores 24 | gzip > $new_path/mi/$s.tsv.gz";
  done > $new_path/mi_jobs.txt
  parallel --jobs 1 --progress < $new_path/mi_jobs.txt

  # Calculate tukey outliers
  #zcat mi/*.tsv.gz | sort -n | uniq | python subset_workflow/spydrpick_tukey.py > subset_workflow/data/subsets/new_subset_total_out/new_weights_2022-06/mi_tukey.txt
  zcat $new_path/mi/*.tsv.gz | sort -n | uniq | python $src/spydrpick_tukey_4_values.py > $new_path/mi_tukey4.txt

  # Extract only outliers to perform final filter with aracne
  lower_t=$(cat $new_path/mi_tukey4.txt | head -1)
  #zcat mi/*.gz | awk -F '\t' -v t="$lower_t" '$3 >= t {print $0}' > subset_workflow/data/subsets/new_subset_total_out/new_weights_2022-06/mi_all_prefilter.txt 
  zcat $new_path/mi/*.gz | awk -F '\t' -v t="$lower_t" '$3 >= t {print $0}' > $new_path/mi_all_prefilter.txt

  # Run aracne
  cat $new_path/mi_all_prefilter.txt | sort -n | uniq | python src/spydrpick_filter.py --cores 1 --outliers $new_path/mi_tukey4.txt > $new_path/tmp.txt
  cat $new_path/tmp.txt | sort -n | uniq | gzip > $new_path/mi_all.tsv.gz

  # Calculate distance between positions
  python3 src/calc_distance.py $new_path/mi_all.tsv.gz
  mv out_spydrpick.txt $new_path/${name}.txt
  mv $i $name
done
