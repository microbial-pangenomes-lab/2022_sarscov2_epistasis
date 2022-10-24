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

# Conserving only positions from 266 to 29768
python3 src/make_pos.py > filtered_dca.fasta.pos

# Possibly add a first step of clustering ? Also with the prefilter and public sequence filter the matrix resulting from the
# further step would be massive.

# Filter sequences removing duplicates and converting it in a numpy matrix for further processing
xzcat prefiltered.fasta.xz | python src/filter_single.py filtered_dca.fasta.pos filtered.npz > filtered_dca.fasta.names


tree='GISAID-hCoV-19-phylogeny-2022-07-05/global.tree'

# Prune tree in order to get only needed leaves for each subset (TreeSwift needed)
python3 src/extract_tree.py --tree $tree --name 'filtered_dca.fasta'

# Use BioPython to calculate weights based on pruned trees (relative sequence distances) for each subset
python3 src/get_weights.py -s filtered_dca.fasta.names --tree filtered_dca.fasta.tree > filtered_dca_tree.fasta.weights

# Get weights based on time points
python3 src/get_time_weight.py ${treedir}/metadata.csv filtered_dca.fasta.names 2022 06 30 > filtered_dca_time.fasta.weights
# Combine the tree and time inferred weights using array product
python src/combine_weights.py product filtered_dca_tree.fasta.weights filtered_dca_time.fasta.weights > product_weights.txt

# Run spydrpick once to calculate mi treshold
python src/spydrpick_alt.py -a filtered.npz -p filtered_dca.fasta.pos -w product_weights.txt > filtered.mi_t

# Create mi folder if it doesn't exist yet
rm -rf mi || true
mkdir -p mi
# Create multiple jobs and run spydrpick for real with chunks of 1000 nucleotides
for s in $(seq 0 1000 $(cat filtered_dca.fasta.pos | tr ',' '\n' | wc -l));
do
  echo "python src/spydrpick_alt.py -a filtered.npz -p filtered_dca.fasta.pos -w product_weights.txt --threshold $(cat filtered.mi_t) --start $s --cores 24 | gzip > mi/$s.tsv.gz"
done > mi_jobs.txt
#parallel --jobs 1 --progress < mi_jobs.txt

# Calculate tukey outliers
zcat mi/*.tsv.gz | sort -n | uniq | python src/spydrpick_tukey_4_values.py > mi_tukey4.txt

# Extract only outliers to perform final filter with aracne
lower_t=$(cat mi_tukey4.txt | head -1)
zcat mi/*.gz | sort -n | uniq | awk -F '\t' -v t="$lower_t" '$3 >= t {print $0}' > mi_all_prefilter.txt
# Run aracne
cat mi_all_prefilter.txt | sort -n | uniq | python src/spydrpick_filter.py --cores 1 --outliers mi_tukey4.txt > tmp.txt
cat tmp.txt | sort -n | uniq | gzip > mi_all.tsv.gz

# Calculate distance between positions
python3 src/calc_distance.py mi_all.tsv.gz

rm tmp.txt
