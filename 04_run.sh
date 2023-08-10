#!/bin/bash
set -eo pipefail

tree=$1
cores=$2
fasta=$3
wdir=$4
# in the format YEAR-MONTH
time=$5

treedir=$(basename $tree .zip)

mkdir -p $wdir

# Conserving only positions from 266 to 29768
python3 src/make_pos.py > $wdir/filtered_dca.fasta.pos

# Filter sequences removing duplicates and converting it in a numpy matrix for further processing
xzcat $fasta | python src/filter_single.py $wdir/filtered_dca.fasta.pos $wdir/filtered.npz > $wdir/filtered_dca.fasta.names

treefile=${treedir}/global.tree

# Prune tree in order to get only needed leaves for each subset (TreeSwift needed)
python3 src/extract_tree.py --tree $treefile --name $wdir/filtered_dca.fasta.names --out $wdir/filtered_dca.fasta.tree

# Use dendropy to calculate weights based on pruned trees (relative sequence distances) for each subset
python3 src/get_weights.py -s $wdir/filtered_dca.fasta.names --tree $wdir/filtered_dca.fasta.tree > $wdir/filtered_dca_tree.fasta.weights

endtime=$(python3 src/get_end_time.py $time)
# Get weights based on time points
python3 src/get_time_weight.py ${treedir}/metadata.csv $wdir/filtered_dca.fasta.names $endtime > $wdir/filtered_dca_time.fasta.weights
# Combine the tree and time inferred weights using array product
python src/combine_weights.py product $wdir/filtered_dca_tree.fasta.weights $wdir/filtered_dca_time.fasta.weights > $wdir/product_weights.txt

# Run spydrpick once to calculate mi treshold
python src/spydrpick_alt.py -a $wdir/filtered.npz -p $wdir/filtered_dca.fasta.pos -w $wdir/product_weights.txt --cores $cores > $wdir/filtered.mi_t

# Create mi folder if it doesn't exist yet
rm -rf $wdir/mi || true
mkdir -p $wdir/mi
# Create multiple jobs and run spydrpick for real with chunks of 1000 nucleotides
for s in $(seq 0 1000 $(cat $wdir/filtered_dca.fasta.pos | tr ',' '\n' | wc -l));
do
  python src/spydrpick_alt.py -a $wdir/filtered.npz -p $wdir/filtered_dca.fasta.pos -w $wdir/product_weights.txt --threshold $(cat $wdir/filtered.mi_t) --start $s --cores $cores | gzip > $wdir/mi/$s.tsv.gz
done

# Calculate tukey outliers
zcat $wdir/mi/*.tsv.gz | sort -n | uniq | python src/spydrpick_tukey_4_values.py > $wdir/mi_tukey4.txt

# Extract only outliers to perform final filter with aracne
lower_t=$(cat $wdir/mi_tukey4.txt | head -1)
zcat $wdir/mi/*.gz | sort -n | uniq | awk -F '\t' -v t="$lower_t" '$3 >= t {print $0}' > $wdir/mi_all_prefilter.txt
# Run aracne
cat $wdir/mi_all_prefilter.txt | sort -n | uniq | python src/spydrpick_filter.py --cores 1 --outliers $wdir/mi_tukey4.txt > $wdir/tmp.txt
cat $wdir/tmp.txt | sort -n | uniq | gzip > $wdir/mi_all.tsv.gz
gzip $wdir/mi_all_prefilter.txt
# annotate the prefiltered version as well
python3 src/spydrpick_annotate_prefilter.py $wdir/mi_all_prefilter.txt.gz $wdir/mi_all_prefilter_distances.txt.gz --outliers $wdir/mi_tukey4.txt
python3 src/annotate_mi.py data/GCF_009858895.2_ASM985889v3_genomic.gff.gz $wdir/mi_all_prefilter_distances.txt.gz $wdir/mi_annotated_prefilter.tsv.gz;

# Calculate distance between positions
python3 src/calc_distance.py $wdir/mi_all.tsv.gz $wdir/mi_all_distances.tsv
gzip $wdir/mi_all_distances.tsv

rm $wdir/tmp.txt

python3 src/annotate_mi.py data/GCF_009858895.2_ASM985889v3_genomic.gff.gz $wdir/mi_all_distances.tsv.gz $wdir/mi_annotated.tsv.gz;
