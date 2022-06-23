#!/bin/bash
set -euo pipefail

tree=$1
mmsa=$2
subset=$3

treedir=$(basename $tree .zip)
mmsadir=$(basename $mmsa .tar.xz)

echo $treedir
rm -rf $treedir || true
rm -rf $mmsadir || true

# Unzip GISAID tree and Multisequence Alignment
unzip $tree
tar -xvf $mmsa

# generate list of public GISAID IDs
curl --silent -L https://github.com/CDCgov/SARS-CoV-2_Sequencing/raw/master/files/epiToPublic.tsv.gz | gzip -d | awk '{print $1}' > data/public_seqs.txt

# If indicated in command line (3rd argument) select random n sequences from filtered gisaid dataset with public sequences
if [ $3 > 0 ]; then
  python3 src/public_seq_filter.py $treedir/metadata.csv data/public_seqs.txt -n $subset > tree_subset.txt
  xzcat $mmsadir/*_masked.fa.xz | python3 src/prefilter.py tree_subset.txt > prefiltered.fasta
else
# prefilter the sequences to keep only those in the metadata file
  python3 src/public_seq_filter.py $treedir/metadata.csv data/public_seqs.txt > tree_subset.txt
  xzcat $mmsadir/*_masked.fa.xz | python3 src/prefilter.py tree_subset.txt > prefiltered.fasta
fi

# only keep those positions with starting from 266 to 29768
python3 src/make_pos.py > filtered_dca.fasta.pos

# deduplicate sequences at the positions that need to be kept
xzcat prefiltered.fasta.xz | python src/filter_single.py filtered_dca.fasta.pos filtered.npz > filtered_dca.fasta.names

# remove gaps from original sequences, keep only those flagged at previous step
xzcat prefiltered.fasta.xz | pypy3 src/post_filter.py filtered_dca.fasta.names > filtered_dca.fasta

# cluster all deduplicated sequences at a defined similarity threshold
# mmseqs easy-linclust filtered_dca.fasta clusters /tmp/mmseqs2 --split-memory-limit 64G --threads 18 --kmer-per-seq-scale 0.3 -c 0.95 --min-seq-id 0.9996
python3 src/run_mmseqs.py
#rm filtered_dca.fasta.xz || true
mv clusters_cluster.tsv cluster.tsv
xz -T 18 filtered_dca.fasta
rm clusters_rep_seq.fasta
rm clusters_all_seqs.fasta

# assign weights to each sequence, proportional to how many sequences are in each cluster
python src/get_weights_linclust.py -s filtered_dca.fasta.names --clusters cluster.tsv > filtered_dca.fasta.weights
#python get_weights.py -s filtered_dca.fasta.names --tree GISAID-hCoV-19-phylogeny-2021-07-08/global.tree > filtered_dca.fasta.weights

# run spydrpick algorithm once to derive a first threshold
python src/spydrpick_alt.py -a filtered.npz -p filtered_dca.fasta.pos -w filtered_dca.fasta.weights --cores 18 > filtered.mi_t

# run spydrpick for real now on all samples and positions
# each job focuses on a single positional window
rm -rf mi || true
mkdir -p mi
for i in $(seq 0 1000 $(cat filtered_dca.fasta.pos | tr ',' '\n' | wc -l));
do
  echo "python src/spydrpick_alt.py -a filtered.npz -p filtered_dca.fasta.pos -w filtered_dca.fasta.weights --cores 18 --threshold $(cat filtered.mi_t) --start $i --chunk 1000 --window 1000 | gzip > mi/$i.tsv.gz";
done > mi_jobs.txt
parallel --jobs 1 --progress < mi_jobs.txt

# compute tukey outlier threshold
zcat mi/*.tsv.gz | python src/spydrpick_tukey.py > mi_tukey.txt
# postprocess results, run ARACNE
zcat mi/*.tsv.gz | python src/spydrpick_filter.py --cores 18 --outliers mi_tukey.txt --chunk 1000 > tmp.txt
# put it all together
sort -n tmp.txt | uniq | gzip > mi_all.tsv.gz

rm tmp.txt
python3 src/calc_distance.py mi_all.tsv.gz
rm -rf $treedir
rm -rf $mmsadir

# Plot and tables from Gwes results + NCBI public SNPs dataset and outbreak.info data
#python3 src/plot_gwes_results.py
