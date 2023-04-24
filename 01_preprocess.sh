#!/bin/bash
set -eo pipefail

tree=$1
mmsa=$2
cores=$3

treedir=$(basename $tree .zip)
mmsadir=$(basename $mmsa .tar.xz)

echo $treedir
rm -rf $treedir || true
rm -rf $mmsadir || true

# Unzip GISAID tree and Multisequence Alignment
unzip $tree
tar -xvf $mmsa
xz -T $cores $mmsadir/*.fa

# generate list of public GISAID IDs
mkdir -p data
curl --silent -L https://hgwdev.gi.ucsc.edu/~angie/epiToPublicAndDate.latest | awk '{print $1}' > data/public_seqs.txt
# Extract the the public sequences from the total tree
python3 src/public_seq_filter.py $treedir/metadata.csv data/public_seqs.txt > data/tree_subset.txt
# Extract public sequences from GISAID multialignment
xzcat $mmsadir/*.fa.xz | python3 src/prefilter.py data/tree_subset.txt > prefiltered.fasta
xz -T $cores prefiltered.fasta
