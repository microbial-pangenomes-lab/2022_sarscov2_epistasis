#!/bin/bash
set -eo pipefail

fasta=$1
idir=$2
odir=$3
cores=$4

mkdir -p data/$idir
mkdir -p out/$idir
mkdir -p $odir

name=$(basename $fasta .fasta.xz)

mkdir -p $odir/$name

xzcat $fasta | python3 src/subset_alignment.py --start 21562 --end 25384 > data/$idir/$name.fasta

plmc/bin/plmc -c $odir/$name.EC -o $odir/$name.params --fast -m 20 -le 20.0 -lh 0.01 -a -AGCT -n $cores -f EPI_ISL_402130 data/$idir/$name.fasta
python3 src/annotate_plmc.py data/GCF_009858895.2_ASM985889v3_genomic.gff.gz $odir/$name.EC.gz $odir/$name/plmc_prefilter.tsv.gz --top-scoring 10000
python3 src/plmc_filter.py $odir/$name/plmc_prefilter.tsv.gz $odir/$name/plmc.tsv.gz
