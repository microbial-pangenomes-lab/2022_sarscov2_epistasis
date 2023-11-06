#!/bin/bash
set -eo pipefail

fasta=$1
odir=$2
cores=$3

mkdir -p data/spikes-10000
mkdir -p out/spikes-10000
mkdir -p $odir

name=$(basename $fasta .fasta.xz)

xzcat $fasta | python3 src/subset_alignment.py --start 21562 --end 25384 > data/spikes-10000/$name.fasta

plmc/bin/plmc -c $odir/$name.EC -o $odir/$name.params --fast -m 20 -le 20.0 -lh 0.01 -a -AGCT -n $cores -f EPI_ISL_402130 data/spikes-10000/$name.fasta
python3 src/annotate_plmc.py data/GCF_009858895.2_ASM985889v3_genomic.gff.gz $odir/$name.EC $odir/plmc_prefilter.tsv.gz
python3 src/plmc_filter.py $odir/plmc_prefilter.tsv.gz $odir/plmc.tsv.gz
