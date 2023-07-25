#!/bin/bash
set -eo pipefail

idir=$1
odir=$2
cores=$3

mkdir -p $odir

nextclade dataset get -n sars-cov-2 -o data/nextclade

for i in $(ls $idir | grep fasta.xz);
do
  echo "nextclade run $idir/$i --output-tsv $odir/$(basename $i .fasta.xz).tsv --input-dataset data/nextclade/ -j 1";
done > nextclade_jobs.txt
parallel -j $cores --progress < nextclade_jobs.txt
rm nextclade_jobs.txt
