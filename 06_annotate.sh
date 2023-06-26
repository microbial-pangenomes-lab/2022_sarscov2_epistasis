#!/bin/bash
set -eo pipefail

odir=$1
cores=$2

rname=/tmp/$(uuidgen | tr -d '-').txt

for i in $(find $odir -type f -name mi_all_distances.tsv.gz);
do
  echo python3 src/annotate_mi.py data/GCF_009858895.2_ASM985889v3_genomic.gff.gz $i $(dirname $i)/mi_annotated.tsv.gz;
done > $rname
parallel -j $cores --progress < $rname
rm $rname
