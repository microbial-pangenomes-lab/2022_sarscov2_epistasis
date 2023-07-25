#!/bin/bash
set -eo pipefail

mmsa=$1
prefiltered=$2
wdir=$3

an=$(xzcat $mmsa | grep -c ">")
pn=$(wc -l data/public_seqs.txt)
tn=$(wc -l data/tree_subset.txt)
pfn=$(xzcat $prefiltered | grep -c ">")
f=$(wc -l wdir/filtered_dca.fasta.names)

echo $an $pn $tn $pfn $f
