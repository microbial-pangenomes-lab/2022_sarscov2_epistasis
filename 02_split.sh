#!/bin/bash
set -eo pipefail

cores=$1

mkdir -p data/subsets
xzcat prefiltered.fasta.xz | python3 src/split_randomly.py data/tree_subset.txt data/subsets 1000 10000 100000 1000000
xz -T $cores data/subsets/*.fasta

mkdir -p data/time
xzcat prefiltered.fasta.xz | python3 src/split_by_time.py data/tree_subset.txt data/time/
xz -T $cores data/time/*.fasta
