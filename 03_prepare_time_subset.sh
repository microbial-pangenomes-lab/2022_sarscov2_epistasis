#!/bin/bash
set -eo pipefail

# in the format YYYY-MM
time=$1
cores=$2

mkdir -p data/time-subsets
xzcat prefiltered.fasta.xz | python3 src/combine_times.py data/time $time data/time-subsets
xz -T $cores data/time-subsets/$time.fasta
