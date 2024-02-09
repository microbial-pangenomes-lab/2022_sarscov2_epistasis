#!/bin/bash
set -eo pipefail

# in the format YYYY-MM
time=$1

mkdir -p data/time-subsets
python3 src/combine_times.py data/time $time data/time-subsets

mkdir -p data/time-filtered-subsets
python3 src/combine_times.py data/time-filtered $time data/time-filtered-subsets

mkdir -p data/time-10000-filtered-subsets
python3 src/combine_times.py data/time-10000-filtered $time data/time-10000-filtered-subsets

mkdir -p data/time-100k-filtered-subsets
python3 src/combine_times.py data/time-100k-filtered $time data/time-100k-filtered-subsets
