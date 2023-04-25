#!/bin/bash
set -eo pipefail

# in the format YYYY-MM
time=$1

mkdir -p data/time-subsets
python3 src/combine_times.py data/time $time data/time-subsets
