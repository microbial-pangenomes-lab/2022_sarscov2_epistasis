#!/bin/bash
set -eo pipefail

cores=$1
metadata=$2
# in the format YYYY-MM
time=$3

mkdir -p data/subsets
xzcat prefiltered.fasta.xz | python3 src/split_randomly.py data/tree_subset.txt data/subsets 1000 10000 100000 1000000
xz -T $cores data/subsets/*.fasta

mkdir -p data/time
xzcat prefiltered.fasta.xz | python3 src/split_by_time.py data/tree_subset.txt data/time/ --last $time
xz -T $cores data/time/*.fasta

mkdir -p data/time-10000
xzcat prefiltered.fasta.xz | python3 src/split_by_time.py data/tree_subset.txt data/time-10000/ --last $time -n 10000
xz -T $cores data/time-10000/*.fasta

nextclade dataset get -n sars-cov-2 -o data/nextclade

mkdir -p out/time-lineages
for i in $(ls data/time/ | grep fasta.xz);
do
  echo "nextclade run data/time/$i --output-tsv out/time-lineages/$(basename $i .fasta.xz).tsv --input-dataset data/nextclade/ -j 1";
done > nextclade_jobs.txt
parallel -j $cores --progress < nextclade_jobs.txt
rm nextclade_jobs.txt

mkdir -p data/time-filtered/
cp covariants/scripts/bad_sequences.py src/
python3 src/filter_times.py data/time $metadata data/time-filtered/

mkdir -p out/time-10000-lineages
for i in $(ls data/time-10000/ | grep fasta.xz);
do
  echo "nextclade run data/time-10000/$i --output-tsv out/time-10000-lineages/$(basename $i .fasta.xz).tsv --input-dataset data/nextclade/ -j 1";
done > nextclade_jobs.txt
parallel -j $cores --progress < nextclade_jobs.txt
rm nextclade_jobs.txt

mkdir -p data/time-10000-filtered/
cp covariants/scripts/bad_sequences.py src/
python3 src/filter_times.py data/time-10000 $metadata data/time-10000-filtered/
