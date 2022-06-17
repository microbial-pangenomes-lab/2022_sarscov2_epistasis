#!/usr/bin/python3

import subprocess

# Run a subprocess and save the output to a file ---> the -c threshold was set to 0.99997 (= 1 - 1 / 29629 , for 1 SNP per sequence)
command = "mmseqs easy-linclust filtered_dca.fasta clusters /tmp/mmseqs2 --split-memory-limit 64G --threads 24 --kmer-per-seq-scale 0.3 -c 0.95 --min-seq-id 0.9989" # The --min-seq-id threshold is gived by considering 30 snps (a min of ~0.00033 dissimilarity between sequences considering the average length of the sequences 29629 and 10 snps for each seq)
# the calculus is 1 - (10/29629) to get the threshold
command2 = "mmseqs easy-linclust subset.fasta clusters /tmp/mmseqs2 --split-memory-limit 64G --threads 18 --kmer-per-seq-scale 0.3 -c 0.9995"
# out file

# start the process
with open('error_mmseqs.log', 'wt') as err:
    child = subprocess.Popen(command, stdin=None, stdout=err, stderr=None, shell=True)
    child.wait()
