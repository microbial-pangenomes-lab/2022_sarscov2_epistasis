#!/usr/bin/python

import os
import sys
import random

# Load tsv with sequences
sequences = open('data/epi.tsv', 'r')
keep = [l.rstrip() for l in sequences.readlines()]

# Load high quality seqs selected for tree in GISAID
gisaid = open(sys.argv[1])
gisaid_tree = {l.rstrip().split(',')[0] for l in gisaid}

# if -n is indicated take the numeric value and
args = len(sys.argv)
#print(len(sys.argv))
if args > 1:
    flag1 = sys.argv[2]
    if flag1 == '-n' or flag1 == '--number':
        try:
            if sys.argv[3].isdigit():
                selected = sys.argv[3]
                sys.stderr.write(f'Number of sequences selected: {selected}\n')
            else:
                #sys.stderr.write('Expecting number, got other type of value\n')
                sys.exit()
        except:
            sys.stderr.write('Expecting numeric value, got None\n')
            sys.exit()
    elif flag1 == '-total':
        pass
    else:
        sys.stderr.write('Indicate -total or --number flags\n')
        sys.exit()
    keep_random = []
    #seq_selected = random.choices(keep, k=int(selected))

good_qual = []
for x in keep:
    if x in gisaid_tree:
        good_qual.append(x)
#print(good_qual)

if 'selected' in locals():
    selected_good = random.choices(good_qual, k=int(selected))
    for x in selected_good:
        print(x)
else:
    sys.stderr.write(f'Number of sequences selected: {len(good_qual)}\n')
    for x in good_qual:
        print(x)
