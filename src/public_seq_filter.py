#!/usr/bin/python

import os
import sys
import random
import argparse

def get_options():
    description = 'Reduce GISAID to public sequences'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('metadata',
                        help='GISAID tree metadata file')
    parser.add_argument('public',
                        help='GISAID public IDs file (one ID per line)')

    parser.add_argument('-n',
                        default=None,
                        type=int,
                        help='Number of sequences to randomly keep '
                             '(default: keep them all)')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_options()

    # Load tsv with sequences
    sequences = open(args.public)
    keep = {l.rstrip() for l in sequences.readlines()}

    # Load high quality seqs selected for tree in GISAID
    gisaid = open(args.metadata)
    gisaid_tree = {l.rstrip().split(',')[0] for l in gisaid
                   if l.rstrip().split(',')[0] in keep}

    if arg.n is not None:
        selected = random.choices(gisaid_tree, k=args.n)
    else:
        selected = gisaid_tree

    for x in selected:
        print(x)
