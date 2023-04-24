#!/usr/bin/python

import os
import sys
import random
import argparse
import pickle


def get_options():
    description = 'Reduce GISAID to public sequences'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('metadata',
                        help='GISAID tree metadata file')
    parser.add_argument('public',
                        default=None,
                        help='GISAID public IDs file (one ID per line)')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_options()

    # Load tsv with sequences
    if args.public is not None:
        sequences = open(args.public)
        keep = {l.rstrip() for l in sequences.readlines()}

    # Load high quality seqs selected for tree in GISAID
    gisaid = open(args.metadata)
    if args.public is not None:
        gisaid_tree = {l.rstrip() for l in gisaid
                       if l.rstrip().split(',')[0] in keep}
    else:
        gisaid_tree = {l.rstrip() for l in gisaid}
    #print(gisaid_tree)
    selected = gisaid_tree
    names_dates = {}
    for x in selected:
        names_dates[x.split(',')[0]] = '-'.join(x.split(',')[2].split('-')[:1])
        print(x)
