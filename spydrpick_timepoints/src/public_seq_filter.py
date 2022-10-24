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
                        help='GISAID public IDs file (one ID per line)')

    parser.add_argument('-n',
                        default=None,
                        type=int,
                        help='Number of sequences to randomly keep '
                             '(default: keep them all)')

    return parser.parse_args()


def draw_random_indexes(num, maxcol):
    indexes_set = set([randint(0, maxcol) for _ in range(0,num)])
    # if duplicates are present, set will be shorter than the number of random numbers wanted: update the set object
    # till the length of set will be the same of the list of random numbers wanted.
    while len(indexes_set) != num:
        diff = num - len(indexes_set)
        indexes_set.update([randint(0,maxcol) for _ in range(0,diff)])
        #print(diff)
    return list(indexes_set)


if __name__ == '__main__':
    args = get_options()

    # Load tsv with sequences
    sequences = open(args.public)
    keep = {l.rstrip() for l in sequences.readlines()}

    # Load high quality seqs selected for tree in GISAID
    gisaid = open(args.metadata)
    gisaid_tree = {l.rstrip() for l in gisaid
                   if l.rstrip().split(',')[0] in keep}
    #print(gisaid_tree)
    selected = gisaid_tree
    names_dates = {}
    for x in selected:
        names_dates[x.split(',')[0]] = '-'.join(x.split(',')[2].split('-')[:1])
        print(x)
    pickle.dump(names_dates, open('names_dates.pkl','wb'))
