#!/usr/bin/python


import sys
import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def get_options():
    description = 'Subset and clean an alignment with a MAF filter'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('input')

    parser.add_argument('--maf', type=float, default=0.01,
                        help='MAF filter threshold '
                             '(%(default).2f)')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_options()

    a = AlignIO.read(args.input, 'fasta')

    keep = set()
    for i in range(len(a[0])):
        d = {}
        seq = a[:, i]
        for c in seq:
            d[c] = d.get(c, 0)
            d[c] += 1
        if sum([v for k, v in d.items() if v != max(d.values())]) / sum(d.values()) <= args.maf:
            keep.add(i)

    for record in a:
        seq = ''.join(record[i] for i in keep)
        print(f'>{record.id}\n{seq}')
