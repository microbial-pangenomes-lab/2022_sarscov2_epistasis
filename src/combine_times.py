#!/usr/bin/python


import os
import sys
import lzma
import argparse


def get_options():
    description = 'Combine time subsets'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('input',
                        help='Input directory (contains fasta files)')
    parser.add_argument('target',
                        help='Target time (format: YYYY-MM)')
    parser.add_argument('output',
                        help='Output directory')

    return parser.parse_args()


def read_fasta(fname):
    ids = set()
    for l in lzma.open(fname, 'rt'):
        if l.startswith('>'):
            ids.add(l.rstrip()[1:])
    return ids


if __name__ == '__main__':
    args = get_options()

    # sort input files
    inputs = sorted([x for x in os.listdir(args.input)
                     if x.endswith('.fasta.xz')],
                    key=lambda x: (int(x.split('-')[0]),
                                   int(x.split('-')[1].split('.')[0])))
    # restrict to desired last time point
    i = inputs.index(f'{args.target}.fasta.xz')
    inputs = inputs[:i+1]

    cmd = 'cat'
    for f in inputs:
        cmd += f' {args.input}/{f}'
    cmd += f' > {args.output}/{args.target}.fasta.xz'
    os.system(cmd)
