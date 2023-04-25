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
    all_ids = set()
    for i, fname in enumerate(inputs):
        sys.stderr.write(f'Reading IDs from {fname}\n')
        ids = read_fasta(f'{args.input}/{fname}')
        for x in ids:
            all_ids.add(x)

    f = open(f'{args.output}/{args.target}.fasta', 'w')
    seqs = 0
    s = []
    sid = ''
    kept = 0
    # Loop on all data
    for l in sys.stdin:
        if l.startswith('>') and len(s) != 0:
            if not seqs % 1000:
                sys.stderr.write(f'Split: {seqs} {kept}\n')
            seqs += 1
            if sid[1:] in all_ids:
                f.write(f'{sid}\n{"".join(s)}\n')
                kept += 1
                s = []
            sid = l.rstrip()
        elif l.startswith('>'):
            sid = l.rstrip()
            continue
        else:
            s.append(l.rstrip())
    if sid[1:] in all_ids and len(s) != 0:
        f.write(f'{sid}\n{"".join(s)}\n')
        kept += 1
    sys.stderr.write(f'Split: {seqs} {kept}\n')
    f.close()
