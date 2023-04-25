#!/usr/bin/python


import os
import sys
import argparse


def get_options():
    description = 'Combine time subsets'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('input',
                        help='Input directory (contains fasta files)')
    parser.add_argument('output',
                        help='Output directory (must be present)')

    return parser.parse_args()


def read_fasta(fname):
    ids = set()
    for l in open(fname):
        if l.startswith('>'):
            ids.add(l.rstrip()[1:])
    return ids


if __name__ == '__main__':
    args = get_options()

    # sort input files
    inputs = sorted([x for x in os.listdir(args.input)
                     if x.endswith('.fasta')],
                    key=lambda x: (int(x.split('-')[0]),
                                   int(x.split('-')[1].split('.')[0])))
    s2dates = {}
    for i, fname in enumerate(inputs):
        sys.stderr.write(f'Reading IDs from {fname}\n')
        ids = read_fasta(f'{args.input}/{fname}')
        for x in ids:
            s2dates[x] = s2dates.get(x, set())
            s2dates[x].add(fname)
            for j in range(i, len(inputs)):
                s2dates[x].add(inputs[j])

    files = {x: open(f'{args.output}/{x}', 'w')
             for x in inputs}
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
            if sid[1:] in s2dates:
                for fname in s2dates[sid[1:]]:
                    f = files[fname]
                    f.write(f'{sid}\n{"".join(s)}\n')
                kept += 1
                s = []
            sid = l.rstrip()
        elif l.startswith('>'):
            sid = l.rstrip()
            continue
        else:
            s.append(l.rstrip())
    if sid[1:] in s2dates and len(s) != 0:
        for fname in s2dates[sid[1:]]:
            f = files[fname]
            f.write(f'{sid}\n{"".join(s)}\n')
        kept += 1
    sys.stderr.write(f'Split: {seqs} {kept}\n')
    for f in files.values():
        f.close()
