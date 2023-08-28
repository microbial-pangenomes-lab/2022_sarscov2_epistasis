#!/usr/bin/python


import os
import sys
import lzma
import argparse
from bad_sequences import bad_seqs


def get_options():
    description = 'Filter time subsets for low quality sequences'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('input',
                        help='Input directory (contains fasta files)')
    parser.add_argument('metadata',
                        help='metadata table')
    parser.add_argument('output',
                        help='Output directory')

    return parser.parse_args()


def read_fasta(fname):
    ids = set()
    for l in lzma.open(fname, 'rt'):
        if l.startswith('>'):
            ids.add(l.rstrip()[1:])
    return ids


def parse_fasta(fname):
    sid = None
    seq = ''
    for l in lzma.open(fname, 'rt'):
        if l.startswith('>'):
            if sid is not None and seq != '':
                yield sid, seq
            seq = ''
            sid = l.rstrip()[1:]
        else:
            seq += l.rstrip()

    if sid is not None and seq != '':
        yield sid, seq


if __name__ == '__main__':
    args = get_options()

    # sort input files
    inputs = sorted([x for x in os.listdir(args.input)
                     if x.endswith('.fasta.xz')
                     and x not in os.listdir(args.output)],
                    key=lambda x: (int(x.split('-')[0]),
                                   int(x.split('-')[1].split('.')[0])))
    all_ids = {y for x in inputs
                 for y in read_fasta(os.path.join(args.input, x))}

    # read metadata
    d = {}
    for l in open(args.metadata):
        sid, sacc, _ = l.rstrip().split(',')
        if sid in all_ids:
            d[sid] = sacc.split('hCoV-19/')[1].replace(' ', '')

    # cycle through each file and write the filtered version
    already = os.listdir(args.output)
    for f in inputs:
        if f in already:
            sys.stderr.write(f'skipping {f}\n')
            continue
        fout = lzma.open(os.path.join(args.output, f), 'wt')
        for sid, seq in parse_fasta(os.path.join(args.input, f)):
            sacc = d.get(sid)
            if sacc is None or sacc in bad_seqs:
                sys.stderr.write(f'skipping {sid} from {f}\n')
                continue
            fout.write(f'>{sid}\n{seq}\n')
        fout.close()
