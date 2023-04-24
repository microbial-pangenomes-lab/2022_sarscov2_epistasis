#!/usr/bin/python


import sys
import argparse
import pandas as pd


def get_options():
    description = 'Split GISAID by time'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('metadata',
                        help='Metadata file (3 columns, last one is the date)')
    parser.add_argument('output',
                        help='Output directory (must be present)')
    parser.add_argument('n',
                        type=int,
                        nargs='+',
                        help='How many sequences to keep')
    parser.add_argument('--seed',
                        type=int,
                        default=100,
                        help='Random integer seed (%(default)d)')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_options()

    # get n random IDs per month
    m = pd.read_csv(args.metadata, header=None)

    ids = {x: set(m.sample(n=x, random_state=args.seed)[0].values)
           for x in args.n}
    all_ids = {x for y in ids.values()
               for x in y}

    files = {x: open(f'{args.output}/{x}.fasta', 'w')
             for x in args.n}
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
                ff = [files[x] for x, y in ids.items()
                      if sid[1:] in y]
                for f in ff:
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
        f = files[dates[sid[1:]]]
        f.write(f'{sid}\n{"".join(s)}\n')
        kept += 1
    sys.stderr.write(f'Split: {seqs} {kept}\n')
    for f in files.values():
        f.close()
