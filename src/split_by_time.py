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
    parser.add_argument('-n',
                        type=int,
                        default=2500,
                        help='How many sequences to keep maximum (%(default)d)')
    parser.add_argument('--initial',
                        default='2020-02',
                        help='Put sequences together from up to this '
                              'month (%(default)s)')
    parser.add_argument('--last',
                        default='2023-03',
                        help='Last month '
                              'month (%(default)s)')
    parser.add_argument('--seed',
                        type=int,
                        default=100,
                        help='Random integer seed (%(default)d)')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_options()

    # get n random IDs per month
    m = pd.read_csv(args.metadata, header=None)
    m['date'] = pd.to_datetime(m[2])
    m['date'] = m['date'].dt.to_period('M')
    m = m[[0, 'date']].groupby('date').apply(
            lambda x: x.sample(n=min(args.n, x.shape[0]),
                random_state=args.seed)[0]
            ).reset_index()[['date', 0]]
    m = m[m['date'] <= args.last]
    ids = {}
    ids[args.initial] = set(m[m['date'] <= args.initial][0].values)
    for date in m[m['date'] > args.initial]['date'].unique():
        ids[str(date)] = set(m[m['date'] == date][0].values)
    dates = {x: y for y, xx in ids.items()
             for x in xx}

    files = {x: open(f'{args.output}/{x}.fasta', 'w')
             for x in ids}
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
            if sid[1:] in dates:
                f = files[dates[sid[1:]]]
                f.write(f'{sid}\n{"".join(s)}\n')
                kept += 1
                s = []
            sid = l.rstrip()
        elif l.startswith('>'):
            sid = l.rstrip()
            continue
        else:
            s.append(l.rstrip())
    if sid[1:] in dates and len(s) != 0:
        f = files[dates[sid[1:]]]
        f.write(f'{sid}\n{"".join(s)}\n')
        kept += 1
    sys.stderr.write(f'Split: {seqs} {kept}\n')
    for f in files.values():
        f.close()
