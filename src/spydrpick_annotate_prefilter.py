#!/usr/bin/env python


import sys
import argparse
import numpy as np
import pandas as pd


def get_options():
    description = 'Post-process the MI values'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('prefilter',
                        help='Prefiltered MI table')
    parser.add_argument('output',
                        help='Output MI table')

    parser.add_argument('--outliers',
                        default=None,
                        help='Outliers file (one value per line'
                             'default: do not report outliers)')

    return parser.parse_args()


def main():
    options = get_options()

    outliers = []
    if options.outliers is not None:
        tmp = []
        for v in open(options.outliers):
            tmp.append(float(v.rstrip()))
        outliers = [(i+1, v) for i, v in enumerate(sorted(tmp))]

    df = pd.read_csv(options.prefilter, sep='\t')
    df['distance'] = abs(df['pos_a'] - df['pos_b'])
    res = []
    for mi in df['mi'].values:
        outl = 0
        for i, v in outliers:
            if mi >= v:
                outl = i
            else:
                break
        res.append(outl)
    df['outlier'] = res
    df[['pos_a', 'pos_b', 'distance', 'outlier', 'mi']].to_csv(options.output, sep=' ', index=False, header=False)


if __name__ == "__main__":
    main()
