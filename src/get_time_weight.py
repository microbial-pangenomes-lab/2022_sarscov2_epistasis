#!/usr/bin/env python


import sys
import argparse
import numpy as np
import pandas as pd


def get_options():
    description = 'Weight sequences based on time'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('metadata',
                        help='GISAID metadata table')
    parser.add_argument('samples',
                        help='Entries to compute weights for (one per line)')
    parser.add_argument('year',
                        nargs='?', default=None,
                        help='Year to consider as last time point')
    parser.add_argument('month',
                        nargs='?', default=None,
                        help='Month to consider as last time point')
    parser.add_argument('day',
                        nargs='?', default=None,
                        help='Day to consider as last time point')

    parser.add_argument('-c',
                        type=int,
                        default=120,
                        help='Day at which to give a weight of 0.5 (default: %(default)d)')
    parser.add_argument('-d',
                        type=float,
                        default=3.,
                        help='Steepness of the Hill curve (default: %(default).1f)')
    return parser.parse_args()


def hill_func(x, a, b, c, d):
    """Hill function

    Args:
        x (numpy.array)
        a (float)
          minimum
        b (float)
          maximum
        c (float)
          value of x where y = 0.5
        d (float)
          "steepness" of the curve

    Returns:
        y (float)
    """
    return a+(b-a)/(1+(x/c)**d)


if __name__ == "__main__":
    options = get_options()

    # time to count days from
    if options.year is None or options.month is None or options.day is None:
        # workaround to implement the analysis without this time correction
        end = None
    else:
        end = pd.to_datetime(f'{options.year}-{options.month}-{options.day}')

    # read entries to consider
    entries = [x.rstrip() for x in open(options.samples)]

    # read metadata, restrict to entries we actually want
    m = pd.read_csv(options.metadata, sep=',', index_col=0)
    m = m.loc[entries]

    if end is not None:
        # convert date to usable format
        m['date'] = pd.to_datetime(m['collection_date'],
                                   format='mixed')

        # compute how many days in the past each sequence is from
        m['delta'] = end - m['date']
        m['days'] = np.float64(m['delta']) / 60 / 60 / 24 / 1_000_000_000

        # check if there are sequences in the future and exit with error
        future = m[m['days'] < 0].shape[0]
        if future > 0:
            sys.stderr.write(f'WARNING: there are {future} entries from the future\n')

        # apply Hill function to derive weigth
        m['weight'] = hill_func(m['days'], 0, 1, options.c, options.d)
        # give a weight of one to the samples from the "future"
        m.loc[m[m['days'] < 0].index, 'weight'] = 1
    else:
        # workaround
        m['weight'] = 1

    first = True
    for w in m['weight'].values:
        if first:
            first = False
            sys.stdout.write(str(w))
            continue
        sys.stdout.write(' ')
        sys.stdout.write(str(w))
