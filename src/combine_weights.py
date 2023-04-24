#!/usr/bin/env python


import argparse
import pandas as pd
import numpy as np


def get_options():
    description = 'Combine weights'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('method',
                        choices=['product', 'mean'],
                        help='Method to combine weights')
    parser.add_argument('weight',
                        nargs='+',
                        help='Weights file (one float per line)')
    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()
    # read all the weights
    weights = []
    for infile in options.weight:
        w = np.loadtxt(infile,
                       delimiter=' ', comments=None)
        weights.append(w)

    if options.method == 'product':
        out = np.prod(np.vstack(weights), axis=0)
    elif options.method == 'mean':
        out = np.mean(np.vstack(weights), axis=0)
    else:
        sys.stderr.write(f'ERROR: unattended case {options.method}\n')
        sys.exit(1)

    for w in out:
        print(w)
