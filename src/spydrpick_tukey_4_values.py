#!/usr/bin/env python

import sys
import argparse
import numpy as np


def tukey_outlier(stream):
    d = {}
    for l in stream:
        s = l.rstrip().split()
        if len(s) == 0:
            continue
        if s[0] == 'pos_a':
            continue
        a, b, mi = s
        a = int(a)
        b = int(b)
        mi = float(mi)
        d[a] = d.get(a, -1)
        d[b] = d.get(b, -1)
        if mi > d[a]:
            d[a] = mi
        if mi > d[b]:
            d[b] = mi
    max_mi = np.fromiter(d.values(), dtype=float)

    Q1, Q3 = np.quantile(max_mi, [0.25, 0.75])

    t1 = Q3 + 1.5 * (Q3 - Q1)
    t2 = Q3 + 3 * (Q3 - Q1)
    t3 = Q3 + 6 * (Q3 - Q1)
    t4 = Q3 + 12 * (Q3 - Q1)
    print(t1)
    print(t2)
    print(t3)
    print(t4)


def get_options():
    description = 'Get the Tukey outlier thresholds for MI values'
    parser = argparse.ArgumentParser(description=description)
    return parser.parse_args()

def main():
    tukey_outlier(sys.stdin)


if __name__ == "__main__":
    main()
