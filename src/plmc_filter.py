#!/usr/bin/env python


import sys
import tqdm
import argparse
import itertools
import numpy as np
import pandas as pd
from multiprocessing import Pool


d = {}
allg = set()


CHUNK = 100


def read_stream(stream):
    for l in stream:
        s = l.rstrip().split()
        if s[0] == 'pos_a':
            continue
        a, b, mi = s
        a = int(a)
        b = int(b)
        mi = float(mi)
        yield a
        yield b
        yield mi


def aracne(a, b, m):
    global d
    global allg

    bad = False
    for c in allg:
        if c == a or c == b: continue
        if (a, c) not in d: continue
        if (b, c) not in d: continue
        if (m < d[(a, c)]) and (m < d[(b, c)]):
            bad = True
            return None

    return a, b, m


def prepare_aracne(hitsA, hitsB, mis):
    global d
    global allg

    allg = set(hitsA) | set(hitsB)

    d = {}
    for a, b, m in zip(hitsA, hitsB, mis):
        d[(a,b)] = m
        d[(b,a)] = m


def iter_pos(hitsA, hitsB, mis):
    for a, b, m in tqdm.tqdm(zip(hitsA, hitsB, mis), total=mis.shape[0]):
        yield a, b, m


def get_options():
    description = 'Post-process the MI values'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('plmc',
                        help='Annotated PLMC table')
    parser.add_argument('output',
                        help='Output table')
    parser.add_argument('--cores',
                        default=1,
                        type=int,
                        help='How many cores should ARACNE use '
                             'default: %(default)d)')
    parser.add_argument('--chunk',
                        default=CHUNK,
                        type=int,
                        help='How many positions to run ARACNE for at any '
                             'given time (affects peak memory usage, '
                             'default: %(default)d)')

    return parser.parse_args()


def main():
    options = get_options()

    cores = options.cores
    chunk = options.chunk

    df = pd.read_csv(options.plmc, sep='\t')

    all_v = df[['pos_a', 'pos_b', 'plmc']].values

    pos_a = all_v[:, 0]
    pos_b = all_v[:, 1]
    mi = all_v[:, 2]

    del all_v

    prepare_aracne(pos_a, pos_b, mi)

    if cores > 1:
        pool = Pool(cores, maxtasksperchild=100)

    iter_f = iter_pos(pos_a, pos_b, mi)
    keep = set()
    if cores == 1:
        for pos in iter_f:
            ret = aracne(*pos)
            if ret is None:
                continue
            a, b, _ = ret
            a = int(a)
            b = int(b)
            keep.add((a, b))
    else:
        while True:
            ret = pool.starmap(aracne, itertools.islice(iter_f, chunk * cores))
            if len(ret) == 0:
                break
            for x in ret:
                if x is None:
                    continue
                a, b, _ = x
                a = int(a)
                b = int(b)
                keep.add((a, b))

    df = df.set_index(['pos_a', 'pos_b']).loc[sorted(keep)].reset_index()

    # remove distinction between a and b by duplicating the table and renaming
    m1 = df[['pos_a', 'pos_b', 'distance', 'outlier', 'plmc',
            'feature_position_a', 'gene_a', 'codon_a',
            'feature_codon_a', 'feature_position_b', 'gene_b',
            'codon_b', 'feature_codon_b', 'codon_distance']].rename(columns={x: x.replace('_a', '_source').replace('_b', '_target')
                                                                             for x in df.columns})
    m2 = df[['pos_b', 'pos_a', 'distance', 'outlier', 'plmc',
            'feature_position_b', 'gene_b', 'codon_b',
            'feature_codon_b', 'feature_position_a', 'gene_a',
            'codon_a', 'feature_codon_a', 'codon_distance']].rename(columns={x: x.replace('_b', '_source').replace('_a', '_target')
                                                                             for x in df.columns})
    df = pd.concat([m1, m2])
    df['interaction'] = np.where(df['gene_source'] == df['gene_target'], 'same gene', 'different gene')

    df.to_csv(options.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
