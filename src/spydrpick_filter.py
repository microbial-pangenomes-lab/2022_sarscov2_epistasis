#!/usr/bin/env python


import sys
import tqdm
import argparse
import itertools
import numpy as np
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

    parser.add_argument('--outliers',
                        default=None,
                        help='Outliers file (one value per line'
                             'default: do not report outliers)')
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

    outliers = []
    if options.outliers is not None:
        tmp = []
        for v in open(options.outliers):
            tmp.append(float(v.rstrip()))
        outliers = [(i+1, v) for i, v in enumerate(sorted(tmp))]

    all_v = np.fromiter(read_stream(sys.stdin), dtype=float)
    all_v.shape = int(all_v.shape[0] / 3), 3

    pos_a = all_v[:, 0]
    pos_b = all_v[:, 1]
    mi = all_v[:, 2]

    del all_v

    prepare_aracne(pos_a, pos_b, mi)

    if cores > 1:
        pool = Pool(cores, maxtasksperchild=100)

    iter_f = iter_pos(pos_a, pos_b, mi)
    print('pos_a\tpos_b\tmi\toutlier')
    if cores == 1:
        for pos in iter_f:
            ret = aracne(*pos)
            if ret is None:
                continue
            a, b, m = ret
            a = int(a)
            b = int(b)
            outl = 0
            for i, v in outliers:
                if m >= v:
                    outl = i
                else:
                    break

            print(f'{a}\t{b}\t{m}\t{outl}')
    else:
        while True:
            ret = pool.starmap(aracne, itertools.islice(iter_f, chunk * cores))
            if len(ret) == 0:
                break
            for x in ret:
                if x is None:
                    continue
                a, b, m = x
                a = int(a)
                b = int(b)
                outl = 0
                for i, v in outliers:
                    if m >= v:
                        outl = i
                    else:
                        break

                print(f'{a}\t{b}\t{m}\t{outl}')


if __name__ == "__main__":
    main()
