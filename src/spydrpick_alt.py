#!/usr/bin/env python


import os
import sys
import argparse
import itertools

CHUNK = 1000
SAMPLES = 100
WINDOW = 1000


def read_positions(pos_filename):
    sys.stderr.write(f'Load: Positions\n')
    pos = np.loadtxt(pos_filename, dtype=np.uint,
                     delimiter=",", comments=None)

    return pos


def read_weights(weights_file):
    sys.stderr.write('Load: Weights\n')
    w = np.loadtxt(weights_file,
                   delimiter=' ', comments=None)

    return w


def read_binary(filename):
    sys.stderr.write(f'Load: Aligment binary\n')
    b = np.load(filename)['arr_0'].T

    return b


def joint_probabilities(j, m, x, y, weights, pos_n, n_eff):
        a = m[pos_n]
        a = (a == x).astype(np.float64)
        asum = np.empty(a.shape[0])
        for i in range(a.shape[0]):
            asum[i] = np.unique(a[i]).shape[0]
        a = weights * a
        for p_m in range(0, m.shape[0], CHUNK):
            b = m[p_m: p_m + CHUNK]
            b = (b == y).astype(np.float64)
            bsum = np.empty(b.shape[0])
            for i in range(b.shape[0]):
                bsum[i] = np.unique(b[i]).shape[0]
            r_x_r_y = (asum * bsum.reshape(-1, 1) * 0.5).T
            d = (np.inner(a, b) + 0.5) / (n_eff + r_x_r_y)
            j[:, p_m: p_m + CHUNK] = d


def compute_mi(m, pos_n, weights, nstates=5):
    n_eff = np.sum(weights)
    sys.stderr.write(f'MI: effective sample size: {n_eff}\n')

    # calculate co-occurence matrix using matrix algebra
    pj = {}
    for x in range(nstates):
        pj[x] = {}
        for y in range(nstates):
            sys.stderr.write(f'MI: compute joint probabilities ({x}, {y})\n')
            pj[x][y] = np.empty((pos_n.shape[0], m.shape[0]))
            joint_probabilities(pj[x][y], m, x, y, weights, pos_n, n_eff)

    sys.stderr.write(f'MI: compute marginal probabilities\n')
    ph = {}
    pi = {}
    for x in range(nstates):
        ph[x] = 0
        pi[x] = 0
        for y in range(nstates):
            ph[x] += pj[x][y]
            pi[x] += pj[y][x]

    sys.stderr.write(f'MI: compute mutual information\n')
    mi = np.empty((pos_n.shape[0], m.shape[0]))
    for x in range(nstates):
        for y in range(nstates):
            mi += pj[x][y] * (np.log(pj[x][y]) - np.log(ph[x]) - np.log(pi[y]))

    return mi


def get_options():
    description = 'Run spydrpick algorithm on an alignment'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-a',
                        '--alignment',
                        required=True,
                        help='Alignment file (uint8 npz format)')
    parser.add_argument('-p',
                        '--positions',
                        required=True,
                        help='Positions file (single line, comma-separated)')

    parser.add_argument('-w',
                        '--weights',
                        default=None,
                        help='Weights file (one value per line)')
    parser.add_argument('-n',
                        '--states',
                        default=5,
                        type=int,
                        help='States (default: %(default)d)')
    parser.add_argument('-q',
                        '--quantile',
                        default=0.9,
                        type=float,
                        help='Quantile to keep for MI threshold '
                             '(default: %(default).2f)')
    parser.add_argument('-t',
                        '--threshold',
                        default=None,
                        type=float,
                        help='MI threshold '
                             '(default: compute it)')
    parser.add_argument('--start',
                        default=0,
                        type=int,
                        help='Alignment column to start from (default: %(default)d)')
    parser.add_argument('--window',
                        default=WINDOW,
                        type=int,
                        help='Alignment region to use (default: start + %(default)d)')
    parser.add_argument('--chunk',
                        default=CHUNK,
                        type=int,
                        help='Window size while computing MI (affects peak memory, '
                             'default: %(default)d)')
    parser.add_argument('--cores',
                        default=1,
                        type=int,
                        help='Cores to use (default: %(default)d)')

    return parser.parse_args()


def main(options):
    pos = read_positions(options.positions)
    aln = read_binary(options.alignment)

    positions, samples = aln.shape

    if options.weights is not None:
        weights = read_weights(options.weights)
    else:
        weights = np.ones(samples)

    if options.threshold is None:
        r_positions = min(SAMPLES, positions)
        sys.stderr.write(f'Computing MI threshold on {r_positions} sample positions\n')
        sample_rows = np.random.randint(0, positions,
                                        size=r_positions)

        mi = compute_mi(aln, sample_rows, weights,
                        nstates=options.states)

        mi_t = np.quantile(mi, options.quantile)
        print(mi_t)
    else:
        if options.start > aln.shape[0]:
            raise ValueError(f'Start position {options.start} is above '
                             'the alignment size')
        if options.start + options.window > aln.shape[0]:
            rows = np.arange(options.start, aln.shape[0])
        else:
            rows = np.arange(options.start, options.start + options.window)
        mi = compute_mi(aln, rows, weights,
                        nstates=options.states)

        np.fill_diagonal(mi, options.threshold - 1)

        hitA, hitB = np.where(mi > options.threshold)
        mi = mi[hitA, hitB]
        hitA += options.start
        keep = hitA!=hitB

        hits = set()
        for a, b, v in zip(hitA[keep], hitB[keep], mi[keep]):
            a = pos[a]
            b = pos[b]
            v = round(v, 5)
            if a <= b:
                hits.add((a, b, v))
            else:
                hits.add((b, a, v))

        print(f'pos_a\tpos_b\tmi')
        for a, b, v in sorted(hits):
            print(f'{a}\t{b}\t{v}')


if __name__ == "__main__":
    options = get_options()

    CHUNK = options.chunk

    cores = str(options.cores)
    # force numpy to use the number of cores we want
    os.environ['OPENBLAS_NUM_THREADS'] = cores
    os.environ['MKL_NUM_THREADS'] = cores
    os.environ['NUMEXPR_NUM_THREADS'] = cores
    os.environ['OMP_NUM_THREADS'] = cores
    #
    import numpy as np

    main(options)
