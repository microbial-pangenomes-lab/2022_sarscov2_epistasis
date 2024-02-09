#!/usr/bin/env python


import os
import math
import random
import warnings
import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats
from sklearn import metrics
from itertools import permutations

import escapecalculator

# RBD mutations

# based on 10.1038/s41467-022-34506-z
AFFINITY = [501, 498]
ESCAPE = [417, 446, 484, 493, 496, 505,
          # these are extra from 10.1126/science.abo7896
          447, 449, 506, 406]
# data from the escape calculator
# commit 5ebb88e
if not os.path.exists('escape.txt'):
    esc = escapecalculator.EscapeCalculator()
    esc = esc.escape_per_site([])[['site', 'original_escape']]
    ESCAPE = sorted(set(ESCAPE).union(esc.loc[esc['original_escape'] > 0.1]['site'].values))
    f = open('escape.txt', 'w')
    f.write('\n'.join([str(x) for x in esc.loc[esc['original_escape'] > 0.1]['site'].values]))
    f.close()
else:
    ESCAPE = sorted(set(ESCAPE).union([int(x.rstrip())
                                       for x in open('escape.txt')]))
# based on outbreak.info
MOI = [484, 18, 417, 439, 452, 477, 494, 501, 681]

RBD_INTERESTING = {x for x in set(MOI).union(AFFINITY).union(ESCAPE)
               if x >= 319 and x <= 540}
PAIRS = {(x, y) for x, y in permutations(RBD_INTERESTING, 2)}
ALL_PAIRS = {(x, y)
             for x, y in permutations(range(319, 541), 2)}

def read_mi(fname, name='all', date=True):
    a = pd.read_csv(fname, sep='\t')

    # filtering
    #
    # 1. end of the genome
    # 2. same gene and same codon
    # 3. different gene and too close
    #
    a = a[((a['pos_source'] <= 29674) & (a['pos_target'] <= 29674)) &
          (((a['gene_source'] == a['gene_target']) &
            (a['feature_codon_source'] != a['feature_codon_target']) &
            (a['codon_distance'] > 1)) |
           ((a['gene_source'] != a['gene_target']) &
            (a['distance'] > 1)))]

    a['name'] = name
    if date:
        a['year'] = int(name.split('-')[0])
        a['month'] = int(name.split('-')[1])

    return a


def get_rbd_mutated(indir, names=None):
    mutated = set()
    for f in os.listdir(indir):
        date = f.split('.')[0]
        if names is not None and date not in names:
            continue
        l = pd.read_csv(f'{indir}/{f}', sep='\t',
                        usecols=['seqName', 'clade', 'Nextclade_pango',
                                 'partiallyAliased', 'clade_nextstrain',
                                 'clade_who', 'clade_display',
                                 'aaSubstitutions', 'aaDeletions', 'aaInsertions'])
        year = int(date.split('-')[0])
        month = int(date.split('-')[1])
    
        d = {}
        for clade, line in l[['clade_display', 'aaSubstitutions']].values:
            if str(line) == 'nan':
                continue
            for x in line.split(','):
                k, v = x.split(':')
                if k != 'S':
                    continue
                while True:
                    v = v[1:]
                    try:
                        v = int(v)
                    except ValueError:
                        v = v[:-1]
                        try:
                            v = int(v)
                        except ValueError:
                            continue
                    break
                d[k] = d.get(k, set())
                d[k].add(v)
    
            for k, vv in d.items():
                for v in vv:
                    if v in range(319, 541):
                        mutated.add(v)
    return mutated


def rbd_network(a, mutated=None, save_gml=None, sample=1, n_random=1000):
    if mutated is None:
        mutated = range(319, 541)

    relevant = {x for x in set(ESCAPE).union(AFFINITY).union(MOI)
                if x >= 319 and x <= 540}
    
    g = nx.Graph()
    for p0, p1 in a[(a['gene_source'] == 'S') &
                    (a['gene_target'] == 'S') &
                    (a['feature_codon_source'] >= 319) &
                    (a['feature_codon_source'] <= 540) &
                    (a['feature_codon_target'] >= 319) &
                    (a['feature_codon_target'] <= 540)][[
            'feature_codon_source',
            'feature_codon_target']].sample(frac=sample).values:
        g.add_edge(int(p0), int(p1))

    for node in g.nodes:
        if node in AFFINITY:
            ntype = 'affinity'
        elif node in ESCAPE:
            ntype = 'escape'
        elif node in MOI:
            ntype = 'MOI/MOC'
        else:
            ntype = 'other'
        g.nodes[node]['type'] = ntype

    if save_gml is not None:
        nx.write_gml(g, save_gml)

    random.seed(100)

    cg = [(x, y) for x, y in g.edges]
    yield (np.array([1 if x in relevant and y in relevant
                     else 0
                     for x, y in cg]),
           'original')

    for _ in range(n_random):
        n_nodes = len(g.nodes)
        n_edges = len(g.edges)
        r = nx.Graph()
        nodes = list({x for x in range(319, 541) if x in mutated})
        if len(nodes) < len(g.nodes):
            nodes = list(g.nodes)
        while len(r.edges) < n_edges:
            n1 = random.choice(nodes)
            n2 = random.choice(nodes)
            if n1 != n2 and abs(n1 - n2) > 1:
                r.add_edge(n1, n2)
        cg = [(x, y) for x, y in r.edges]
        yield (np.array([1 if x in relevant and y in relevant
                         else 0
                         for x, y in cg]),
               'random')


def enrichment_score(*g):
    g = g[0]
    
    RBD_LENGTH = 540 - 319 + 1
    relevant = {x for x in set(ESCAPE).union(AFFINITY).union(MOI)
                if x >= 319 and x <= 540}
    
    all_possible = math.factorial(RBD_LENGTH) / (math.factorial(2) * math.factorial(RBD_LENGTH - 2))
    int_rel = g.sum()
    int_not = g.shape[0] - int_rel
    all_possible_rel = math.factorial(len(relevant)) / (math.factorial(2) * math.factorial(len(relevant) - 2))
    not_int_rel = all_possible_rel - int_rel
    not_int_not_rel = all_possible - not_int_rel - int_not - int_rel

    table = [[int_rel, not_int_rel],
             [int_not, not_int_not_rel]]
    odds_ratio, pvalue = stats.fisher_exact(table,
                                            alternative='greater')

    return odds_ratio


def enrichment(a, mutated=None, save_gml=None, sample=1, n_random=1000):
    res = []
    for g, ntype in rbd_network(a, save_gml=save_gml,
                                mutated=mutated, sample=sample,
                                n_random=n_random):
        odds_ratio = enrichment_score(g)
        if ntype == 'original':
            ci = stats.bootstrap((g,), enrichment_score, n_resamples=999,
                                 confidence_level=0.95)
            high = ci.confidence_interval.high
            low = ci.confidence_interval.low
        else:
            high = np.nan
            low = np.nan

        res.append((ntype, odds_ratio, low, high))

    df = pd.DataFrame(res,
                      columns=['type', 'odds-ratio', 'low', 'high'])

    return df


def get_ml_data(a, mutated=None, shuffle=False, col='mi'):
    if mutated is None:
        mutated = range(319, 541)
    else:
        mutated = sorted(mutated)

    tl = a.groupby('outlier')[col].min().to_dict()

    a = a[(a['gene_source'] == 'S') &
          (a['gene_target'] == 'S') &
          (a['feature_codon_source'].isin(range(319, 541))) &
          (a['feature_codon_target'].isin(range(319, 541)))]

    if shuffle:
        a = a.copy()
        values = []
        observed = sorted(set(a['feature_codon_source']).union(a['feature_codon_target']))
        if len(observed) > len(mutated):
            mutated = observed
        i = 0
        while len(values) < a.shape[0]:
            p1 = random.choice(mutated)
            p2 = random.choice(mutated)
            if p1 != p2 and abs(p1 - p2) > 1 and (p1, p2) not in values:
                values.append((p1, p2))
            i += 1
            if i > 100 * len(mutated):
                # give up
                return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        a['feature_codon_source'] = [x for x, _ in values]
        a['feature_codon_target'] = [x for _, x in values]
    
    da = a.set_index(['feature_codon_source', 'feature_codon_target'])[col].to_dict()
    rda = {}
    for (p1, p2), v in da.items():
        rda[(int(p1), int(p2))] = v
    
    y_true = np.array([1 if k in PAIRS
                       else 0
                       for k in ALL_PAIRS
                       if k in rda])
    
    y_values = np.array([rda.get(k, np.nan) for k in ALL_PAIRS if k in rda])

    v1, v2, v3, v4 = np.nan, np.nan, np.nan, np.nan
    for thr, var in zip(range(1, 5),
                        [v1, v2, v3, v4]):
        try:
            var = y_values > tl[thr]
        except KeyError:
            var = np.nan

    return y_true, y_values, v1, v2, v3, v4


def get_specificity(y_true, y_values):
    try:
        tn, fp, fn, tp = metrics.confusion_matrix(y_true, y_values).ravel()
    except ValueError:
        return np.nan
    return tn / (tn + fp)


def get_sensitivity(y_true, y_values):
    try:
        tn, fp, fn, tp = metrics.confusion_matrix(y_true, y_values).ravel()
    except ValueError:
        return np.nan
    return tp / (tp + fn)


def get_roc_auc(y_true, y_values):
    x, y, _ = metrics.roc_curve(y_true, y_values)
    return metrics.auc(x, y)


def get_pr_auc(y_true, y_values):
    y, x, _ = metrics.precision_recall_curve(y_true, y_values)
    return metrics.auc(x, y)


def ml_metrics(a, mutated=None, shuffle=False, col='mi'):
    y_true, y_values, v1, v2, v3, v4 = get_ml_data(a,
                                                   mutated=mutated,
                                                   shuffle=shuffle,
                                                   col=col)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if np.isnan(v1):
            f1 = np.nan
        else:
            f1 = metrics.f1_score(y_true, v1)
        if np.isnan(v2):
            f2 = np.nan
        else:
            f2 = metrics.f1_score(y_true, v2)
        if np.isnan(v3):
            f3 = np.nan
        else:
            f3 = metrics.f1_score(y_true, v3)
        if np.isnan(v4):
            f4 = np.nan
        else:
            f4 = metrics.f1_score(y_true, v4)
    
        if not shuffle:
            try:
                ci = stats.bootstrap((y_true, v1), metrics.f1_score, paired=True,
                                     n_resamples=999)
                f1l, f1h = ci.confidence_interval.low, ci.confidence_interval.high
            except ValueError:
                f1l, f1h = np.nan, np.nan
            try:
                ci = stats.bootstrap((y_true, v2), metrics.f1_score, paired=True,
                                     n_resamples=999)
                f2l, f2h = ci.confidence_interval.low, ci.confidence_interval.high
            except ValueError:
                f2l, f2h = np.nan, np.nan
            try:
                ci = stats.bootstrap((y_true, v3), metrics.f1_score, paired=True,
                                     n_resamples=999)
                f3l, f3h = ci.confidence_interval.low, ci.confidence_interval.high
            except ValueError:
                f3l, f3h = np.nan, np.nan
            try:
                ci = stats.bootstrap((y_true, v4), metrics.f1_score, paired=True,
                                     n_resamples=999)
                f4l, f4h = ci.confidence_interval.low, ci.confidence_interval.high
            except ValueError:
                f4l, f4h = np.nan, np.nan
        else:
            f1l, f1h = np.nan, np.nan
            f2l, f2h = np.nan, np.nan
            f3l, f3h = np.nan, np.nan
            f4l, f4h = np.nan, np.nan

        spec1 = get_specificity(y_true, v1)
        sens1 = get_sensitivity(y_true, v1)
        spec2 = get_specificity(y_true, v2)
        sens2 = get_sensitivity(y_true, v2)
        spec3 = get_specificity(y_true, v3)
        sens3 = get_sensitivity(y_true, v3)
        spec4 = get_specificity(y_true, v4)
        sens4 = get_sensitivity(y_true, v4)
        if not shuffle:
            try:
                ci = stats.bootstrap((y_true, v1), get_specificity, paired=True,
                                     n_resamples=999)
                sp1l, sp1h = ci.confidence_interval.low, ci.confidence_interval.high
            except ValueError:
                sp1l, sp1h = np.nan, np.nan
            try:
                ci = stats.bootstrap((y_true, v2), get_specificity, paired=True,
                                     n_resamples=999)
                sp2l, sp2h = ci.confidence_interval.low, ci.confidence_interval.high
            except ValueError:
                sp2l, sp2h = np.nan, np.nan
            try:
                ci = stats.bootstrap((y_true, v3), get_specificity, paired=True,
                                     n_resamples=999)
                sp3l, sp3h = ci.confidence_interval.low, ci.confidence_interval.high
            except ValueError:
                sp3l, sp3h = np.nan, np.nan
            try:
                ci = stats.bootstrap((y_true, v4), get_specificity, paired=True,
                                     n_resamples=999)
                sp4l, sp4h = ci.confidence_interval.low, ci.confidence_interval.high
            except ValueError:
                sp4l, sp4h = np.nan, np.nan

            try:
                ci = stats.bootstrap((y_true, v1), get_sensitivity, paired=True,
                                     n_resamples=999)
                se1l, se1h = ci.confidence_interval.low, ci.confidence_interval.high
            except ValueError:
                se1l, se1h = np.nan, np.nan
            try:
                ci = stats.bootstrap((y_true, v2), get_sensitivity, paired=True,
                                     n_resamples=999)
                se2l, se2h = ci.confidence_interval.low, ci.confidence_interval.high
            except ValueError:
                se2l, se2h = np.nan, np.nan
            try:
                ci = stats.bootstrap((y_true, v3), get_sensitivity, paired=True,
                                     n_resamples=999)
                se3l, se3h = ci.confidence_interval.low, ci.confidence_interval.high
            except ValueError:
                se3l, se3h = np.nan, np.nan
            try:
                ci = stats.bootstrap((y_true, v4), get_sensitivity, paired=True,
                                     n_resamples=999)
                se4l, se4h = ci.confidence_interval.low, ci.confidence_interval.high
            except ValueError:
                se4l, se4h = np.nan, np.nan
        else:
            sp1l, sp1h = np.nan, np.nan
            sp2l, sp2h = np.nan, np.nan
            sp3l, sp3h = np.nan, np.nan
            sp4l, sp4h = np.nan, np.nan

            se1l, se1h = np.nan, np.nan
            se2l, se2h = np.nan, np.nan
            se3l, se3h = np.nan, np.nan
            se4l, se4h = np.nan, np.nan

        # roc_auc = get_roc_auc(y_true, y_values)
        # pr_auc = get_pr_auc(y_true, y_values)
    
        # if not shuffle:
            # ci = stats.bootstrap((y_true, y_values), get_roc_auc, paired=True, n_resamples=999)
        #     rl, rh = ci.confidence_interval.low, ci.confidence_interval.high
        #     ci = stats.bootstrap((y_true, y_values), get_pr_auc, paired=True, n_resamples=999)
        #     pl, ph = ci.confidence_interval.low, ci.confidence_interval.high
        # else:
        #     rl, rh = np.nan, np.nan
        #     pl, ph = np.nan, np.nan
    
    return pd.DataFrame((
                         ('f1', 1, f1, f1l, f1h),
                         ('f1', 2, f2, f2l, f2h),
                         ('f1', 3, f3, f3l, f3h),
                         ('f1', 4, f4, f4l, f4h),
                         ('specificity', 1, spec1, sp1l, sp1h),
                         ('specificity', 2, spec2, sp2l, sp2h),
                         ('specificity', 3, spec3, sp3l, sp4h),
                         ('specificity', 4, spec4, sp3l, sp4h),
                         ('sensitivity', 1, sens1, se1l, se1h),
                         ('sensitivity', 2, sens2, se2l, se2h),
                         ('sensitivity', 3, sens3, se3l, se4h),
                         ('sensitivity', 4, sens4, se3l, se4h),
                         # ('pr_auc', np.nan, pr_auc, pl, ph),
                         # ('roc_auc', np.nan, roc_auc, rl, rh)
                         ),
                        columns=['metric', 'outlier',
                                 'value', 'low', 'high']
                       )


if __name__ == "__main__":
    pass
