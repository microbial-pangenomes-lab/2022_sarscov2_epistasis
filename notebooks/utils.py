#!/usr/bin/env python


import os
import math
import random
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
esc = escapecalculator.EscapeCalculator()
esc = esc.escape_per_site([])[['site', 'original_escape']]
ESCAPE = sorted(set(ESCAPE).union(esc.loc[esc['original_escape'] > 0.1]['site'].values))
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


def enrichment(a, mutated=None, save_gml=None):
    if mutated is None:
        mutated = range(319, 541)
    
    g = nx.Graph()
    for p0, p1 in a[(a['gene_source'] == 'S') &
                    (a['gene_target'] == 'S') &
                    (a['feature_codon_source'] >= 319) &
                    (a['feature_codon_source'] <= 540) &
                    (a['feature_codon_target'] >= 319) &
                    (a['feature_codon_target'] <= 540)][['feature_codon_source',
                                                         'feature_codon_target']].values:
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

    # Simulate if enriched compared to random
    relevant = {x for x in set(ESCAPE).union(AFFINITY).union(MOI)
                if x >= 319 and x <= 540}

    RBD_LENGTH = 540 - 319 + 1

    res = []

    # randomization #1: any position in the RBD
    # randomization #2: same positions as in the original graph
    rgraphs1 = []
    rgraphs2 = []
    for _ in range(1000):
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
        rgraphs1.append(r)
        node_mapping = dict(zip(g.nodes(),
                            sorted(g.nodes(), key=lambda k: random.random())))
        r = nx.relabel_nodes(g, node_mapping)
        # remove connections between subsequent nodes
        # change them to another random edge
        adj = [(x, y) for x, y in r.edges
               if abs(x - y) == 1]
        if len(adj) > 1:
            orig = len(r.edges)
            nodes = list(r.nodes)
            r.remove_edges_from(adj)
            while len(r.edges) < orig:
                n1 = random.choice(nodes)
                n2 = random.choice(nodes)
                if n1 != n2 and abs(n1 - n2) > 1:
                    r.add_edge(n1, n2)
        #
        rgraphs2.append(r)

    for i, cg in enumerate([g] + rgraphs1):
        all_possible = math.factorial(RBD_LENGTH) / (math.factorial(2) * math.factorial(RBD_LENGTH - 2))
        int_rel = len([(x, y) for x, y in cg.edges if x in relevant and y in relevant])
        int_not = len(cg.edges) - int_rel
        all_possible_rel = math.factorial(len(relevant)) / (math.factorial(2) * math.factorial(len(relevant) - 2))
        not_int_rel = all_possible_rel - int_rel
        not_int_not_rel = all_possible - not_int_rel - int_not - int_rel

        table = [[int_rel, not_int_rel],
                 [int_not, not_int_not_rel]]
        odds_ratio, pvalue = stats.fisher_exact(table,
                                                alternative='greater')

        if i == 0:
            gtype = 'original'
            i = np.nan
        else:
            gtype = 'random'

        res.append((gtype, i, odds_ratio, pvalue, 'any_positions'))

    # for i, cg in enumerate([g] + rgraphs2):
    #     all_possible = math.factorial(len(cg.nodes)) / (math.factorial(2) * math.factorial(len(cg.nodes) - 2))
    #     int_rel = len([(x, y) for x, y in cg.edges if x in relevant and y in relevant])
    #     int_not = len(cg.edges) - int_rel
    #     all_possible_rel = math.factorial(len(relevant)) / (math.factorial(2) * math.factorial(len(relevant) - 2))
    #     not_int_rel = all_possible_rel - int_rel
    #     not_int_not_rel = all_possible - not_int_rel - int_not - int_rel

    #     table = [[int_rel, not_int_rel],
    #              [int_not, not_int_not_rel]]
    #     print(table, all_possible)
    #     odds_ratio, pvalue = stats.fisher_exact(table,
    #                                             alternative='greater')

    #     if i == 0:
    #         gtype = 'original'
    #         i = np.nan
    #     else:
    #         gtype = 'random'

    #     res.append((gtype, i, odds_ratio, pvalue, 'fixed_positions'))

    df = pd.DataFrame(res,
                      columns=['type', 'round', 'odds-ratio', 'p-value', 'randomization'])

    return df


def ml_metrics(a, mutated=None, shuffle=False):
    if mutated is None:
        mutated = range(319, 541)
    else:
        mutated = sorted(mutated)

    tl = a.groupby('outlier')['mi'].min().to_dict()

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
    
    da = a.set_index(['feature_codon_source', 'feature_codon_target'])['mi'].to_dict()
    rda = {}
    for (p1, p2), v in da.items():
        rda[(int(p1), int(p2))] = v
    
    y_true = np.array([1 if k in PAIRS
                       else 0
                       for k in ALL_PAIRS
                       if k in rda])
    
    y_values = np.array([rda.get(k, np.nan) for k in ALL_PAIRS if k in rda])

    f1 = metrics.f1_score(y_true, y_values > tl[1])
    f2 = metrics.f1_score(y_true, y_values > tl[2])
    f3 = metrics.f1_score(y_true, y_values > tl[3])
    f4 = metrics.f1_score(y_true, y_values > tl[4])
    
    pr_auc = metrics.auc(metrics.precision_recall_curve(y_true, y_values)[1],
                         metrics.precision_recall_curve(y_true, y_values)[0])
    
    roc_auc = metrics.auc(metrics.roc_curve(y_true, y_values)[0],
                          metrics.roc_curve(y_true, y_values)[1])

    return f1, f2, f3, f4, pr_auc, roc_auc


if __name__ == "__main__":
    pass
