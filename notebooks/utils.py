#!/usr/bin/env python


import math
import random
import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats

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


def enrichment(a, save_gml=None):
    g = nx.Graph()
    for p0, p1 in a[(a['gene_source'] == 'S') &
                    (a['gene_target'] == 'S') &
                       (a['pos_source'] > 22519) &
                       (a['pos_source'] < 23186) &
                       (a['pos_target'] > 22519) &
                       (a['pos_target'] < 23186)][['feature_codon_source',
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
    relevant = set(ESCAPE).union(AFFINITY).union(MOI)

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
        nodes = list(range(319, 541))
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


if __name__ == "__main__":
    pass
