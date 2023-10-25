#!/usr/bin/env python


import gzip
import argparse
import numpy as np
import pandas as pd

from collections import namedtuple

Feature = namedtuple('Feature', ['id',
                                 'ftype',
                                 'chromosome',
                                 'start',
                                 'end',
                                 'strand',
                                 'gene',
                                 'product'])


def parse_gff(file_name):
    # output dict
    # key: feature ID
    # value: Feature NamedTuple
    features = {}

    with gzip.open(file_name, 'rt') as gff:
        for line in gff:
            if line.lstrip().startswith('##FASTA'):
                # start of FASTA entries, end of file
                break

            elif line.lstrip().startswith('#'):
                # comment, ignore
                continue

            # should be a valid GFF3 line
            entries = line.split('\t')

            try:
                ftype = entries[2]

                chrom = entries[0]
                start = int(entries[3])
                end = int(entries[4])
                strand = entries[6]

                # integer takes up less space
                if strand == '+':
                    strand = 1
                else:
                    strand = -1

                # fetch the feature ID from the last field
                ID = None
                for entry in entries[8].split(';'):
                    if entry.startswith('ID=') and '=' in entry:
                        ID = entry.split('=')[1]

                # could not find it, skip this entry
                if ID is None:
                    continue

                # fetch the gene name
                gene = np.nan
                for entry in entries[8].split(';'):
                    if entry.startswith('gene=') and '=' in entry:
                        gene = entry.split('=')[1]

                product = np.nan
                for entry in entries[8].split(';'):
                    if entry.startswith('product=') and '=' in entry:
                        product = entry.split('=')[1]

                # save the relevant details
                features[ID] = Feature(ID, ftype, chrom, start, end, strand, gene, product)

            except Exception as e:
                # not distinguishing between exceptions
                # not great behaviour
                logger.warning(f'{e}, skipping line "{line.rstrip()}" from {file_name}')
                continue

    res = []
    for feat in features.values():
        for i, position in enumerate(range(feat.start, feat.end+1)):
            res.append( (feat.chromosome, feat.ftype, position, i+1, feat.strand, feat.id, feat.gene, feat.product) )
    r = pd.DataFrame(res,
                     columns=['chromosome', 'ftype', 'position', 'feature_position', 'strand', 'id', 'gene', 'product'])

    r['codon'] = np.nan
    r['codon'] = (r['feature_position'] % 3).astype(int)
    r.loc[r[r['codon'] == 0].index, 'codon'] = 3
    r.loc[r[r['ftype'] != 'CDS'].index, 'codon'] = np.nan

    r['feature_codon'] = (r['feature_position'] / 3).astype(int)+1
    r.loc[r[r['ftype'] != 'CDS'].index, 'feature_codon'] = np.nan
    r.loc[r[r['codon'] == 3.0].index, 'feature_codon'] -= 1

    r = r[r['ftype'].isin(['five_prime_UTR', 'three_prime_UTR', 'CDS'])]

    return r


def get_options():
    description = 'Annotate MI table using a GFF file'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('gff',
                        help='GFF file')
    parser.add_argument('plmc',
                        help='PLMC table (with extension .EC)')
    parser.add_argument('output',
                        help='Output table')
    parser.add_argument('--keep-all',
                        action='store_true',
                        default=False,
                        help='Keep all values, not just outliers')
    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    r = parse_gff(options.gff)

    s = pd.read_csv(options.plmc, sep=' ', header=None)
    s.columns = ['pos_a', 'base_a',
                 'pos_b', 'base_b',
                 '_', 'plmc']
    s['pos_a'] += 21562
    s['pos_b'] += 21562
    s['distance'] = abs(s['pos_b'] - s['pos_a'])

    Q1, Q3 = np.quantile(s['plmc'], [0.25, 0.75])
    t1 = Q3 + 6 * (Q3 - Q1)
    t2 = Q3 + 9 * (Q3 - Q1)
    t3 = Q3 + 13.5 * (Q3 - Q1)
    t4 = Q3 + 20.25 * (Q3 - Q1)

    s.loc[s[s['plmc'] < t1].index, 'outlier'] = 0
    s.loc[s[s['plmc'] >= t1].index, 'outlier'] = 1
    s.loc[s[s['plmc'] >= t2].index, 'outlier'] = 2
    s.loc[s[s['plmc'] >= t3].index, 'outlier'] = 3
    s.loc[s[s['plmc'] >= t4].index, 'outlier'] = 4

    if not options.keep_all:
        s = s[s['outlier'] > 0]
    s['outlier'] = s['outlier'].astype(np.int64)

    s1 = pd.merge(s, r, left_on='pos_a', right_on='position', how='left')
    s1 = s1.rename(columns={x: f'{x}_a'
                          for x in ['ftype', 'feature_position', 'feature_codon', 'id',
                                    'gene', 'product', 'codon']})
    s1 = s1.drop(columns=['position', 'strand', 'chromosome'])
    s2 = pd.merge(s1, r, left_on='pos_b', right_on='position', how='left')
    s2 = s2.rename(columns={x: f'{x}_b'
                          for x in ['ftype', 'feature_position', 'feature_codon', 'id',
                                    'gene', 'product', 'codon']})
    s2 = s2.drop(columns=['position', 'strand', 'chromosome'])
    # add aminoacid distance
    idx = s2[s2['gene_a'] == s2['gene_b']].index

    s2.loc[idx, 'codon_distance'] = (s2.loc[idx, 'feature_codon_a'] - s2.loc[idx, 'feature_codon_b']).abs()
    s2.loc[idx, 'codon_distance']

    # remove distinction between a and b by duplicating the table and renaming
    m1 = s2[['pos_a', 'pos_b', 'distance', 'outlier', 'plmc',
            'feature_position_a', 'gene_a', 'codon_a',
            'feature_codon_a', 'feature_position_b', 'gene_b',
            'codon_b', 'feature_codon_b', 'codon_distance']].rename(columns={x: x.replace('_a', '_source').replace('_b', '_target')
                                                                             for x in s2.columns})
    m2 = s2[['pos_b', 'pos_a', 'distance', 'outlier', 'plmc',
            'feature_position_b', 'gene_b', 'codon_b',
            'feature_codon_b', 'feature_position_a', 'gene_a',
            'codon_a', 'feature_codon_a', 'codon_distance']].rename(columns={x: x.replace('_b', '_source').replace('_a', '_target')
                                                                             for x in s2.columns})
    m = pd.concat([m1, m2]).reset_index()
    m['interaction'] = np.where(m['gene_source'] == m['gene_target'], 'same gene', 'different gene')

    m.to_csv(options.output, sep='\t', index=False)
