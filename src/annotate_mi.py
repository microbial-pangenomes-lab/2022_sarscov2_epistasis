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
    parser.add_argument('mi',
                        help='MI table')
    parser.add_argument('output',
                        help='Output table')
    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    r = parse_gff(options.gff)

    s = pd.read_csv(options.mi, sep=' ', header=None)
    s.columns = ['pos_a', 'pos_b', 'distance', 'outlier', 'mi']

    s1 = pd.merge(s, r, left_on='pos_a', right_on='position', how='left')
    s1 = s1.drop(columns=['position', 'strand', 'chromosome'])
    s1 = s1.rename(columns={x: f'{x}_a'
                          for x in ['ftype', 'feature_position', 'feature_codon', 'id',
                                    'gene', 'product', 'codon']})
    s1
    s2 = pd.merge(s1, r, left_on='pos_b', right_on='position', how='left')
    s2 = s2.drop(columns=['position', 'strand', 'chromosome'])
    s2 = s2.rename(columns={x: f'{x}_b'
                          for x in ['ftype', 'feature_position', 'feature_codon', 'id',
                                    'gene', 'product', 'codon']})
    # add aminoacid distance
    idx = s2[s2['gene_a'] == s2['gene_b']].index

    s2.loc[idx, 'codon_distance'] = (s2.loc[idx, 'feature_codon_a'] - s2.loc[idx, 'feature_codon_b']).abs()
    s2.loc[idx, 'codon_distance']

    # remove distinction between a and b by duplicating the table and renaming
    m1 = s2[['pos_a', 'pos_b', 'distance', 'outlier', 'mi',
            'feature_position_a', 'gene_a', 'codon_a',
            'feature_codon_a', 'feature_position_b', 'gene_b',
            'codon_b', 'feature_codon_b', 'codon_distance']].rename(columns={x: x.replace('_a', '_source').replace('_b', '_target')
                                                                             for x in s2.columns})
    m2 = s2[['pos_b', 'pos_a', 'distance', 'outlier', 'mi',
            'feature_position_b', 'gene_b', 'codon_b',
            'feature_codon_b', 'feature_position_a', 'gene_a',
            'codon_a', 'feature_codon_a', 'codon_distance']].rename(columns={x: x.replace('_b', '_source').replace('_a', '_target')
                                                                             for x in s2.columns})
    m = pd.concat([m1, m2]).reset_index()

    m['interaction'] = np.where(m['gene_source'] == m['gene_target'], 'same gene', 'different gene')

    m.drop(columns=['index']).to_csv(options.output, sep='\t', index=False)
