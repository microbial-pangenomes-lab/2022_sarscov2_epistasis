#!/usr/bin/python3
import gzip
import os
import pandas as pd
import numpy as np
import sys

os.mkdir('results')

from collections import namedtuple

Feature = namedtuple('Feature', ['id',
                                 'ftype',
                                 'chromosome',
                                 'start',
                                 'end',
                                 'strand',
                                 'gene',
                                 'product'])

# Take a gff and extract relevant columns and make a feature table
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
                sys.stderr.write(f'{e}, skipping line "{line.rstrip()}" from {file_name}')
                continue

    return features

# Add codon feature
def add_codon(r):
    r['codon'] = np.nan
    r['codon'] = (r['feature_position'] % 3).astype(int)
    r.loc[r[r['codon'] == 0].index, 'codon'] = 3
    r.loc[r[r['ftype'] != 'CDS'].index, 'codon'] = np.nan
    r['feature_codon'] = (r['feature_position'] / 3).astype(int) + 1
    r.loc[r[r['ftype'] != 'CDS'].index, 'feature_codon'] = np.nan
    # r.loc[r['codon'] == 3.0]
    r = r[r['ftype'].isin(['five_prime_UTR', 'three_prime_UTR', 'CDS'])]
    return r


features = parse_gff('data/GCF_009858895.2_ASM985889v3_genomic.gff.gz')

res = []
for feat in features.values():
    for i, position in enumerate(range(feat.start, feat.end+1)):
        res.append((feat.chromosome, feat.ftype, position, i+1, feat.strand, feat.id, feat.gene, feat.product))

r = pd.DataFrame(res,
                 columns=['chromosome', 'ftype', 'position', 'feature_position', 'strand', 'id', 'gene', 'product'])

r = add_codon(r)

# Read spydrpick output
s = pd.read_csv('out_spydrpick.txt', sep=' ', header=None)
s.columns = ['pos_a', 'pos_b', 'distance', 'outlier', 'mi', '_']
s = s.drop(columns='_')

# Merge spydrpick output with the gff feature table
s = pd.merge(s, r, left_on='pos_a', right_on='position', how='left')
s = s.drop(columns=['position', 'strand', 'chromosome'])
s = s.rename(columns={x: f'{x}_a'
                      for x in ['ftype', 'feature_position', 'feature_codon', 'id',
                                'gene', 'product', 'codon']})
s = pd.merge(s, r, left_on='pos_b', right_on='position', how='left')
s = s.drop(columns=['position', 'strand', 'chromosome'])
s = s.rename(columns={x: f'{x}_b'
                      for x in ['ftype', 'feature_position', 'feature_codon', 'id',
                                'gene', 'product', 'codon']})
# Save to results folder
s.to_csv('results/annotated_mi_results.tsv', sep='\t', index=False)


# Load NCBI SNPs dataset  at  https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/scov2_snp to label data
mutations = pd.read_csv('data/mutations_public_ncbi.tsv', sep='\t', header=None)
mutations.columns = ['Protein', 'Protein Change', 'Count', 'Genomic Location', 'Codon Change', 'Protein Change Type', 'Collection Location']

# Change genomic location due to one error in NCBI table: 1 position behind
gen_loc=list(mutations['Genomic Location'])
gen_loc_new = []
for x in gen_loc:
    gen_loc_new.append(int(x) -1)

# Merge our data with the parsed NCBI data_table
merged = pd.merge(s, mutations, left_on='pos_a', right_on='Genomic Location', how='left')
merged2 = pd.merge(s, mutations, left_on='pos_b', right_on='Genomic Location', how='left')
total = pd.merge(merged, merged2, how='outer')


# select only mutations which are flagged as non-synonymous in the mutations file, that are mi outliers (1 or 2), and are not in the same codon. Order them firstly for 'mi' value, then for 'Count' value
relevant_mi_table = total[(total['Protein Change Type'] == 'non_synonymous') & (total['outlier'] > 0) & ((total['gene_a'] == total ['gene_b']) & (total['feature_codon_a'] != total['feature_codon_b']))].sort_values(by=['mi','Count'], ascending=False)

# And mutations coevolving in different genes
relevant_intergenic = total[(total['Protein Change Type'] == 'non_synonymous') & (total['outlier'] >0) & (total['gene_a'] != total['gene_b'])]
# And save it in results
relevant_intergenic.to_csv('results/relevant_intergenic.csv', sep='\t',header=True)

# Add outbreak.info labels
# Load outbreak.info table
outbreak_snps = pd.read_csv('data/outbreak_mutations.txt',sep='\t', header=None)
outbreak_snps.columns = ['mutation', 'type of mutation']
# Consider in this step only Spike gene
relevant_S= relevant_mi_table[(relevant_mi_table['gene_a'] == 'S') |(relevant_mi_table['gene_b'] == 'S') ]
relevant_S.to_csv('results/relevant_s.csv', sep='\t')
# Use outbreak snps dataset to label relevant_S table
outbreak_table = pd.merge(relevant_S, outbreak_snps, left_on='Protein Change', right_on='mutation', how='left')
outbreak_table.to_csv('results/relevant_mutations_outbreak_merged.tsv', sep='\t')
# filter out positions without outbreak labels
outbreak_filtered = outbreak_table[outbreak_table['mutation'].notnull() == True].reset_index()
outbreak_filtered.to_csv('results/selected_outbreak_mutations_only.tsv', sep='\t')

###### PLOTS ######
# Manhattan plot on position a (only above threshold)
import matplotlib as plt
import seaborn as sns
# Using seaborn pos_a plot
fig, ax = plt.subplots(1,1,figsize=(20, 10))
sns.scatterplot(data=relevant_mi_table,
                x = 'pos_a',
                y='mi',
                hue='gene_a',
                hue_order= ['ORF1ab', 'S', 'ORF3a', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF10'],
                size = 'gene_a',
                sizes=(20,20),
                linewidth=0)
# Add interval lines for RBD (Receptor binding domain)
plt.axvline(x=23186, color='r', ls='--', label='RBD')
plt.axvline(x=22520, ls='--', color='r')
p2 = plt.legend(loc = 'lower right', bbox_to_anchor=(1.1, 0.25), ncol=1);
# Save in results
plt.savefig('results/relevant_snps_a.pdf')

# GWES plot
s['pos_a'] = pd.to_numeric(s['pos_a'])
s['pos_b'] = pd.to_numeric(s['pos_b'])
plt.figure(figsize=(15,10))
sns.scatterplot(data=s, x='pos_a', y='mi', linewidth=0, alpha = 0.7, s=10, palette="deep") # if different colors for
                                                                                            # each outlier are desired:
                                                                                            # add hue="outlier"
plt.axhline(0.015875, c='red',linestyle='dashed', linewidth=1);
plt.axhline(0.02504, c='red',linestyle='dashed', linewidth=1);
plt.savefig('results/gwesplot_total_dataset.pdf')
