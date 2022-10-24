#!/usr/bin/python3

import pandas as pd
from random import randint
import numpy as np

# Load the metadata csv into a pandas df
df = pd.read_csv('GISAID-hCoV-19-phylogeny-2022-07-05/metadata.csv', sep=',')
df['collection_date'] = pd.to_datetime(df['collection_date'])
#print(df)
# sort
df.sort_values(by='collection_date', inplace=True)

# Select first interval
first = (df['collection_date'] > '2019-12-01') & (df['collection_date'] <= '2020-04-30')
subdf1 = df.loc[first].reset_index()
#print(subdf1)
nrows=len(subdf1)
print(nrows)

# Obtain 1000 random indexes, without duplicates
def draw_random_indexes(num, maxcol):
    indexes_set = set([randint(0, maxcol) for _ in range(0,num)])
    # if duplicates are present, set will be shorter than the number of random numbers wanted: update the set object
    # till the length of set will be the same of the list of random numbers wanted.
    while len(indexes_set) != num:
        diff = num - len(indexes_set)
        indexes_set.update([randint(0,maxcol) for _ in range(0,diff)])
        #print(diff)
    return list(indexes_set)


# First subset
first = (df['collection_date'] > '2019-12-01') & (df['collection_date'] <= '2020-04-30')
subdf1 = df.loc[first].reset_index()
nrows = len(subdf1)
# Select 1000 random sequences
list_indexes = draw_random_indexes(1000, nrows)
subdf1_1000 = subdf1.iloc[list_indexes]

# save to a file
subdf1_1000[['accession_id', 'virus_name','collection_date']].to_csv('subset_workflow/subsets/subset1.csv',
                                                                     sep=',',
                                                                    index=False)

# Work with dates as indexes
df_date = df.set_index('collection_date')
df_date

# Create
#  Now create all the other subsets
total = subdf1_1000[['accession_id', 'virus_name','collection_date']]
# Make subsets
for y in [2020, 2021, 2022]:
    for month in range(1,13):
        if y == 2022 and month > 6:
            break
        print(month)
        # exclude months till april
        if y == 2020 and month <= 4:
            pass
        else:
            indexes = np.where((df_date.index.year == y) & (df_date.index.month == month))
            indexed_df = df.loc[indexes].reset_index()
            # random 1000 sequences
            list_indexes = draw_random_indexes(1000, len(indexed_df)-1)
            subsetted_indexed_df = indexed_df.loc[list_indexes]
            # concatenate the dfs
            total = pd.concat([total, subsetted_indexed_df], ignore_index=True)
            name_file = 'subset_workflow/subsets/subset_' + str(y) + '-' + str(month) + '.csv'
            total[['accession_id', 'virus_name','collection_date']].to_csv(name_file, sep=',', index=False)
