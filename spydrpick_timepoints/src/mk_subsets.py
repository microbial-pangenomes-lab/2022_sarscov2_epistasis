import random
import pickle
import sys
import numpy
import pandas as pd
#import matplotlib.dates as mdates

# Generate list of dates
count=0
list_dates = []

for i in range(4):
    year = str(2019 + i)
    print(year)
    if year=='2019':
        month = '12'
        year_month = year + '-' + month
        list_dates.append(year_month)
    elif year == '2022':
        for i in range(1,7):
            month = '0' + str(i)
            year_month= year + '-' + month
            list_dates.append(year_month)
    else:
        for i in range(1,13):
            if i < 10:
                month = '0' + str(i)
                year_month= year + '-' + month
            else:
                month = str(i)
                year_month = year + '-' + month
            list_dates.append(year_month)

# Dataframe for metadata (sort per date)
meta=pd.read_csv('../data/tree_subset.txt', sep=',',header=None)
meta.columns = ['name', 'x1', 'date']
#meta[['date']] += '-01'
#print(meta)
# Bool mask for the dates which have no all 3 spots filled yyyy/mm/dd
bool_mask = [len(x.split('-')) == 3 for x in meta['date']]
meta = meta[bool_mask]

# Convert date to datetime object to sort properly then reconvert it to string
meta['date'] = pd.to_datetime(meta['date'])
meta['date'] = meta['date'].dt.to_period('M')
meta=meta.sort_values(by='date')
meta['date'] = meta['date'].astype(str)
#print(meta)

# FIRST: take the first 4 months
meta2 = meta[(meta['date'] == '2019-12') | (meta['date'] == '2020-01') | (meta['date'] == '2020-02') | (meta['date'] == '2020-03') | (meta['date'] == '2020-04')]
choices = {}

#  -------- Or Load it directly from table ----------
# Extract 1000 random sequences from this first subset
#meta2=pd.read_csv('tmp.csv', header=None, sep='\t')
#meta2.columns = ['name', 'x1', 'date']
# create a random generator from numpy.random
meta2 = meta2.reset_index()
#print(meta2)
rng = numpy.random.default_rng()
selected_idx = rng.choice(meta2.index,size=1000,replace=False)
#print(selected_idx)
meta3 = meta2.iloc[selected_idx,:]
#print(meta3)
#print(set(list(meta3['date'])))
# Create a dict with key = 'name' and value = 'date'
sequences = dict(zip(meta3['name'],meta3['date']))
#print(len(sequences))


# Create all the other subsets ---> if the code works then start using only this loop without creating 1st subset and 2 subset (I will merge the first 5 months after)
for i in list_dates[5:]:
    #print(i)
    abc = meta[meta['date'] == i].reset_index()
    print(abc.info())
    rng = numpy.random.default_rng()
    selected_idx = rng.choice(abc.index,size=1000, replace=False)
    meta_final = abc.iloc[selected_idx]
    #print(meta_final)
    sequences.update(dict(zip(meta_final['name'],meta_final['date'])))
    print(set(list(meta_final['date'])))
#print(len(sequences))
pickle.dump(sequences,open('../data/choices_new.pkl','wb'))


#print(choices)
choices = pickle.load(open('../data/choices_new.pkl','rb'))

# Filter again sequences using the dict to select them and create a file for each inside subset folder
# Filter again sequences using the dict to select them and create a file for each inside subset folder
seqs = 0
s = []
sid = ''
kept = 0
# Loop on all data
for l in sys.stdin:
    if l.startswith('>') and len(s) != 0:
#        print(choices[l[1:].rstrip()])
        if sid[1:] in choices.keys():
            namefile= choices[sid[1:]] + '.fasta'
            print(namefile)
            with open(namefile,'a') as fasta:
                line1=''.join([sid,'\n'])
                line2 = ''.join(s)
                fasta.write(line1)
                fasta.write(line2)
                fasta.write('\n')
                kept += 1
            s = []
            sid = l.rstrip()
            seqs += 1
    elif l.startswith('>'):
        sid = l.rstrip()
        print(sid)
        continue
    elif sid[1:] in choices.keys() and len(s) == 0:
        s.append(l.rstrip())
        #print(s)
if sid[1:] in choices.keys():
    line1=''.join([sid,'\n'])
    line2 = ''.join(s)
    fasta.write(line1)
    fasta.write(line2)
    fasta.write('\n')
    kept += 1

