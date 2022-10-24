#!/usr/bin/python3

from treeswift import read_tree_newick
import sys

# The command should be launched as:
# python3 --tree $treefile --basename $path_to_sample

#tree = read_tree_newick('~/spydrpick/2022_sarscov2_epistasis/GISAID-hCoV-19-phylogeny-2022-07-05/global.tree')
tree =read_tree_newick(sys.argv[2])
#print(tree)
# Extract names of each subset
names = set()
path_to_filtered_names = '/home/gabriel/spydrpick/2022_sarscov2_epistasis/subset_workflow/data/subsets/new_subset_total_out/' + sys.argv[4] + '/filtered_dca.fasta.names'
print('Where to find filtered files:', path_to_filtered_names)

with open(path_to_filtered_names, 'rt') as file:
    for n in file:
        names.add(n.rstrip())

# Prune tree with the selected names
tree2 = tree.extract_tree_with(names)

'''
# write tree to file
# into new folder
path_to_new_folder = '/home/gabriel/spydrpick/2022_sarscov2_epistasis/subset_workflow/data/subsets/new_subset_total_out/new_out_weights'
name_new_file = sys.argv[4] + '.tree'
path_to_new_file= path_to_new_folder + '/' + sys.argv[4] + '/' + name_new_file
print('Where to write new file tree: ', path_to_new_file)
tree2.write_tree_newick(path_to_new_file)
'''
