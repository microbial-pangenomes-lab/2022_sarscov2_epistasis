#!/usr/bin/python3


import argparse
from treeswift import read_tree_newick
import sys


def get_options():
    description = 'Prune a tree'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--tree',
                        required=True,
                        help='Phylogenetic tree')
    parser.add_argument('--name',
                        required=True,
                        help='Sample names to keep')
    parser.add_argument('--out',
                        required=True,
                        help='Output file for the resulting tree')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_options()

    print(f'Reading {args.tree}')
    tree = read_tree_newick(args.tree)
    # Extract names of each subset
    names = set()
    print('Where to find filtered files:', args.name)

    with open(args.name, 'rt') as file:
        for n in file:
            names.add(n.rstrip())

    # Prune tree with the selected names
    print('Pruning tree')
    tree2 = tree.extract_tree_with(names)

    print('Where to write new file tree:', args.out)
    tree2.write_tree_newick(args.out)
