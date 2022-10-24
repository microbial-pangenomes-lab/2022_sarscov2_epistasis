import sys
sys.setrecursionlimit(99999)
import numpy as np
import os
import argparse
__version__ = '0.0.1'
import dendropy as dpy

def read_samples(samples_filename):
    sys.stderr.write(f'Load: Samples\n')
    sample_names = np.array([x.rstrip()
                             for x in open(samples_filename)])

    return sample_names

def get_weights_phylogeny(tree_file, sample_names):
    sys.stderr.write('Load: Tree\n')
    # read in tree
    tree = dpy.Tree.get(path=tree_file,
                        schema="newick",
                        preserve_underscores=True)

    # check sample names match up
    #tip_labels = [t.label for t in tree.taxon_namespace]

    #if len(set(tip_labels).intersection(set(sample_names))) != len(sample_names):
    #    raise ValueError(
    #        "Sample name mismatch between tree file and presence/absence matrix"
    #    )

    sys.stderr.write('Load: Tree distances\n')
    # iterate top down to weight by edge length.
    for node in tree.preorder_node_iter():
        if node.parent_node is not None:
            node.total_weight = node.parent_node.total_weight + float(
                node.edge.length) / len(node.leaf_nodes())
        else:
            node.total_weight = 0

    # iterate through leaves getting weights
    weight_dict = {}
    for node in tree.leaf_node_iter():
        weight_dict[node.taxon.label] = node.total_weight

    weights = np.zeros(len(sample_names))
    for i, name in enumerate(sample_names):
        weights[i] = weight_dict[name]

    # normalise (not really necessary)
#    weights = weights / np.sum(weights)

    return weights


def get_options():
    import argparse

    description = 'Extract weights from a phylogenetic tree'
    parser = argparse.ArgumentParser(description=description)

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "-s",
        "--sample",
        dest="sample_file",
        required=True,
        help="Samples in the alignment")
    io_opts.add_argument(
        "--tree",
        dest="tree_file",
        default=None,
        help=("phylogeny in newick format for weigting samples to" +
              " control for population structure"))

    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()

    return (args)


def main():
    args = get_options()

    sample_names = read_samples(args.sample_file)

    # get weights
    weights = get_weights_phylogeny(args.tree_file, sample_names)
    
    for w in weights:
        print(w)

    return


if __name__ == '__main__':
    main()
