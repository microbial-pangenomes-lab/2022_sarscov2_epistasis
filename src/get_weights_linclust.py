import sys
import numpy as np
import os
import argparse
__version__ = '0.0.1'

def read_samples(samples_filename):
    sys.stderr.write(f'Load: Samples\n')
    sample_names = np.array([x.rstrip()
                             for x in open(samples_filename)])

    return sample_names

def get_weights_clusters(clusters_file, sample_names):
    sys.stderr.write('Load: Clusters\n')
    
    # from cluster file, create a dictionary made by:
    # Cluster name as key
    # Set of all sequences into each cluster as value
    d = {}
    for l in open(clusters_file):
        cname, name = l.rstrip().split()
        d[cname] = d.get(cname, set())
        d[cname].add(name)
    
    # iterate through previous dictionary and create a new one made by:
    # Each sequence name is a key
    # Relative cluster is the value
    r = {}
    for c, names in d.items():
        for name in names:
            r[name] = c

    # Create an array in which zeros as long as sample names table
    # Fill it with weight which is calculated as 1 / $(lenght of each value)
    weights = np.zeros(len(sample_names))
    for i, name in enumerate(sample_names):
        weights[i] = 1 / len(d[r[name]])

    # normalise (not really necessary)
    #weights = weights / np.sum(weights)

    return weights


def get_options():
    import argparse

    description = 'Extract weights from the output of mmseqs linclust'
    parser = argparse.ArgumentParser(description=description)

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "-s",
        "--sample",
        dest="sample_file",
        required=True,
        help="Samples in the alignment")
    io_opts.add_argument(
        "--clusters",
        dest="clusters_file",
        default=None,
        help="Cluster file from linclust")

    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()

    return (args)


def main():
    args = get_options()

    sample_names = read_samples(args.sample_file)

    # get weights
    weights = get_weights_clusters(args.clusters_file, sample_names)
    
    for w in weights:
        print(w)

    return


if __name__ == '__main__':
    main()
