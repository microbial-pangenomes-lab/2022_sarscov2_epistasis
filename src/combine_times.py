#!/usr/bin/python


import os
import sys
import argparse


def get_options():
    description = 'Combine time subsets'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('input',
                        help='Input directory (contains fasta files)')
    parser.add_argument('output',
                        help='Output directory (must be present)')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_options()

    # sort input files
    inputs = sorted([x for x in os.listdir(args.input)
                     if x.endswith('.fasta')],
                    key=lambda x: (int(x.split('-')[0]),
                                   int(x.split('-')[1].split('.')[0])))
    for i in range(len(inputs)):
        fname = inputs[i]
        sys.stderr.write(f'Writing {fname}\n')
        f = open(f'{args.output}/{fname}', 'w')
        for j in range(i+1):
            for l in open(f'{args.input}/{inputs[j]}'):
                f.write(l)
                if not l.endswith('\n'):
                    f.write('\n')
        f.close()
