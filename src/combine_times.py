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
        cmd = 'cat'
        for j in range(i+1):
            cmd += ' '
            cmd += f'{args.input}/{inputs[j]}'
        cmd += f' > {args.output}/{fname}'
        os.system(cmd)
