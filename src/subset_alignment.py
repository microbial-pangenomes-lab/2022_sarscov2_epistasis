#!/usr/bin/python


import sys
import argparse


def get_options():
    description = 'Subset and clean an alignment to a specified region'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--start', type=int, default=21562,
                        help='First base to consider '
                             '(0-based index, default %(default)d)')
    parser.add_argument('--end', type=int, default=25384,
                        help='Last base to consider '
                             '(0-based index, default %(default)d)')

    return parser.parse_args()


if __name__ == '__main__':
    args = get_options()

    seq = ''
    allowed  = '-AGCTagct'
    for l in sys.stdin:
        if l.startswith('>'):
            if seq != '':
                seq = ''.join([x.upper()
                               if x in allowed else '-'
                               for x in seq])
                seq = seq[args.start: args.end]
                print(f'>{sid}\n{seq}')
                seq = ''
            sid = l[1:].rstrip().split()[0]
        else:
            seq += l.strip()
