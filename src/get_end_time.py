#!/usr/bin/python


import os
import sys
import argparse


def get_options():
    description = 'Get exact date for the time weight thing'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('yearmonth', nargs='?', default=None)

    return parser.parse_args()


if __name__ == '__main__':
    args = get_options()

    if args.yearmonth is None:
        print('')
        sys.exit(0)

    year, month = args.yearmonth.split('-')
    month = int(month)

    # 30 di conta novembre con april giugno e settembre
    if month in (11, 4, 6, 9):
        day = 30
    # di 28 ce n'e' uno
    elif month == 2:
        day = 28
    # tuti gli altri ne han 31
    else:
        day = 31
    print(f'{year} {month} {day}')
