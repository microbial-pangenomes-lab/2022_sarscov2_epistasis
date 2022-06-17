#!/usr/bin/python3

import sys
import os
import gzip

# take file
def main():
    with open('out_spydrpick.txt', 'wt') as o:
        with gzip.open(sys.argv[1], 'rt') as f:
            for l in f.readlines():
                if not l.startswith('pos'):
                    A = l.rstrip().split('\t')[0]
                    B = l.rstrip().split('\t')[1]
                    C =l.rstrip().split('\t')[2]
                    D = l.rstrip().split('\t')[3]
                    absdiff = abs(int(A) - int(B))
                    line = ' '.join([A,B,str(absdiff),D,C]) + ' \n'
                    o.write(line)

if __name__ == '__main__':
    main()
