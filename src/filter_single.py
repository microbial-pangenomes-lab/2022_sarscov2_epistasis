import os
import sys
import hashlib
import binascii
import numpy as np
from numba import jit


@jit(nopython=True)
def get_reduced_sequence(s, keep):
    a = np.frombuffer(s, dtype=np.uint8).copy()
    a[a == 65] = 1
    a[a == 71] = 2
    a[a == 67] = 3
    a[a == 84] = 4
    a[a > 4] = 0
    return a[keep]


if __name__ == '__main__':
    a = {}

    keep = {int(x) - 1
            for l in open(sys.argv[1])
            for x in l.rstrip().split(',')}
    nkeep = np.array(sorted(keep))

    sys.stderr.write(str(len(keep)))
    sys.stderr.write('\n')

    seqs = 0
    s = ''
    sid = ''

    unique = set()
    for l in sys.stdin:
        if l.startswith('>') and s != '':
            if not seqs % 100:
                sys.stderr.write(f'Filter and transform: {seqs} {len(unique)}\n')
            b = get_reduced_sequence(s.upper().encode(), nkeep)
            b = b.reshape((1, b.shape[0]))
            h = binascii.b2a_base64(hashlib.md5(b).digest())
            if h not in unique:
                a[len(unique)] = b
                unique.add(h)
                print(sid)
            s = ''
            sid = l.rstrip()[1:]
            seqs += 1
        elif l.startswith('>'):
            sid = l.rstrip()[1:]
            continue
        else:
            s += l.rstrip()
    b = get_reduced_sequence(s.upper().encode(), nkeep)
    b = b.reshape((1, b.shape[0]))
    h = binascii.b2a_base64(hashlib.md5(b).digest())
    if h not in unique:
        a[len(unique)] = b
        unique.add(h)
        print(sid)
    seqs += 1
    sys.stderr.write(f'Filter and transform (finished): {seqs} {len(unique)}\n')
    a = np.concatenate(tuple(a.values()), axis=0)
    np.savez_compressed(sys.argv[2], a)
