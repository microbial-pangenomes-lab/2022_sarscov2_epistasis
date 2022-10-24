import os
import sys


if __name__ == '__main__':
    a = {}

    names = {l.rstrip()
             for l in open(sys.argv[1])}

    sys.stderr.write(str(len(names)))
    sys.stderr.write('\n')

    seqs = 0
    s = ''
    sid = ''

    kept = 0
    for l in sys.stdin:
        if l.startswith('>') and s != '':
            if not seqs % 100:
                sys.stderr.write(f'Keep and de-gap: {seqs} {kept}\n')
            if sid in names:
                s = s.replace('-', '')
                print(f'>{sid}\n{s}')
                kept += 1
            s = '' 
            sid = l.rstrip()[1:]
            seqs += 1
        elif l.startswith('>'):
            sid = l.rstrip()[1:]
            continue
        else:
            s += l.rstrip()
    if sid in names:
        s = s.replace('-', '')
        print(f'>{sid}\n{s}')
        kept += 1
    seqs += 1
    sys.stderr.write(f'Keep and de-gap: {seqs} {kept}\n')
