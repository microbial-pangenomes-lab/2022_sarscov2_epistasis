import sys


if __name__ == '__main__':
    keep = {l.rstrip().split(',')[0]
            for l in open(sys.argv[1])}

    sys.stderr.write(str(len(keep)))
    sys.stderr.write('\n')

    seqs = 0
    s = []
    sid = ''
    kept = 0
    for l in sys.stdin:
        if l.startswith('>') and len(s) != 0:
            if not seqs % 100:
                sys.stderr.write(f'Prefilter: {seqs} {kept}\n')
            if sid[1:] in keep:
                print(sid)
                print(''.join(s))
                kept += 1
            s = []
            sid = l.rstrip()
            seqs += 1
        elif l.startswith('>'):
            sid = l.rstrip()
            continue
        elif sid[1:] in keep or len(s) == 0:
            s.append(l.rstrip())
    if sid[1:] in keep:
        print(sid)
        print(''.join(s))
        kept += 1
    seqs += 1
    sys.stderr.write(f'Prefilter: {seqs} {kept}\n')
