#!/usr/bin/env python

import sys, os, re

infilef = sys.argv[1]
infile = open(infilef)

cutoff=0.75
counter = 0
for l in infile:
    if '>' in l and not counter == 0:
        m=re.sub("-",'',upperseq)
        lenali=len(m)
        fraction=float(lenali)/float(length)
        if fraction>cutoff:
            sys.stdout.write('\n>sequence{0:07d}/1-100\n'.format(counter))
            sys.stdout.write(upperseq)
            counter += 1
        upperseq=''
    elif '>' not in l and counter >1:
        l = l.strip()
        upperseq = ''.join([c for c in l if not c.islower()])
        upperseq = upperseq.replace('X', '-')
    elif '>' in l and counter == 0:
        sys.stdout.write('>target/1-100\n')
        counter += 1
    elif '>' not in l and counter == 1:
        l = l.strip()
        upperseq = ''.join([c for c in l if not c.islower()])
        upperseq = upperseq.replace('X', '-')
        length=len(upperseq)

sys.stdout.write('\n')
