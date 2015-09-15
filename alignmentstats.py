#!/usr/bin/env python

import sys, subprocess, os, multiprocessing
from distutils.spawn import find_executable

execdir = os.path.dirname(os.path.realpath(__file__))
cdhit = find_executable('cd-hit')

if not cdhit:
    cdhit = find_executable('cdhit')
if not cdhit:
    for d in ['/home/skwarkmj/sw/cd-hit', '.']:
        if os.path.exists(d + '/cd-hit'):
            cdhit = d + '/cd-hit'

if not cdhit:
    sys.stderr.write('Cannot find CD-HIT!\n')
    sys.exit(0)

def cluster(a, threshold):
    if threshold < 0.7:
        b = subprocess.check_output(cdhit + ' -n 3 -i ' + a + ' -c {:6.4f}'.format(threshold) + ' -o ' + a + '.cd{:d}'.format(int(threshold*100)), shell=True)
        subprocess.call('rm ' + a + '.cd{:d}*'.format(int(threshold*100)), shell=True)
    else:
        b = subprocess.check_output(cdhit +' -i ' + a + ' -c {:6.4f}'.format(threshold) + ' -o ' + a + '.cd{:d}'.format(int(threshold*100)), shell=True)
        subprocess.call('rm ' + a + '.cd{:d}*'.format(int(threshold*100)), shell=True)
    return int(b[b.rfind('finished'):].split()[1])

for f in sys.argv[1:]:
    stem = f[:f.rfind('.')]
    if os.path.exists(stem + '.stats'):
        a = subprocess.check_output('wc ' + stem + '.stats', shell=True).split()[0]
        if int(a[0]) != 6:
            pass
        else:
            continue
    sys.stderr.write('Reading in the alignment and analysing\n')
    a2 = open(f).readlines()
    a = set()
    for l in a2:
        if l.find('>') > -1:
            continue
        else:
            if len(l.strip()) < 2:
                continue
            a.add(l.strip())
    a = list(a)
    a1 = ''.join(a)
    g = open(stem + '.stats', 'w')
    g.write('Length: {:d}\n'.format(len(a[0])))
    g.write('Count: {:d}\n'.format(len(a)))
    g.write('Fraction of gaps: {:6.4f}\n'.format(a1.count('-')/(float(len(a) * len(a[0])))))
    sys.stderr.write('Clustering to 90% identity\n')
    g.write('Count at 90% identity: {:d}\n'.format(cluster(f, 0.9)))
    sys.stderr.write('Clustering to 70% identity\n')
    g.write('Count at 70% identity: {:d}\n'.format(cluster(f, 0.7)))
    sys.stderr.write('Clustering to 50% identity\n')
    g.write('Count at 50% identity: {:d}\n'.format(cluster(f, 0.5))) 
