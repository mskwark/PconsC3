#!/usr/bin/env python

import sys, subprocess, os, multiprocessing
from distutils.spawn import find_executable

execdir = os.path.dirname(os.path.realpath(__file__))
julia = find_executable('julia')
if not julia:
    for d in ['/home/skwarkmj/sw/julia/bin', '.', '~/julia']:
        if os.path.exists(d + '/julia'):
            julia = d + '/julia'
if not julia:
    sys.stderr.write('Cannot find Julia!\n')
    sys.exit(0)


print "Using: ",julia

cpus = multiprocessing.cpu_count()

alignment = sys.argv[1]

if not os.path.exists(alignment):
    print "Error alignment do not exist"
    sys.exit(0)

stem = alignment[:alignment.rfind('.')]

if os.path.exists(stem + '.gdca'):
    print "Error output file do exist"
    sys.exit(0)

a = subprocess.check_output('{0:s}  {1:s}/rungDCA.jl {2:s} {3:s}.gdca'.format(julia, execdir, alignment, stem), shell=True)

f = open(stem + '.gneff', 'w')
f.write(a)
f.close()
