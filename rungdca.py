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


print julia

cpus = multiprocessing.cpu_count()

alignment = sys.argv[1]

if not os.path.exists(alignment):
    sys.stderr.write('Input alignment {:s} does not exist.\n'.format(alignment))
	sys.exit(0)

stem = alignment[:alignment.rfind('.')]
strength = sys.argv[2]

if os.path.exists(stem + '.gdca'):
    sys.stderr.write('Output file {:s}.gdca exists.\n'.format(stem))
    sys.exit(0)

a = subprocess.check_output('{0:s} -p {1:d} {2:s}/rungDCA.jl {3:s} {4:s}..gdca'.format(julia, cpus, execdir, alignment, stem), shell=True)

f = open(stem + '.gneff', 'w')
f.write(a)
f.close()
