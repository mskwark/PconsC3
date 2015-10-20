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


cpus = multiprocessing.cpu_count()

cpus = 4
print 'Using {:d} CPUs'.format(cpus)

alignment = sys.argv[1]
strength = "0.02"

if len(sys.argv) > 2:
    strength = sys.argv[2]

if not os.path.exists(alignment):
	sys.exit(0)

stem = alignment[:alignment.rfind('.')]

if os.path.exists(stem + '.{:s}.plm20'.format(strength)):
    sys.exit(0)

a = subprocess.check_output('{0:s} -p {1:d} {2:s}/runplm.jl {3:s} {4:s} {5:s}.{4:s}.plm20'.format(julia, cpus, execdir, alignment, strength, stem), shell=True)
