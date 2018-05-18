#!/usr/bin/env python

import sys, subprocess, os, math
import string as s
from localconfig import *
from datetime import datetime
import multiprocessing

def check_output(command):
    print ' '.join(command)
    return subprocess.Popen(command, stdout=subprocess.PIPE).communicate()[0]


sys.stderr.write('\nTesting dependencies...\n')

### Check jackhmmers ###
### Check Jackhmmer ###
try:
    f = open(os.devnull, "w")
    x  = subprocess.call([jackhmmer, '-h'], stdout=f, stderr=f)
    f.close()
except Exception as e:
    sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
    sys.stderr.write(jackhmmer + ' -h \n\n')
    sys.stderr.write('Chosen jackhmmer binary does not seem to work!\n')
    sys.exit(1)


if '-c' in sys.argv:
    idx = sys.argv.index('-c')
    try:
        n_cores = int(sys.argv[idx+1])
    except:
        print 'Number of cores -c must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_cores)
        sys.exit(1)
    del sys.argv[idx]
    del sys.argv[idx]
else:
    n_cores = 1

if '-e' in sys.argv:
    idx = sys.argv.index('-e')
    evalue = sys.argv[idx+1]
    del sys.argv[idx]
    del sys.argv[idx]
else:
    evalue = 0.001

if '-db' in sys.argv:
    idx = sys.argv.index('-db')
    database = sys.argv[idx+1]
    del sys.argv[idx]
    del sys.argv[idx]
else:
    database = uniref

if '-name' in sys.argv:
    idx = sys.argv.index('-name')
    name = sys.argv[idx+1]
    del sys.argv[idx]
    del sys.argv[idx]
else:
    name = 'JH'+str(evalue)


if len(sys.argv) != 2:
    sys.stderr.write('Incorrect number of command line arguments.\n')
    sys.stderr.write('Usage: ' + sys.argv[0] + ' [-c <n_cores>] [-e E-value-cutoff]  [-db database]  <sequence file>\n\n')
    sys.exit(0)

seqfile = os.path.abspath(sys.argv[1])

if not os.path.exists(database):
    sys.stderr.write('\n' + database + 'does not exist\n')
    sys.exit(1)
if not os.path.exists(seqfile):
    sys.stderr.write('\n' + seqfile + 'does not exist\n')
    sys.exit(0)


f = open(seqfile).read()

if os.path.exists(seqfile + '.fasta'):
    subprocess.call(['mv', seqfile + '.fasta', seqfile +'.bak'])

f2 = open(seqfile +'.fasta', 'w')
if f[0] != '>':
    f2.write('>target\n' + f +'\n')
else:
    x = f.split('\n')
    if len(x[0]) > 6:
        target = x[0][1:5] + x[0][6]
    f2.write('>target\n' + "".join(x[1:]) + '\n')
f2.close()

sys.stderr.write(str(datetime.now()) + ' ' + name + ': generating Jackhmmer alignment\nThis may take quite a few minutes!\n ')
t = check_output([jackhmmer, '--cpu', str(n_cores), '-N', '5', '--incE', str(evalue), '-E', str(evalue), '-A', seqfile +'.' + name + '.sto', seqfile + '.fasta', database])
check_output([reformat, 'sto', 'a3m', seqfile + '.' + name + '.sto', seqfile + '.' + name + '.a3m'])
#check_output(['rm', seqfile + '.' + name + '.ali'])

