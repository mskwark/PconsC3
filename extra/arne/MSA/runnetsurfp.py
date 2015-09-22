#!/usr/bin/env python
from localconfig import *
from datetime import datetime
import sys, subprocess, os, math

def check_output(command):
    print ' '.join(command)
    return subprocess.Popen(command, stdout=subprocess.PIPE).communicate()[0]


seqfile=sys.argv[1]
sys.stderr.write(str(datetime.now()) + ': running NetSurfP\nThis may take quite a few minutes!\n')
t = check_output([netsurf, '-i', seqfile, '-a'])
f = open(seqfile + '.rsa', 'w')
f.write(t)
f.close()
