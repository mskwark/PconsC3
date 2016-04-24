#!/usr/bin/env python
from localconfig import *
from datetime import datetime
import sys, subprocess, os, math
import tempfile

def check_output(command):
    print ' '.join(command)
    return subprocess.Popen(command, stdout=subprocess.PIPE).communicate()[0]


seqfile=sys.argv[1]
sys.stderr.write(str(datetime.now()) + ': running NetSurfP\nThis may take quite a few minutes!\n')
tmp=tempfile.NamedTemporaryFile(mode='w',suffix='.fa',delete=False)
tmpfile=tmp.name
counter =0
with open(seqfile,'r') as seq:
    for line in seq:
        if '>' in line:
            tmp.write('>target\n')
        elif '>' not in line:
            tmp.write(line)
tmp.close()
t = check_output([netsurf, '-i', tmpfile,'-a','-v','-k'])
f = open(seqfile + '.rsa', 'w')
f.write(t)
f.close()
#os.remove(tmpfile)


