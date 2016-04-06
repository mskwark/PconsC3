#!/usr/bin/env python
import os
import sys
import multiprocessing
import subprocess
from os.path import expanduser
home = expanduser("~")


if __name__ == '__main__':
    print 'Please do not run me! Use run_pconsc.py'
    print '\n\tYours sincerely,\n\n\t', sys.argv[0]
    sys.exit(0)

# Directory where PconsC in the distributable package is located
#root = os.path.dirname(os.path.abspath(sys.argv[0])) + '/'
#root = home +'/git/PconsC3//'
#print root

# Directory to PconsC3 (i.e. this one)
PconsC3 = home +'/git/PconsC3//'

########################################################
### Please adjust the following paths to your system ###
########################################################

if os.path.exists("/proj/bioinfo/"):
    print " We are at triolith"
    root = "/proj/bioinfo/software/PconsC2-extra/hhsuite-2.0.16/"
### Jackhmmer executable ###
#jackhmmer = root + 'dependencies/hmmer-3.0/src/jackhmmer'
    jackhmmer = '/proj/bioinfo/software/PconsC2-extra/hmmer-3.1b1-linux-intel-x86_64/binaries/jackhmmer'
    uniref = '/proj/bioinfo/data/uniref90.fasta'

### HHblits executable ###
#hhblits = root + 'dependencies/hhsuite-2.0.16/bin/hhblits'
    hhblits = '/software/apps/hhsuite/2.0.16/gcc01/hhsuite-2.0.16/bin/hhblits'
    hhdatabase = '/proj/bioinfo/data/hhsuite/hhsuite_dbs/uniprot20_2013_03/uniprot20_2013_03'

### PSICOV executable ###
#psicov = root + 'dependencies/psicov-1.11/psicov'
#psicov = '/usr/local/bin/psicov'
    psicov= '/proj/bioinfo/software/PconsC2-extra/psicov2/psicov2'

### NetSurfP executable ###
#netsurf = root + 'dependencies/netsurfp-1.0/netsurfp'
    netsurf = '/proj/bioinfo/software/PconsC2-extra/netsurfp-1.0/netsurfp'

### PSIPRED executable ###
#psipred = root + 'dependencies/psipred/runpsipred'
#psipred = '/scratch/arne/PconsC2-extra/psipred/runpsipred'
    psipred = '/proj/bioinfo/software/PconsC2-extra/psipred/bin/psipred'

### MATLAB executable ###
# Please set this variable to None if you don't have access to matlab. 
# PconsFold will then try to use the compiled version. 
#matlab = '/sw/apps/matlab/x86_64/8.1/bin/matlab'
    matlab = None
#matlab= '/pdc/vol/matlab/r2012a/bin/matlab'

### Path to MATLAB compiler ###
# Only needed if matlab is not available.
    matlabdir = '' 

# Directory to PconsC3 scripts (i.e. this one)
    PconsC3 = '/proj/bioinfo/software/PconsC3/'




elif os.path.exists("/pfs/nobackup/home/"):
    print "We are at HPC2N"
    root = "/pfs/nobackup/home/a/arnee/Software/PconsC2-extra/hhsuite-2.0.16-linux-x86_64/"
### Jackhmmer executable ###
#jackhmmer = root + 'dependencies/hmmer-3.0/src/jackhmmer'
    jackhmmer = '/home/a/arnee/bin/jackhmmer'
    uniref = '/pfs/nobackup/home/a/arnee/data/uniref90.fasta'

### HHblits executable ###
#hhblits = root + 'dependencies/hhsuite-2.0.16/bin/hhblits'
    hhblits = '/pfs/nobackup/home/a/arnee/Software/PconsC2-extra/hhsuite-2.0.16-linux-x86_64/bin/hhblits'
    hhdatabase = '/pfs/nobackup/home/a/arnee/data/hhsuite_dbs/uniprot20_2013_03/uniprot20_2013_03'

### PSICOV executable ###
#psicov = root + 'dependencies/psicov-1.11/psicov'
#psicov = '/usr/local/bin/psicov'
    psicov= '/pfs/nobackup/home/a/arnee/Software/PconsC2-extra/psicov2/psicov2'

### NetSurfP executable ###
#netsurf = root + 'dependencies/netsurfp-1.0/netsurfp'
    netsurf = '/pfs/nobackup/home/a/arnee/Software/PconsC2-extra/netsurfp-1.0/netsurfp'

### PSIPRED executable ###
#psipred = root + 'dependencies/psipred/runpsipred'
#psipred = '/scratch/arne/PconsC2-extra/psipred/runpsipred'
    psipred = '/pfs/nobackup/home/a/arnee/Software/PconsC2-extra/psipred/bin/psipred'

### MATLAB executable ###
# Please set this variable to None if you don't have access to matlab. 
# PconsFold will then try to use the compiled version. 
#matlab = '/sw/apps/matlab/x86_64/8.1/bin/matlab'
    matlab = None
#matlab= '/pdc/vol/matlab/r2012a/bin/matlab'

### Path to MATLAB compiler ###
# Only needed if matlab is not available.
    matlabdir = '' 

# Directory to PconsC3 scripts (i.e. this one)
    PconsC3 = '/home/a/arnee/git/PconsC3/'

elif os.path.exists("/scratch/arne/"):
    print "We are at local machine"
    root = "/scratch/arne/PconsC2-extra/hhsuite-2.0.16-linux-x86_64/"
### Jackhmmer executable ###
#jackhmmer = root + 'dependencies/hmmer-3.0/src/jackhmmer'
    jackhmmer = '/usr/local/bin/jackhmmer'
    uniref = '/scratch/data/uniref90.fasta'

### HHblits executable ###
#hhblits = root + 'dependencies/hhsuite-2.0.16/bin/hhblits'
    hhblits = '/scratch/arne/PconsC2-extra/hhsuite-2.0.16-linux-x86_64/bin/hhblits'
    hhdatabase = '/scratch/data/hhsuite/hhsuite_dbs/uniprot20_2013_03/uniprot20_2013_03'

### PSICOV executable ###
#psicov = root + 'dependencies/psicov-1.11/psicov'
    psicov = '/usr/local/bin/psicov2'
#    psicov= '/pfs/nobackup/home/a/arnee/Software/PconsC2-extra/psicov2/psicov2'

### NetSurfP executable ###
#netsurf = root + 'dependencies/netsurfp-1.0/netsurfp'
#    netsurf = '/pfs/nobackup/home/a/arnee/Software/PconsC2-extra/netsurfp-1.0/netsurfp'
    netsurf = '/usr/local/bin/netsurfp'
### PSIPRED executable ###
#psipred = root + 'dependencies/psipred/runpsipred'
#psipred = '/scratch/arne/PconsC2-extra/psipred/runpsipred'
    psipred = '/usr/local/bin/psipred/bin/psipred'

### MATLAB executable ###
# Please set this variable to None if you don't have access to matlab. 
# PconsFold will then try to use the compiled version. 
#matlab = '/sw/apps/matlab/x86_64/8.1/bin/matlab'
#    matlab = None
    matlab= '/pdc/vol/matlab/r2012a/bin/matlab'

### Path to MATLAB compiler ###
# Only needed if matlab is not available.
    matlabdir = '' 

# Directory to PconsC3 scripts (i.e. this one)
    PconsC3 = '/home/arnee/git/PconsC3/'

########################################################
###  Please do not change anything below this line   ###
########################################################

# Paths to included scripts
trim2jones = root + 'scripts/a3mToJones.py'
trim2trimmed = root + 'scripts/a3mToTrimmed.py'

# Reformat script scavenged from HHsuite. Please cite the HHblits paper!
reformat = root + 'scripts/reformat.pl'
#reformat = '/pfs/nobackup/home/a/arnee/Software/PconsC2-extra/hhsuite-2.0.16-linux-x86_64/scripts/reformat.pl'

# Maximum amount of cores to use per default
n_cores = multiprocessing.cpu_count()

# Enable work-around for PSICOV not handling low complexity alignments
psicovfail = True

# Adjust plmdca path to either standalone or compiled, 
# depending on presence of matlab.
if matlab:
    plmdca = None # matlab licence present: do not use compiled version
#    plmdcapath = root + 'dependencies/plmDCA_symmetric_v2'
    plmdcapath = '/scratch/arne/PconsC2-extra/plmDCA_asymmetric_v2'
else:
    plmdca = os.path.join(PconsC3, 'dependencies/plmdca_standalone/2012/build01/bin/plmdca')
    plmdcapath = None

