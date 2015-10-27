#!/usr/bin/env perl

# Find all contacts beween domains..

import sys, os, re, string
import argparse
from os.path import expanduser

home = expanduser("~")
sys.path.append(home + '/bioinfo-toolbox/parsing')
sys.path.append(home + '/git/bioinfo-toolbox/parsing')

import parse_contacts
import numpy as np
import matplotlib
matplotlib.use('Agg')


sep=5

contacts = parse_contacts.parse(open(c_filename, 'r'), sep)
contacts_np = parse_contacts.get_numpy_cmap(contacts)
contacts_np = contacts_np[start:end,start:end]

for i in range(len(contacts)):
    score = contacts[i][0]
    c_x = contacts[i][1] - 1
    c_y = contacts[i][2] - 1
        
    # only look at contacts within given range
    # default: take full sequence range into account
    if c_x < start or c_x >= end:
        continue
    if c_y < start or c_y >= end:
        continue
    if c_y-c_x < start or c_y >= end:
        continue
    
    if c_x < domain
        
        pos_diff = abs(c_x - c_y)
        too_close = pos_diff < 5
    
if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Plot protein residue contact maps.')
    p.add_argument('-t', '--threshold', default=-1, type=float)
    p.add_argument('--start', default=0, type=int)
    p.add_argument('--end', default=-1, type=int)
    p.add_argument('--sep', default=5, type=int)
    p.add_argument('--domain', default=-1, type=int)
