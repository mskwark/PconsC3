#!/usr/bin/env python

# Trim multiple sequence alignment to first sequecne in file

import sys, os

infilef = sys.argv[1]

if len(sys.argv) > 2:
        lookfor = sys.argv[2]
else:
	lookfor = 'target'

onetonumber = 'ARNDCEQGHILKMFPSTWYV-'

infile = open(infilef).read().split('\n')

sequences = {}
key = 'n/a'
#trimto = 'target'
trimto = 'FirstLine'
#Using  firs sequence

for l in infile:
    if l.find('>') == 0:
        key = l[1:]
        continue
    if trimto == "FirstLine":
        trimto = key
    if key in sequences:
        sequences[key] = sequences[key] + l
    else:
        sequences[key] = l


keep = []

for i in range(len(sequences[trimto])):
	if sequences[trimto][i] == '-':
		continue
	else:
		keep.append(i)

seq = ''
print ">"+trimto
for i in range(len(keep)):
    seq = seq + sequences[trimto][keep[i]]
print seq


sequence = ""
key = 'n/a'
for key in sequences:
    if key == trimto:
        continue
    print ">"+key
    seq = ""
    sequence = ""
    for i in range(len(keep)):
        seq = seq + sequences[key][keep[i]]
    print seq
    sequence = sequence + l.replace('X', '-')


