#!/usr/bin/env python

import os, sys, subprocess

def cluster(a, threshold):
	if threshold < 0.7:
		b = subprocess.check_output('/triton/ics/work/mjs/cd-hit-v4.5.4-2011-03-07/cd-hit -n 3 -i ' + a + ' -c {:6.4f}'.format(threshold) + ' -o ' + a + '.cd{:d}'.format(int(threshold*100)), shell=True)
		subprocess.call('rm ' + a + '.cd{:d}'.format(int(threshold*100)), shell=True)
	else:
		b = subprocess.check_output('/triton/ics/work/mjs/cd-hit-v4.5.4-2011-03-07/cd-hit -i ' + a + ' -c {:6.4f}'.format(threshold) + ' -o ' + a + '.cd{:d}'.format(int(threshold*100)), shell=True)
		subprocess.call('rm ' + a + '.cd{:d}'.format(int(threshold*100)), shell=True)
	return int(b[b.rfind('finished'):].split()[1])
	
	result = set()
	for s in a:
		add = True
		length  = len(s)
		for r in result:
			c = 0
			if len(r) != len(s):
				continue
			limit = threshold*min(len(s.replace('-', '')), len(r.replace('-', '')))
			for i in range(length):
				if s[i] == r[i]:
					c+=1
				if c >= limit:
					add = False
					break
			if not add:
				break
		if add:
			result.add(s)
	return result


for f in sys.argv[1:]:
	if f.find('.trimmed') != len(f) - len('.trimmed'):
		continue
	if os.path.exists(sys.argv[1] + '/' + f[:f.rfind('.trimmed')] + '.stats'):
		a = subprocess.check_output('wc ' + sys.argv[1] + '/' + f[:f.rfind('.trimmed')] + '.stats', shell=True).split()[0]
		if int(a[0]) != 6:
			pass
		else:
			continue
	a2 = open(f).readlines()
	a1 = ''
	a = []
	for l in a2:
		if l.find('>') > -1:
			continue
		else:
			if len(l.strip()) < 2:
				continue
			a1 += l
			a.append(l.strip())
	g = open(f[:f.rfind('.trimmed')] + '.stats', 'w')
	g.write('Length: {:d}\n'.format(len(a[0])))
	g.write('Count: {:d}\n'.format(len(a)))
	g.write('Fraction of gaps: {:6.4f}\n'.format(a1.count('-')/(float(len(a) * len(a[0])))))
	g.write('Count at 90% identity: {:d}\n'.format(cluster(f[:-len('.trimmed')] + '.trimmed', 0.9)))
	g.write('Count at 70% identity: {:d}\n'.format(cluster(f[:-len('.trimmed')] + '.trimmed', 0.7)))
	g.write('Count at 50% identity: {:d}\n'.format(cluster(f[:-len('.trimmed')] + '.trimmed', 0.5)))
	
