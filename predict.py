#!/usr/bin/env python2.7

import joblib
import sys
import os, random, time
import numpy as np

forestlocation = '/dev/shm/'

# maximum time per layer
maxtime = pow(10,6)

# fraction of trees to use (prediction time scales linearly with the number of trees, 
# while expected precision is roughly the same for values > 0.3

treefraction = 1.0

if not os.path.exists(forestlocation + '/tlayer0'):
    forestlocation =  os.path.dirname(os.path.realpath(__file__))

for i in range(5):
    abort = False
    if not os.path.exists(forestlocation + '/tlayer{:d}/tree.list'.format(i)):
        sys.stderr.write('Forest data for layer {:d} is missing.\n'.format(i))
        abort = True
    if abort:
        sys.exit(0)

firststart = time.time()
if len(sys.argv) != 9:
    print 'Usage: ' + sys.argv[0] + ' <files>'
    print 'That is:'
    print ' GaussDCA prediction'
    print ' plmDCA.jl predictions'
    print ' ML contact predictions'
    print ' NetSurf RSA'
    print ' SS file'
    print ' Alignment stats'
    print ' Alignment'
    print ' Outfile'
    sys.exit(1)

def parsePSIPRED(f):
    SSdict = {}
    try:
        x = open(f).read().split('\n')
    except:
        return SSdict
    for l in x:
        y = l.split()
        if len(y) != 6:
            continue
        i = int(y[0])
        SSdict[i] = [float(y[3]), float(y[4]), float(y[5])]
    return SSdict
        
def parseNetSurfP(f):
    netSurfdict = {}
    for l in open(f).readlines():
        al = []
        x = l.split()
        if l.find('#') == 0:
            continue
        if l[0] not in ['B', 'E']:
            y = ['E']
            y.extend(x)
            x = y
        for y in [4,6,7, 8, 9]:
            al.append(float(x[y]) )
        netSurfdict[ int(x[3] )] = al
    return netSurfdict

def parsePSSM(alignment):
    pssm = {}
    one2number = 'ARNDCEQGHILKMFPSTWYV-'
    bi = [ 0.0825, 0.0553, 0.0406, 0.0545, 0.0137, 0.0393, 0.0675, 0.0707, 0.0227, 0.0595, 0.0966, 0.0584, 0.0242, 0.0386, 0.0470, 0.0657, 0.0534, 0.0108, 0.0292, 0.0687 ]
    b = {}
    for i in one2number[:-1]:
        b[i] = bi[one2number.find(i)] 
    freqs = {}
    seqcount = 0.
    gapcount = 0
    coverage = []
    for l in open(alignment):
        if l.find('>') > -1:
            continue
        x = l.strip()
        if len(x) < 3:
            continue
        seqcount += 1
        coverage.append( (len(x) - x.count('-'))/float(len(x)))
        for i in range(len(x)):
            try:
                freqs[i][x[i]] += 1 
            except:
                try:
                    freqs[i][x[i]] = 1 
                except:
                    freqs[i] = {} 
                    freqs[i][x[i]] = 1 
            if x[i] == '-':
                gapcount += 1
    b['-'] = gapcount/(seqcount * len(freqs.keys()))
    entropy = []
    for i in sorted(freqs.keys()):
        q = []
        for l in one2number:
            try:
                q.append(np.log( freqs[i][l]/(b[l] * seqcount)))
                q.append(freqs[i][l]/(b[l] * seqcount) * np.log( freqs[i][l]/(b[l] * seqcount)))
                entropy.append(freqs[i][l]/(b[l] * seqcount) * np.log( freqs[i][l]/(b[l] * seqcount)))
            except Exception as e:
                q.append(np.log( 0.1/(b[l] * seqcount)))
                q.append(0)
                entropy.append(0)
        pssm[i+1] = q
    return (pssm, np.mean(entropy), [np.min(coverage), np.max(coverage), np.mean(coverage), np.median(coverage)])

def parseStats(f):
    stats = []
    ff = open(f).readlines()
    if len(ff) != 6:
        sys.stderr.write(f + ' has incorrect format!\n')
        return [-1, -1, -1, -1, -1, -1]
    for l in ff:
        stats.append(float(l.split()[-1]))
    return stats

def parseContacts(f):
    contacts = set()
    for l in open(f):
        x = l.split()
        if len(x) != 3:
            sys.stderr.write('Incorrect format for ' + f)
            sys.exit(1)
        if float(x[-1]) < 8:
            contacts.add( (int(x[0]), int(x[1])) )
    return contacts
files = sys.argv[1:]

selected = set()
contacts = {}
X = []
Y = []
maxres = -1
outfile = files[7]
if os.path.exists(outfile):
    pass

sys.stderr.write('Doing ' + outfile + '\n')
accessibility = parseNetSurfP(files[3])
SSdict = parsePSIPRED(files[4])
stats  = parseStats(files[5])
pssm = parsePSSM(files[6])
entropy = pssm[1]
coverage = pssm[2]
pssm = pssm[0]

selected = set()
for index in range(3):
    contacts[index] = {}
    d = files[index]
    r = []
    if not os.path.exists(d):
        sys.stderr.write(d + ' does not exist!\n')
        continue
    infile = open(d).readlines()
    for m in infile:
        if d.find('gdca') > -1:
            x = m.split()
            c = 2
        elif d.find('.plm') > -1:
            x = m.split(',')
            if len(x) != 3:
                print d + ' has wrong format!'
                sys.exit(1)
        else:
            x = m.split()
            if len(x) < 3 or x[2] != '0' or x[3] != '8':
                continue
            c = -1
        if len(x) < 3:
            continue
        aa1 = int(x[0])
        aa2 = int(x[1])
        if aa1 > maxres:
            maxres = aa1
        if aa2 > maxres:
            maxres = aa2    
        if x[c].find('nan') > -1:
            score = -3
        else:
            score = float(x[c])
        contacts[index][(aa1, aa2)] = score
        contacts[index][(aa2, aa1)] = score
        selected.add((aa1,aa2))

clist = []
for c in contacts[0].keys():
    q = [ c ]
    for i in contacts.keys():
        try:
            q.append( contacts[i][c] )
        except:
            q.append( -3 )
    clist.append(q)

selected2 = set()
for i in contacts.keys():
    clist.sort(key = lambda x: -x[i+1])
    counter = -1
    c = 0
    while counter < maxres:
        j = clist[c]
        selected2.add(j[0])
        c+=1
        if abs(j[0][0] - j[0][1]) > 4:
            counter += 1

maxscores = []
meantop = []
stdtop = []
for index in range(3):
        maxscores.append(max(contacts[index].values()))
        q = []
        for s in list(selected2):
                try:
                        q.append(contacts[index][s])
                except:
                        pass
        meantop.append(np.mean(q))
        stdtop.append(np.std(q))

selected = list(selected)

selected.sort()
lastseeny = -1

X = []
Y = []
sys.stderr.write('Reading in data\n')
sys.stderr.flush()
count = 0
allcount = len(selected)
start = time.time()

for s in selected:
    count += 1
    if count % 100 == 0:
        sys.stderr.write('\rProgress: [' + '#' * (80*count/allcount) + ' ' * (80*(allcount-count)/allcount) + ']')
        now = time.time()
        sys.stderr.write('Time remaining: {:7.1f}s'.format( (allcount-count) * (now-start)/count ) )
    q = []
    q.append(abs(s[0]-s[1]))

    for ss in stats:
        q.append(ss)
    q.append(entropy)
    q.extend(coverage)
        
    q.append(maxscores[0])
    q.append(maxscores[1])
    q.append(maxscores[2])

    for i in range(-5, 6):
        for j in range(-5, 6):
            for index in range(3):
                try:
                    q.append(contacts[index][(s[0]+i, s[1]+j)])
                    q.append((contacts[index][(s[0]+i, s[1]+j)] - meantop[index])/stdtop[index])
                except:
                    q.append(0)
                    q.append(0)

        for i in range(-4, 5):
                try:
                        q.extend(SSdict[s[0]+i] )
                except:
                        q.extend([0,0,0])

        for i in range(-4, 5):
                try:
                        q.extend(SSdict[s[1]+i] )
                except:
                        q.extend([0,0,0])
    
    for i in range(-4, 5):
        try:
            q.extend(accessibility[s[0]+i] )
        except:
            q.extend([0,0,0,0,0])
    
    for i in range(-4, 5):
        try:
            q.extend(accessibility[s[1]+i] )
        except:
            q.extend([0,0,0,0,0])

    q.extend(pssm[s[0]])
    q.extend(pssm[s[1]])
    
    X.append(q)

sys.stderr.write('\n')
sys.stderr.flush()

def predict(dir, X):
    if not os.path.exists(dir + '/tree.list'):
        sys.stderr.write('Directory {:s} does not contain proper random forest!\n'.format(dir))
        sys.exit(0)
    trees = open(dir + '/tree.list').read().strip().split('\n')
    random.shuffle(trees)
    trees = trees[:int(len(trees)*treefraction)]
    predictions = np.zeros(len(X))
    start = time.time()
    count = 0
    allcount = len(trees)
    ccc = 0.
    for t in trees:
        try:
            t = joblib.load(dir + '/' + t.split('/')[-1])
        except:
            continue
        for i in range(len(X)):
            rrr = predict_tree(t, X[i])
            predictions[i] += rrr
        count += 1
        sys.stderr.write('\rProgress: [' + '#' * (80*count/allcount) + ' ' * (80- (80*count)/allcount) + ']')
        now = time.time()
        sys.stderr.write(' {:6.1f} trees/s. Time remaining: {:7.1f}s'.format( count/(now-start), min( (allcount-count) * (now-start)/count, maxtime-(now-start)) ))
        if now - start > maxtime:
            break
    predictions = predictions/count
    return predictions

def predict_tree(tree, q):
    v = [0,0]
    i = 0
    while i >= 0:
        if q[tree[0][i]] <= tree[1][i]:
            j = tree[2][i]
        else:
            j = tree[3][i]
        if j < 0:
            v = tree[4][i][0]
        i = j
    return float(v[1])/sum(v)

# first layer
sys.stderr.write('\nPredicting base layer:\n')
p = predict(forestlocation + '/tlayer0/', X)
of = open(outfile + '.l0', 'w')
previouslayer = {} 

for t in range(len(p)):
    of.write('{:d} {:d} {:7.5f}\n'.format(selected[t][0], selected[t][1], p[t]))
    try:
        previouslayer[selected[t][0]][selected[t][1]] = p[t]
    except:
        previouslayer[selected[t][0]] = {}
        previouslayer[selected[t][0]][selected[t][1]] = p[t]
of.close()

Xp = X
Yp = selected
for layer in range(1,5):
    X = []
    sys.stderr.write('\nPredicting convolution layer {:d}:\n'.format(layer))
    for p in range(len(Xp)):
        y = Yp[p]
        q = list(Xp[p])
        for i in range(-5,6):
            for j in range(-5,6):
                try:
                    q.append(previouslayer[y[0]+i][y[1] + j])
                except:     
                    q.append(-3)
        X.append(q)

    p = predict(forestlocation + '/tlayer{:d}/'.format(layer), X)

    previouslayer = {}
    of = open(outfile + '.l{:d}'.format(layer), 'w')
    for t in range(len(p)):
        of.write('{:d} {:d} {:7.5f}\n'.format(Yp[t][0], Yp[t][1], p[t]))
        try:
            previouslayer[Yp[t][0]][Yp[t][1]] = p[t]
        except:
            previouslayer[Yp[t][0]] = {}
            previouslayer[Yp[t][0]][Yp[t][1]] = p[t]
    of.close()

sys.stderr.write('\n\nSuccesfully completed in {:7.1f} seconds\n'.format(time.time() - firststart))
