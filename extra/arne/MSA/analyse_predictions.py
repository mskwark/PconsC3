#!/usr/bin/env python
from __future__ import division
import sys, os, re, string
import argparse
from math import *

# on UPPMAX only
sys.path.append('/sw/apps/bioinfo/biopython/1.59/tintin/lib/python')

from Bio import pairwise2

import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#from mpl_toolkits.axes_grid1 import make_axes_locatable

from os.path import expanduser
home = expanduser("~")
sys.path.append(home + '/bioinfo-toolbox/parsing')
sys.path.append(home + '/git/bioinfo-toolbox/parsing')

import parse_contacts
import parse_psipred
import parse_fasta
import parse_iupred
import parse_pdb


def s_score(d, d0):
    return 1/(1+pow(d/d0, 2))

def s_score_vec(d, d0):
    f = np.vectorize(s_score, otypes=[np.float])
    return f(d, d0)


def get_min_dist(res1, res2):
    
    min_dist = float('inf')

    for atm1 in res1:
        for atm2 in res2:
            diff_vec = atm1 - atm2
            dist = np.sqrt(np.sum(diff_vec * diff_vec))
            if dist < min_dist:
                min_dist = dist

    return min_dist


def get_heavy_contacts(gapped_res_lst):

    seqlen = len(gapped_res_lst)
    dist_mat = np.zeros((seqlen, seqlen), np.float)
    dist_mat.fill(float('inf'))
    
    for i, res1 in enumerate(gapped_res_lst):
        if res1 == '-':
            continue
        for j, res2 in enumerate(gapped_res_lst):
            if res2 == '-':
                continue
            dist_mat[i,j] = get_min_dist(res1[1], res2[1])
    return dist_mat


def get_cb_contacts(gapped_cb_lst):

    seqlen = len(gapped_cb_lst)
    dist_mat = np.zeros((seqlen, seqlen), np.float)
    dist_mat.fill(float('inf'))
    
    #offset = 0
    #first_i = gapped_cb_lst[0].keys()[0]
    #if first_i < 0:
    #    offset = abs(first_i)

    for i, cb1 in enumerate(gapped_cb_lst):
        if str(cb1) == str('-'):
            continue
        for j, cb2 in enumerate(gapped_cb_lst):
            if str(cb2) == str('-'):
                continue
            diff_vec = cb1 - cb2
            #dist_mat[i+offset,j+offset] = np.sqrt(np.sum(diff_vec * diff_vec))
            dist_mat[i,j] = np.sqrt(np.sum(diff_vec * diff_vec))
    return dist_mat



def print_contacts(fasta_filename,score,contacts_x, contacts_y, ref_contact_map, atom_seq_ali, ref_len, factor,disorder):
    #    for num_c in range(min(len(contacts_x), int(ceil(ref_len * factor))) + 1)[1:]:
    TP = 0.0
    FP = 0.0
    disoTP = 0.0
    disoFP = 0.0
    mixTP = 0.0
    mixFP = 0.0
    for i in range(len(contacts_x) ):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali[c_x] == '-':
            continue
        if atom_seq_ali[c_y] == '-':
            continue
        if len(disorder)>i:
            print "DATA: ",fasta_filename,i,scores[i],c_x,c_y,disorder[c_x],disorder[c_y]
        else:
            print "DATA: ",fasta_filename,i,scores[i],c_x,c_y



def get_ppvs(contacts_x, contacts_y, ref_contact_map, atom_seq_ali, ref_len, factor,disorder):

    PPVs = []
    TPs = []
    FPs = []
    disoPPVs = []
    disoTPs = []
    disoFPs = []
    mixPPVs = []
    mixTPs = []
    mixFPs = []
    disocount = 0
    mixcount = 0
    #    for num_c in range(min(len(contacts_x), int(ceil(ref_len * factor))) + 1)[1:]:
    TP = 0.0
    FP = 0.0
    disoTP = 0.0
    disoFP = 0.0
    mixTP = 0.0
    mixFP = 0.0
    for i in range(len(contacts_x) ):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali[c_x] == '-':
            continue
        if atom_seq_ali[c_y] == '-':
            continue
        if (disorder[c_x] > 0.5 and disorder[c_y] > 0.5):
            if (disocount < ref_len * factor):
                disocount+=1
            if ref_contact_map[c_x, c_y] > 0:
                disoTP += 1.0 
            else:
                disoFP += 1.0 
            disoPPVs.append(disoTP / (disoTP + disoFP))
            disoTPs.append(disoTP/disocount)
            disoFPs.append(disoFP/disocount)

        elif (disorder[c_x] > 0.5 or disorder[c_y] > 0.5):
            if (mixcount < ref_len * factor):
                mixcount+=1
            if ref_contact_map[c_x, c_y] > 0:
                mixTP += 1.0 
            else:
                mixFP += 1.0 
            mixPPVs.append(mixTP / (mixTP + mixFP))
            mixTPs.append(mixTP/mixcount)
            mixFPs.append(mixFP/mixcount)
        if (i < ref_len * factor):
            if ref_contact_map[c_x, c_y] > 0:
                TP += 1.0 / (ref_len*factor)
            else:
                FP += 1.0 / (ref_len*factor)
            PPVs.append(TP / (TP + FP))
            TPs.append(TP)
            FPs.append(FP)


    if len(PPVs) == 0:
        PPVs.append(0.0)
    if len(TPs) == 0:
        TPs.append(0.0)
    if len(FPs) == 0:
        FPs.append(0.0)
    if len(mixPPVs) == 0:
        mixPPVs.append(0.0)
    if len(mixTPs) == 0:
        mixTPs.append(0.0)
    if len(mixFPs) == 0:
        mixFPs.append(0.0)
    if len(disoPPVs) == 0:
        disoPPVs.append(0.0)
    if len(disoTPs) == 0:
        disoTPs.append(0.0)
    if len(disoFPs) == 0:
        disoFPs.append(0.0)

    return PPVs, TPs, FPs,mixPPVs, mixTPs, mixFPs,disoPPVs, disoTPs, disoFPs


def get_tp_colors(contacts_x, contacts_y, ref_contact_map, atom_seq_ali):

    tp_colors = []

    for i in range(len(contacts_x)):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali[c_x] == '-':
            #tp_colors.append('green')
            tp_colors.append('red')
            continue
        if atom_seq_ali[c_y] == '-':
            #tp_colors.append('green')
            tp_colors.append('red')
            continue
        if ref_contact_map[c_x, c_y] > 0:
            tp_colors.append('blue')
        else:
            tp_colors.append('red')

    return tp_colors
 

def get_colors(contacts_np, ref_contact_map=[], atom_seq_ali=[], th=0.5):

    N = contacts_np.shape[0]
    img = np.ones((N,N,4))

    for i in xrange(N):
        for j in xrange(N):
            if j < (i+5):
                continue
            sc = contacts_np[i,j]
            if len(ref_contact_map) > 0:
                assert N == ref_contact_map.shape[0]
                # FN
                if sc <= th and ref_contact_map[i,j] < 8:
                    #img[i,j] = [0.5,0.5,1,1]
                    img[i,j] = [0.5,0.5,0.5,1]
                    img[j,i] = [0.5,0.5,0.5,1]
                # TP
                elif sc > th and ref_contact_map[i,j] < 8:
                    img[i,j] = [0,1,0,1]
                    img[j,i] = [0,1,0,1]
                    #img[j,i] = [1-sc,1-sc,1-sc,1]
                # FP
                #elif contacts_np[i,j] > th and ref_contact_map[i,j] >= 8:
                elif sc > th and ref_contact_map[i,j] >= 12:
                    img[i,j] = [1,0,0,1]
                    img[j,i] = [1,0,0,1]
                    #img[j,i] = [1-sc,1-sc,1-sc,1]
                # grey zone between 8 and 12 Angstroem
                elif sc > th and (ref_contact_map[i,j] < 12 or ref_contact_map[i,j] >= 8):
                    val = (ref_contact_map[i,j] - 8)/(12 - 8)
                    img[i,j] = [0.5+val/2,1-val/2,0,1]
                    img[j,i] = [0.5+val/2,1-val/2,0,1]
                    #img[j,i] = [1-sc,1-sc,1-sc,1]
            else:
                if sc > th:
                    img[i,j] = [0.5-sc/2,0.5-sc/2,1,1]
                    img[j,i] = [0.5-sc/2,0.5-sc/2,1,1]

    return img


def get_seqlen(filename):
    alifile = open(filename, 'r')
    l = alifile.readline()
    if l.startswith('>'):
        L = len(alifile.readline().strip())
    else:
        L = len(l.strip())
    alifile.close()
    return L


def get_ali_coverage(filename):
    L = get_seqlen(filename)
    N = 0
    alifile = open(filename, 'r')
    coverage = dict.fromkeys(range(L), 0)
    for line in alifile:
        # ignore headers
        if line.startswith('>'):
            continue
        # remove possible inserts
        line = line.translate(None,string.ascii_lowercase)
        pos_lst = [m.start() for m in re.finditer('-', line)]
        for pos in pos_lst:
            coverage[pos] += 1
        N += 1
    alifile.close()
    #coverage_lst = [1-(coverage[i]/float(N)) for i in range(L)]
    coverage_lst = [N - coverage[i] for i in range(L)]
    return coverage_lst


def get_meff_coverage(meff_file):
    meff_lst = []
    with open(meff_file) as f:
        for l in f:
            if 'MeffPerPos' in l:
                meff_lst = l.split()[-1].strip('[]').split(',')
                meff_lst = map(float, meff_lst)
    return meff_lst

            
#    if outfilename:
#        out
#    else:
#        outfile=%S

    #    ax = fig.add_subplot(111)#, aspect='auto')
    #    ax.set_adjustable('box-forced')
    #    ax.tick_params(direction='out', right='off', top='off')
    #    ax.set_xlim([-unit,ref_len])
    #    ax.set_ylim([-unit,ref_len])

def contactanalysis(fasta_filename, c_filename, factor=1.0, cutoff=9999.99, th=-1, c2_filename='', psipred_horiz_fname='', psipred_vert_fname='',iupred_fname='', pdb_filename='', is_heavy=False, chain='', sep=',', outfilename='', ali_filename='',  meff_filename='', name='', start=0, end=-1):  
    #acc = c_filename.split('.')[0]
    #acc = fasta_filename.split('.')[0][:4]
    if name == '': 
        acc = '.'.join(os.path.basename(fasta_filename).split('.')[:-1])
    else:
        acc = name

    ### get sequence
    seq = parse_fasta.read_fasta(open(fasta_filename, 'r')).values()[0][0]
    ref_len = len(seq)

    ### trim sequence according to given positions
    ### default: take full sequence
    if end == -1:
        end = ref_len
    seq = seq[start:end]
    ref_len = len(seq)
    unit = (ref_len/50.0)

    if ali_filename:
        coverage_lst = get_ali_coverage(ali_filename)
        max_cover = max(coverage_lst)
    elif meff_filename:
        coverage_lst = get_meff_coverage(meff_filename)
        max_cover = max(coverage_lst)
    else:
        max_cover = 0
    average_disorder=0.
    fraction_disorder=0.
    if iupred_fname:
        disorder = parse_iupred.pred(open(iupred_fname, 'r'))
    else:
        disorder=np.zeros(ref_len)
    average_disorder = np.sum(disorder)/ref_len
    fraction_disorder = 0.0
    for i in disorder:
        if (i>0.5):
            fraction_disorder +=1/ref_len
            
    ### get top "factor" * "ref_len" predicted contacts
    contacts = parse_contacts.parse(open(c_filename, 'r'), sep,1)
    contacts_np = parse_contacts.get_numpy_cmap(contacts)
    contacts_np = contacts_np[start:end,start:end]

    contacts_x = []
    contacts_y = []
    scores = []
    mixscores = []
    disoscores = []
    tooclose = []    
    contact_dict = {}

    count = 0
    mixcount = 0
    disocount = 0
    highscore = 0
    numbins=20
    sum=0.0
    disosum=0.0
    mixsum=0.0
    average=0.0
    mixaverage=0.0
    disoaverage=0.0
    histo=np.zeros(numbins)
    disotop=0
    doubletop=0
    mixcount=0
    mixtop=0


    # We actually divide the analysis into three groups (ordered,disordered and mixed)
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
        
        pos_diff = abs(c_x - c_y)
        too_close = pos_diff < 5
        if not too_close:
            if score > cutoff:
                contacts_x.append(c_x - start)
                contacts_y.append(c_y - start)
                if (disorder[c_x] > 0.5 and disorder[c_y] > 0.5):
                    disocount += 1
                    disoscores.append(score)
                    if (disocount <= ref_len * factor):
                        disosum += score
                        disoaverage=disosum/disocount
                elif (disorder[c_x] > 0.5 or disorder[c_y] > 0.5):
                    mixcount += 1
                    mixscores.append(score)
                    if (mixcount <= ref_len * factor):
                        mixsum += score
                        mixaverage=mixsum/mixcount
                count += 1
                scores.append(score)
                if (count <= ref_len * factor):
                    sum += score
                    average=sum/count


        else:
            tooclose.append(score)
           
                
    line="Highs: %.1f (%.1f%%) (%.1f%%)\t average:  %.2f (%.2f) (%.2f)\t Meff: %.0f\t Diso: %.1f%% \t" % (count/ref_len,100*mixcount/count,100*disocount/count,average,mixaverage,disoaverage,max_cover,100*fraction_disorder)
    fig = plt.figure(figsize=(8, 8), dpi=96, facecolor='w')
    plt.hist((tooclose,scores), numbins,range=(0,1), histtype='bar',
             normed=(numbins,numbins), alpha=0.75,
             label=['Too_Close','Contacts'])
    plt.xlabel('Score')
    plt.ylabel('Normalized count')
    fig.suptitle('%s\n%s\n' %  (c_filename,line))

    
    ### Calculate reference contacts in the background if given
    if pdb_filename:
        chain='*'
        # We try to get all chains...
        res_lst = parse_pdb.get_coordinates(open(pdb_filename, 'r'), chain)
#        cb_lst = parse_pdb.get_ca_coordinates(open(pdb_filename, 'r'), chain)
        cb_lst = parse_pdb.get_cb_coordinates(open(pdb_filename, 'r'), chain)
        atom_seq = parse_pdb.get_atom_seq(open(pdb_filename, 'r'), chain)
                
#        align = pairwise2.align.globalms(atom_seq, seq, 2, -1, -0.5, -0.1)
        align = pairwise2.align.globalms(atom_seq, seq, 2, -1, -10.5, -10.1)

        atom_seq_ali = align[-1][0]
        seq_ali = align[-1][1]

        print atom_seq_ali
        print seq_ali

        j = 0
        gapped_res_lst = []
        gapped_cb_lst = []

        for i in xrange(len(atom_seq_ali)):
#            print i,j
            if atom_seq_ali[i] == '-':
                gapped_res_lst.append('-')
                gapped_cb_lst.append('-')
            elif seq_ali[i] == '-':
                j += 1
                continue
            else:
                gapped_res_lst.append(res_lst[j])
                gapped_cb_lst.append(cb_lst[j])
                j += 1

        if is_heavy:
            dist_mat = get_heavy_contacts(gapped_res_lst)
            heavy_cutoff = 5
            ref_contact_map = dist_mat < heavy_cutoff
            ref_contacts = np.where(dist_mat < heavy_cutoff)
        else:
            dist_mat = get_cb_contacts(gapped_cb_lst)
            cb_cutoff = 8
            ref_contact_map = dist_mat < cb_cutoff
            ref_contacts = np.where(dist_mat < cb_cutoff)
            #ref_contacts = np.where(np.ma.array(dist_mat, mask=np.tri(dist_mat.shape[0]), fill_value=float("inf")) < cb_cutoff)
        
        ref_contacts_x = ref_contacts[0]
        ref_contacts_y = ref_contacts[1]

        PPVs, TPs, FPs,mixPPVs, mixTPs, mixFPs,disoPPVs, disoTPs, disoFPs = get_ppvs(contacts_x, contacts_y, ref_contact_map, atom_seq_ali, ref_len, factor,disorder)


           
#        ppv='PPV: %.2f %.2f %.2f' % (float(PPVs[-1]), float(TPs[-1]), float(FPs[-1]))
        ppv='PPV: %.2f (%d)\t%.2f (%d)\t%.2f (%d) ' % (float(PPVs[-1]),len(PPVs), float(mixPPVs[-1]),len(mixPPVs), float(disoPPVs[-1]),len(disoPPVs))
        print "STATs: %s \t %s\t%s\t" % (fasta_filename,line,ppv)
    else:
        print "STATs: %s \t %s\t" % (fasta_filename,line)


    # We should print statistics for each residue..

    

    if outfilename:
        if outfilename.endswith('.pdf'):
            pp = PdfPages(outfilename+"_statistics.pdf")
            pp.savefig(fig)
            pp.close()
        elif outfilename.endswith('.png'):
            plt.savefig(outfilename+"_statistics.png")
        else:
            pp = PdfPages('%s_statistics.pdf' % outfilename)
            pp.savefig(fig)
            pp.close()
    else:
        pp = PdfPages('%s_statistics.pdf' % c_filename)
        pp.savefig(fig)
        pp.close()
        
    
if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Plot protein residue contact maps.')
    p.add_argument('fasta_file')#, required=True)
    p.add_argument('contact_file')#, required=True)
    p.add_argument('-o', '--outfile', default='')
    p.add_argument('-f', '--factor', default=1.0, type=float)
    p.add_argument('-c', '--cutoff', default=0.3, type=float)
    p.add_argument('-t', '--threshold', default=-1, type=float)
    p.add_argument('--c2', default='')
    p.add_argument('--psipred_horiz', default='')
    p.add_argument('--psipred_vert', default='')
    p.add_argument('--iupred', default='')
    p.add_argument('--pdb', default='')
    p.add_argument('--heavy', action='store_true')
    p.add_argument('--chain', default='')
    p.add_argument('--alignment', default='')
    p.add_argument('--meff', default='')
    p.add_argument('--name', default='')
    p.add_argument('--start', default=0, type=int)
    p.add_argument('--end', default=-1, type=int)

    args = vars(p.parse_args(sys.argv[1:]))

    fasta_filename = args['fasta_file']
    c_filename = args['contact_file']
    psipred_filename = args['psipred_horiz']
    iupred_filename = args['iupred']

    # guessing separator of constraint file
    line = open(c_filename,'r').readline()
    if len(line.split(',')) != 1:
        sep = ','
    elif len(line.split(' ')) != 1:
        sep = ' '
    else:
        sep = '\t'
    

    contactanalysis(args['fasta_file'], args['contact_file'], factor=args['factor'], cutoff=args['cutoff'], th=args['threshold'], c2_filename=args['c2'], psipred_horiz_fname=args['psipred_horiz'], psipred_vert_fname=args['psipred_vert'], iupred_fname=args['iupred'], pdb_filename=args['pdb'], is_heavy=args['heavy'], chain=args['chain'], sep=sep, outfilename=args['outfile'], ali_filename=args['alignment'], meff_filename=args['meff'], name=args['name'], start=args['start'], end=args['end'])
