#!/bin/bash -x

# Basically a verion of run_pconsC3 that starts with an alignment (useful for pairwise predictions)

#make sure everything is in the path
export PATH=$PATH:$HOME/git/PconsC3/
#export PATH=$PATH:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64/bin:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64_patch/bin:/scratch/arne/PconsC3/bin/../dependencies/netsurfp-1.0/bin:/scratch/arne/PconsC3/bin/../dependencies/phycmap.release/bin:/scratch/arne/PconsC3/bin/../dependencies/psipred/bin:/scratch/arne/PconsC3/bin/../dependencies/blast:/scratch/arne/PconsC3/bin/../dependencies/cd-hit-v4.5.4-2011-03-07:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64_patch:/scratch/arne/PconsC3/bin/../dependencies/hmmer-3.1b2-linux-intel-x86_64:/scratch/arne/PconsC3/bin/../dependencies/netsurfp-1.0:/scratch/arne/PconsC3/bin/../dependencies/phycmap.release:/scratch/arne/PconsC3/bin/../dependencies/plmDCA_asymmetric_v2:/scratch/arne/PconsC3/bin/../dependencies/psipred

name=`basename $1 .aln`
dir=`dirname $1`
seqfile=$dir"/"$name

bin="/scratch/arne/contactpreds/chaperonin/bin/"
PconsC3=$HOME/git/PconsC3/

# First sequence needs to be trimmed
$bin/trimsequence.py $1 > $seqfile.trimmed
head -2 $seqfile.trimmed > $seqfile.fasta

# Run svmcon (the best thing I could find)
$bin/predict_map_pair.sh $seqfile.fasta $seqfile.svmcon.rr $seqfile.trimmed.aln

# Add secondary structure (replaces by sspro..)
#$bin/addss.pl $seqfile.fas $seqfile.addss -fas
# Surface are prediction Run Netsurfp
#$bin/runnetsurfp.py $seqfile.fas
#

# Run plmdca
#$PconsC3/runplmdca.py $seqfile.trimmed
# Run gdca
#$PconsC3/rungdca.py $seqfile.trimmed

# 
