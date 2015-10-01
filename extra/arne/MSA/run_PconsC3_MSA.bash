#!/bin/bash -x

# A script to run PconsC3 and all necessary pre-calculations
seqfile=$1  # Required: Sequence file
alnfile=$2  # Optional: Alignment file in a3m format
workdir=$3 # Optional: Working directory

# Basically a verion of run_pconsC3 that starts with an alignment (useful for pairwise predictions)

#make sure everything is in the path
export PATH=$PATH:$HOME/git/PconsC3/
#export PATH=$PATH:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64/bin:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64_patch/bin:/scratch/arne/PconsC3/bin/../dependencies/netsurfp-1.0/bin:/scratch/arne/PconsC3/bin/../dependencies/phycmap.release/bin:/scratch/arne/PconsC3/bin/../dependencies/psipred/bin:/scratch/arne/PconsC3/bin/../dependencies/blast:/scratch/arne/PconsC3/bin/../dependencies/cd-hit-v4.5.4-2011-03-07:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64_patch:/scratch/arne/PconsC3/bin/../dependencies/hmmer-3.1b2-linux-intel-x86_64:/scratch/arne/PconsC3/bin/../dependencies/netsurfp-1.0:/scratch/arne/PconsC3/bin/../dependencies/phycmap.release:/scratch/arne/PconsC3/bin/../dependencies/plmDCA_asymmetric_v2:/scratch/arne/PconsC3/bin/../dependencies/psipred

bin="/scratch/arne/contactpreds/chaperonin/bin/"
PconsC3=$HOME/git/PconsC3/


currdir=`pwd -P`
if [[ $workdir == '' ]];then 
    workdir=`echo '$dir=int(rand(100000));$dir=".PconsC3.$dir.$seqfile";if(-d $dir||-e $dir){ }else{print $dir}' | perl - $pdbid` 
    workdir=`pwd`/$workdir
fi
if [[ ! -d $workdir ]] 
then
    mkdir -p $workdir
    if [ $? -ne 0 ];then echo "ERROR cannot make $workdir" ; exit $? ; fi
fi
cp $seqfile $workdir

seqfile=`basename $1 .aln`
#seqfile=$dir"/"$name
homedir=`pwd -P`
cd $dir


# First sequence needs to be trimmed
$bin/trimsequence.py $1 > $seqfile.trimmed
head -2 $seqfile.trimmed > $seqfile.fasta

# Add secondary structure (Could perhaps be replaces by SSPRO but format is not the same)
$bin/addss.pl $seqfile.fas $seqfile.addss -fas
# Surface are prediction Run Netsurfp (could perhaps be replace by SSPRO)
$bin/runnetsurfp.py $seqfile.fas
#

# Run svmcon (the best thing I could find)
$bin/predict_map.sh $seqfile.fasta $seqfile.svmcon.rr $seqfile.trimmed.aln


# Run plmdca
$PconsC3/runplmdca.py $seqfile.trimmed
# Run gdca
$PconsC3/rungdca.py $seqfile.trimmed

# and now run PconsC3
$PconsC3/predict.py $seqfile.gdca $seqfile.0.02.plm20 $seqfile.rr $seqfile.rsa $seqfile.ss2 $seqfile.gneff $seqfile.trimmed $seqfile.PconsC3


# 
cd $homedir
