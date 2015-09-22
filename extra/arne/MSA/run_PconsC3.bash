#!/bin/bash -x

# A script to run PconsC3 and all necessary pre-calculations
seqfile=$1  # Required: Sequence or multiple sequence (ALN) file 
workdir=$2 # Optional: Working directory

# THis version uses SVMCON (and not phyCMAP) as we have not enabled
# the use of precalculated alignments for PhyCMAP and also due to the
# unavailability of PhyCMAP. Please note that the SVMCON is "hacked"
# so that it uses SCRATCH-1D predictions and SSPRO4 (as this does not
# run on a modern linux box).
#


#
# A few variables to set.
#make sure everything is in the path
export PATH=$PATH:$HOME/git/PconsC3/
export PATH=$PATH:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64/bin:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64_patch/bin:/scratch/arne/PconsC3/bin/../dependencies/netsurfp-1.0/bin:/scratch/arne/PconsC3/bin/../dependencies/phycmap.release/bin:/scratch/arne/PconsC3/bin/../dependencies/psipred/bin:/scratch/arne/PconsC3/bin/../dependencies/blast:/scratch/arne/PconsC3/bin/../dependencies/cd-hit-v4.5.4-2011-03-07:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64:/scratch/arne/PconsC3/bin/../dependencies/hhsuite-2.0.16-linux-x86_64_patch:/scratch/arne/PconsC3/bin/../dependencies/hmmer-3.1b2-linux-intel-x86_64:/scratch/arne/PconsC3/bin/../dependencies/netsurfp-1.0:/scratch/arne/PconsC3/bin/../dependencies/phycmap.release:/scratch/arne/PconsC3/bin/../dependencies/plmDCA_asymmetric_v2:/scratch/arne/PconsC3/bin/../dependencies/psipred

bin="/scratch/arne/contactpreds/chaperonin/bin/"
PconsC3=$HOME/git/PconsC3/
#HHLIB=/scratch/arne/PconsC2-extra/hhsuite-2.0.16-linux-x86_64/lib/hh/
export HHLIB=/usr/local/lib/hh/
# ----------

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





seqname=`basename  $seqfile`
rootname=`echo $seqname | sed -E "s/\..*//"`
cp $seqfile $workdir
cd $workdir

#check if we have a muktiple sequence alignment or a single sequence
numseq=`grep -c \> $seqname`
if [[ $numseq -gt 2 ]]
then
# First sequence needs to be trimmed
    $bin/trimsequence.py $seqname > $rootname.trimmed
    head -2 $rootname.trimmed > $rootname.fasta
    SEQ=$rootname.fasta
    ALN=$rootname.trimmed
else
# We need to run hhblits
    i='hhE0'
    
    if [[ ! -s $rootname.trimmed ]] 
    then 
	$bin/runhhblits.py -c $cpu -name $i -e 1 $seqname
	$bin/a3mToTrimmed.py $seqname.$i.a3m > $rootname.trimmed 
    fi
    SEQ=$seqname
    ALN=$rootname.trimmed
fi

# Now we should have the sequences.

# Predict secondary structure (Could perhaps be replaces by SSPRO but format is not the same)
if [ ! -s $rootname.ss2 ]
then
    $bin/addss.pl $ALN $ALN.addss -fas 
fi
# Surface are prediction Run Netsurfp (could perhaps be replace by SSPRO)
if [ ! -s $SEQ.rsa ]
then
    $bin/runnetsurfp.py $SEQ
fi

# Run svmcon (the best thing I could find)
if [ ! -s $SEQ.rr ]
then
    echo $numseq > $ALN.raw
    grep -v \> $ALN >> $ALN.raw
    $bin/predict_map.sh $SEQ $SEQ.rr $ALN.raw
fi


# Run plmdca
if [ ! -s $ALN.0.02.plm20 ]
then
    $PconsC3/runplm.py $ALN
fi
# Run gdca
if [ ! -s $ALN.gdca ]
then
    $PconsC3/rungdca.py $ALN
fo
# and now run PconsC3
$PconsC3/predict.py $ALN.gdca $ALN.0.02.plm20 $SEQ.rr $SEQ.rsa $rootname.ss2 $ALN.gneff $ALN $seqfile.PconsC3
cp *l5 $curddir/

# 
cd $currdir
