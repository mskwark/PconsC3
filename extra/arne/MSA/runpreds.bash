#!/bin/bash -x

for i in */*aln 
do 
    echo $i 
    j=`basename $i` 
    k=`dirname $i` 
    cd $k 
    ../bin/trimsequence.py $j > $k.trimmed
    head -2 $k.trimmed > $k.seq
    if [[ ! -e $k.gdca ]]
    then
	/scratch/arne/git/PconsC3/rungdca.py $k.trimmed
    fi
    if [[ ! -e $k.0.02.plm20 ]]
    then
	/scratch/arne/git/PconsC3/runplm.py $k.trimmed 
    fi
    cd .. 
done
