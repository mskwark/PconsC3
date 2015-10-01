#!/bin/bash -x

for i in */*aln 
do 
    echo $i 
    j=`basename $i` 
    k=`dirname $i` 
    cd $k 
    ../bin/trimsequence.py $j > $k.trimmed
    cp $k.trimmed $k.seq.a3m
    head -2 $k.trimmed > $k.seq
    if [[ ! -e $k.gdca ]]
    then
	 /proj/bioinfo/software/PconsC3/rungdca.py $k.trimmed
    fi
    if [[ ! -e $k.0.02.plm20 ]]
    then
	srun -A snic2015-10-12 --time=24:00:00 -n 1 -c 1 /proj/bioinfo/software/PconsC3/runplm.py $k.trimmed &> runplm.out  &
    fi
    if [[ ! -e $k.rr ]]
    then
	srun -A snic2015-10-12 --time=24:00:00 -n 1 -c 1 /proj/bioinfo/users/x_arnel/contactpreds/chaperonin/bin/runPhyCMAP-msa.bash $k.seq $k.seq.a3m -cpu 1 &> runphycmap-msa.out &
    fi
    cd .. 
done
