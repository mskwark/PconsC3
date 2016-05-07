#!/bin/bash -x
#SBATCH -A snic2015-10-12 
# We actually start 6 jobs in parallel.
# Probably more efficient than running with 5 thread.
#SBATCH -n 6 
#SBATCH -c 1 
#SBATCH --time=04:00:00 
#SBATCH -J RunJackhmmer 
#SBATCH --output=out/PconsC3.%J.out
#SBATCH --error=err/PconsC3.%J.err 

for i in "$@"
do
    if [ -s $i.JH0.001.trimmed ]
    then
	i=$i.JH0.001.trimmed
    fi
    j=`basename $i .trimmed`
    m=`echo $i | sed "s/\.fasta.*//" | sed "s/\.seq.*//"  | sed "s/\.fa.*//"`
    n=$m.fasta
    if  [ -s $n ]
    then
	echo "Using $n"
    else
	n=$m.fa
	if  [ -s $n ]
	then
	    echo "Using $n"
	else
	    n=$m.fasta
	    if  [ -s $n ]
	    then
		echo "Using $n"
		n=$m.fa
	    else
		n=$m.seq
		if  [ -s $n ]
		then
		    echo "Using $n"
		else
		    echo "ERROR: Not found FASTa file $i $m $n"
#		    exit 0
		fi
	    fi
	fi
    fi
    
    
    if [ !  -s $j.PconsC3.l5 ] && [ -s $j.trimmed ]
    then
  	if [ -s $j.gdca ]  && [ -s $j.0.02.plm20 ]   && [ -s $j.rr ]  && [ -s $n.rsa ]  && [ -s  $j.ss2 ] && [ -s  $j.gneff ]  && [ -s $j.trimmed ]
  	then
	    echo " srun -A snic2015-10-12 --time=$longtime -n 1 -c 6 $HOME/git/PconsC3/predict.py $j.gdca $j.0.02.plm20 $j.rr $n.rsa $j.ss2 $j.gneff $j.trimmed $HOME/git/PconsC3/ -1 $j.PconsC3  &> $j-predict.out &"
	    $HOME/git/PconsC3/predict-parallel.py $j.gdca $j.0.02.plm20 $j.rr $n.rsa $j.ss2 $j.gneff $j.trimmed $HOME/git/PconsC3/ -1 $j.PconsC3  &> $j-predict.out 
  	else 
  	    ls -l $j.gdca $j.0.02.plm20 $j.rr $n.rsa $j.ss2 $j.gneff $j.trimmed 
  	fi
    fi
done
