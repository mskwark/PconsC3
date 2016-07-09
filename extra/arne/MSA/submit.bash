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

snic=snic2016-10-22

for i in "$@"
do
    if [ -s $i.JH0.001.trimmed ]
    then
	i=$i.JH0.001.trimmed
    fi
    j=`basename $i .trimmed`
    m=`echo $j | sed "s/\.fasta.*//" | sed "s/\.seq.*//"  | sed "s/\.fa.*//"`
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
    
    minitime="01:00:00"
    shorttime="12:00:00"
    longtime="24:00:00"
    mem="60GB"
    
    if [ ! -s $j.trimmed ] && [ ! -s $i.JH0.001.trimmed  ]&& [ ! -s $i.HH0.001.trimmed  ]
    then
	if [ -s $j.a3m ]
	then
	    ~/git/PconsC3/extra/arne/MSA/a3mToTrimmed.py $j.a3m > $j.trimmed
	elif [ -s $j.JH0.001.a3m ]
	then
	    ~/git/PconsC3/extra/arne/MSA/a3mToTrimmed.py $j.JH0.001.a3m > $j.JH0.001.trimmed 

	else
	     srun -A $snic --time=$longtime -n 1 -c 6 ~/git/PconsC3/extra/arne/MSA/runjackhmmer.py $i > $j-runjackhmmer.out &
	     srun -A $snic --time=$longtime -n 1 -c 6 ~/git/PconsC3/extra/arne/MSA/runhhblits.py $i > $j-runhhblits.out &
	     srun -A $snic --time=$longtime -n 1 -c 6 ~/git/PconsC3/extra/arne/MSA/runjackhmmer.py -e 1  $i > $j-runjackhmmer-e1.out &
	     srun -A $snic --time=$longtime -n 1 -c 6 ~/git/PconsC3/extra/arne/MSA/runhhblits.py -e 1 $i > $j-runhhblits-e1.out &
	     srun -A $snic --time=$longtime -n 1 -c 6 ~/git/PconsC3/extra/arne/MSA/runjackhmmer.py -e 1.e-10  $i > $j-runjackhmmer-e-10.out &
	     srun -A $snic --time=$longtime -n 1 -c 6 ~/git/PconsC3/extra/arne/MSA/runhhblits.py -e 1.e.10 $i > $j-runhhblits-e-10.out &
	     srun -A $snic --time=$longtime -n 1 -c 6 ~/git/PconsC3/extra/arne/MSA/runjackhmmer.py -e 1.e-40  $i > $j-runjackhmmer-e-40.out &
	     srun -A $snic --time=$longtime -n 1 -c 6 ~/git/PconsC3/extra/arne/MSA/runhhblits.py -e 1.e.40 $i > $j-runhhblits-e-40.out &
	fi
    else
# #    if [ !  -s $k.out ]
# #    then
# #        srun -A $snic --time=04:00:00 -n 1 -c 8 ~/git/PconsC3/extra/arne/MSA/run_PconsC3.bash  $i &> $j.out &
# #    fi
	length=`head -2 $j.trimmed | wc -c `
	
	if [ $length -gt 500 ] 
	then
	    minitime="08:00:00"
	    shorttime="24:00:00"
	    longtime="48:00:00"
	    mem="60GB"
	fi

	if [ ! -s $j.rr ] && [ -s $m.rr ]
	then
	    ln -s $m.rr $j.rr 
	fi
	if [ ! -s $j.fa.rsa ] && [ -s $m.rsa ]
	then
	    ln -s $m.rsa $j.fa.rsa 
	fi

	if [ !  -s $j.rr ] && [ ! -s $m.rr ]
	then
            srun --mem $mem -A $snic --time=$longtime -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/runPhyCMAP.bash $n &> $m-phycmap.out &
  	    echo "srun --mem $mem -A $snic --time=$longtime -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/runPhyCMAP.bash $n &> $m-phycmap.out &"
	fi
	if [ !  -s $n.rsa ] && [ ! -s $i.rsa ]
	then
            srun -A $snic --time=$minitime -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/runnetsurfp.py $n &> $j-rsa.out &
            echo "srun -A $snic --time=$minitime -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/runnetsurfp.py $n &> $j-rsa.out &"
 	    #$HOME/git/PconsC3/extra/arne/MSA/runnetsurfp.py $n &> $j-rsa.out 
 	    #ln -s $m.rsa $i.rsa
	fi
	if [ !  -s $j.ss2 ] && [ -s $j.trimmed  ]
	then
            srun  -A $snic --time=$shorttime -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/addss.pl $i &> $j-adss.out &
            echo "srun -A $snic --time=$shorttime -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/addss.pl $i &> $j-adss.out &"
	fi
	if [ !  -s $j.gdca ] && [ -s $j.trimmed  ]
	then
            srun --mem $mem -A $snic --time=$longtime -n 1 -c 6 $HOME/git/PconsC3/rungdca.py $i &> $j-gdca.out &
            echo "srun -A $snic --time=$longtime -n 1 -c 6 $HOME/git/PconsC3/rungdca.py $i &> $j-gdca.out &"
	fi
	if [ !  -s $j.0.02.plm20 ]
	then
            srun -A $snic --mem $mem --time=$longtime -n 1 -c 6 $HOME/git/PconsC3/runplm.py $i &> $j-plm.out &
            echo "srun -A $snic --time=$longtime -n 1 -c 6 $HOME/git/PconsC3/runplm.py $i &> $j-plm.out &"
	fi
	
	if [ !  -s $j.PconsC3.l5 ] && [ -s $j.trimmed ]
	then
  	    if [ -s $j.gdca ]  && [ -s $j.0.02.plm20 ]   && [ -s $j.rr ]  && [ -s $n.rsa ]  && [ -s  $j.ss2 ] && [ -s  $j.gneff ]  && [ -s $j.trimmed ]
  	    then
		echo " srun -A $snic --time=$minitime -n 1 -c 6 $HOME/git/PconsC3/predict.py $j.gdca $j.0.02.plm20 $j.rr $n.rsa $j.ss2 $j.gneff $j.trimmed $HOME/git/PconsC3/ -1 $j.PconsC3  &> $j-predict.out &"
		srun -A $snic --mem $mem  --time=$minitime -n 1 -c 6 $HOME/git/PconsC3/predict-parallel.py $j.gdca $j.0.02.plm20 $j.rr $n.rsa $j.ss2 $j.gneff $j.trimmed $HOME/git/PconsC3/ -1 $j.PconsC3  &> $j-predict.out &
  	    else 
  		ls -l $j.gdca $j.0.02.plm20 $j.rr $n.rsa $j.ss2 $j.gneff $j.trimmed 
  	    fi
	fi
    fi
    sleep 10
done
