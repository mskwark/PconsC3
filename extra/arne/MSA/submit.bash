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
    sleep 1
    k=`basename $i .fa`
    l=`basename $k .fasta`
    j=`basename $l .trimmed`
    m=`echo $j | sed "s/.fasta.*//"`
    n=$m.fasta
    if  [ -e $n ]
    then
	echo "Using $n"
    else
	m=`echo $j | sed "s/.fa.*//"`
	n=$m.fa
	if  [ -e $n ]
	then
	    echo "Using $n"
	else
	    echo "ERROR: Not found FASTa file $i $m $n"
	    exit 0
	fi
    fi
    
    shorttime="04:00:00"
    longtime="24:00:00"
    mem="120GB"

    if [ ! -e $j.trimmed ]
    then
	if [ -e $j.a3m ] 
	then
	    ~/git/PconsC3/extra/arne/MSA/a3mToTrimmed.py $j.a3m > $j.trimmed
	else
	     srun -A snic2015-10-12 --time=$shorttime -n 1 -c 6 ~/git/PconsC3/extra/arne/MSA/runjackhmmer.py $i > $j-runjackhmmer.out &
	fi
    else
# #    if [ !  -e $k.out ]
# #    then
# #        srun -A snic2015-10-12 --time=04:00:00 -n 1 -c 8 ~/git/PconsC3/extra/arne/MSA/run_PconsC3.bash  $i &> $j.out &
# #    fi
	length=`head -2 $j.trimmed | wc -c `
	
	if [ $length -gt 500 ] 
	then
	    shorttime="06:00:00"
	    longtime="24:00:00"
	    mem="120GB"
	else
	    shorttime="01:00:00"
	    longtime="06:00:00"
	    mem="64GB"
	fi

	if [ !  -e $j.rr ] && [ ! -e $m.rr ]
	then
            srun --mem 120GB -A snic2015-10-12 --time=$longtime -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/runPhyCMAP.bash $n &> $m-phycmap.out &
  	    echo "srun --mem $mem -A snic2015-10-12 --time=$longtime -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/runPhyCMAP.bash $n &> $m-phycmap.out &"
	fi
	if [ !  -e $m.rsa ]
	then
            srun -A snic2015-10-12 --time=$shorttime -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/runnetsurfp.py $m &> $j-rsa.out &
            echo "srun -A snic2015-10-12 --time=00:15:00 -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/runnetsurfp.py $m &> $j-rsa.out &"
 	    # $HOME/git/PconsC3/extra/arne/MSA/runnetsurfp.py $m &> $j-rsa.out 
 	    ln -s $m.rsa $j.fa.rsa
	fi
	if [ !  -e $j.ss2 ] && [ -e $j.trimmed  ]
	then
            srun --mem $mem -A snic2015-10-12 --time=$shorttime -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/addss.pl $i &> $j-adss.out &
            echo "srun -A snic2015-10-12 --time=$shorttime -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/addss.pl $i &> $j-adss.out &"
	fi
	if [ !  -e $j.gdca ] && [ -e $j.trimmed  ]
	then
            srun --mem $mem -A snic2015-10-12 --time=$shorttime -n 1 -c 6 $HOME/git/PconsC3/rungdca.py $i &> $j-gdca.out &
            echo "srun -A snic2015-10-12 --time=$shorttime -n 1 -c 6 $HOME/git/PconsC3/rungdca.py $i &> $j-gdca.out &"
	fi
	if [ !  -e $j.0.02.plm20 ]
	then
            srun -A snic2015-10-12 --mem $mem --time=$shorttime -n 1 -c 6 $HOME/git/PconsC3/runplm.py $i &> $j-plm.out &
            echo "srun -A snic2015-10-12 --time=$longtime -n 1 -c 6 $HOME/git/PconsC3/runplm.py $i &> $j-plm.out &"
	fi
	if [ ! -e $j.rr ] && [ -e $m.rr ]
	then
	    ln -s $m.rr $j.rr 
	fi
	if [ ! -e $j.fa.rsa ] && [ -e $m.rsa ]
	then
	    ln -s $m.rsa $j.fa.rsa 
	fi
	
	if [ !  -e $j.PconsC3.l5 ] && [ -e $j.trimmed ]
	then
  	    if [ -e $j.gdca ]  && [ -e $j.0.02.plm20 ]   && [ -e $j.rr ]  && [ -e $j.fa.rsa ]  && [ -e  $j.ss2 ] && [ -e  $j.gneff ]  && [ -e $j.trimmed ]
  	    then
		echo " srun -A snic2015-10-12 --time=$longtime -n 1 -c 6 $HOME/git/PconsC3/predict.py $j.gdca $j.0.02.plm20 $j.rr $j.fa.rsa $j.ss2 $j.gneff $j.trimmed $HOME/git/PconsC3/ -1 $j.PconsC3  &> $j-predict.out &"
		srun -A snic2015-10-12 --mem $mem  --time=$longtime -n 1 -c 6 $HOME/git/PconsC3/predict.py $j.gdca $j.0.02.plm20 $j.rr $j.fa.rsa $j.ss2 $j.gneff $j.trimmed $HOME/git/PconsC3/ -1 $j.PconsC3  &> $j-predict.out &
  	    else 
  		ls -l $j.gdca $j.0.02.plm20 $j.rr $j.fa.rsa $j.ss2 $j.gneff $j.trimmed 
  	    fi
	fi
    fi
done
