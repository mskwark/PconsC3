#!/bin/bash                                                                                                       
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
    k=`basename $i .fa`
    l=`basename $k .fasta`
    j=`basename $l .trimmed`
#    if [ !  -e $k.out ]
#    then
#        #srun -A snic2015-10-12 --time=04:00:00 -n 1 -c 8 ~/git/PconsC3/extra/arne/MSA/run_PconsC3.bash  $i &> $j.out &
#    fi
     if [ !  -e $j.rr ]
     then
         #srun --mem 96GB -A snic2015-10-12 --time=24:00:00 -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/runPhyCMAP.bash $j.fa &> $j-phycmap.out &
	 echo "#srun --mem 96GB -A snic2015-10-12 --time=24:00:00 -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/runPhyCMAP.bash $j.fa &> $j-phycmap.out &"
     fi
     if [ !  -e $j.fa.rsa ]
     then
         #srun -A snic2015-10-12 --time=01:30:00 -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/runnetsurfp.py $j.fa &> $j-rsa.out &
         echo "#srun -A snic2015-10-12 --time=01:30:00 -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/runnetsurfp.py $j.fa &> $j-rsa.out &"
     fi
     if [ !  -e $j.ss2 ]
     then
         #srun -A snic2015-10-12 --time=05:30:00 -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/addss.pl $i &> $j-adss.out &
         echo "#srun -A snic2015-10-12 --time=01:30:00 -n 1 -c 6 $HOME/git/PconsC3/extra/arne/MSA/addss.pl $i &> $j-adss.out &"
     fi
     if [ !  -e $j.gdca ]
     then
         #srun -A snic2015-10-12 --time=04:00:00 -n 1 -c 6 $HOME/git/PconsC3/rungdca.py $i &> $j-gdca.out &
         echo "#srun -A snic2015-10-12 --time=04:00:00 -n 1 -c 6 $HOME/git/PconsC3/rungdca.py $i &> $j-gdca.out &"
     fi
     if [ !  -e $j.0.02.plm20 ]
     then
         #srun -A snic2015-10-12 --time=24:00:00 -n 1 -c 6 $HOME/git/PconsC3/runplm.py $i &> $j-plm.out &
         echo "#srun -A snic2015-10-12 --time=24:00:00 -n 1 -c 6 $HOME/git/PconsC3/runplm.py $i &> $j-plm.out &"
     fi
    if [ !  -e $j.PconsC3.l5 ]
    then
	if [ -e $j.gdca ]  && [ -e $j.0.02.plm20 ]   && [ -e $j.rr ]  && [ -e $j.fa.rsa ]  && [ -e  $j.ss2 ] && [ -e  $j.gneff ]  && [ -e $j.trimmed ]
	then
            echo " #srun -A snic2015-10-12 --time=24:00:00 -n 1 -c 6 $HOME/git/PconsC3/predict.py $j.gdca $j.0.02.plm20 $j.rr $j.fa.rsa $j.ss2 $j.gneff $j.trimmed $HOME/git/PconsC3/ 5 $j.PconsC3  &> $j-predict.out &"
            #srun -A snic2015-10-12 --time=24:00:00 -n 1 -c 6 $HOME/git/PconsC3/predict.py $j.gdca $j.0.02.plm20 $j.rr $j.fa.rsa $j.ss2 $j.gneff $j.trimmed $HOME/git/PconsC3/ 5 $j.PconsC3  &> $j-predict.out &
	else 
	    ls -l $j.gdca $j.0.02.plm20 $j.rr $j.fa.rsa $j.ss2 $j.gneff $j.trimmed 
	fi
    fi
done
