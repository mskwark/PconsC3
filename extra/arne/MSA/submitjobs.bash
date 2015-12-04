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
    j=`basename $i`
    if [ !  -e $j.out ]
    then
        srun -A snic2015-10-12 --time=12:00:00 -n 1 -c 8 ~/git/PconsC3/extra/arne/MSA/run_PconsC3.bash $i &> $j.out &
    fi
done

