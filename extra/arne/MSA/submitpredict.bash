#!/bin/bash 
#SBATCH -A snic2015-10-12                                                                                     
#SBATCH -n 1                                                                                                  
#SBATCH -c 1                                                                                                  
#SBATCH --time=120:00:00
#SBATCH -J PconsC3-Predict
#SBATCH --output=out/PconsC3.%J.out
#SBATCH --error=err/PconsC3.%J.err

rsync -arv /proj/bioinfo/software/PconsC3.julia/trees/ /dev/shm/

~/git/PconsC3/predict.py alpha_beta.gdca alpha_beta.0.02.plm20 alpha_beta.rr alpha_beta.fasta.rsa alpha_beta.ss2 alpha_beta.gneff alpha_beta.trimmed /dev/shm/ PconsC3-shm.out
