#!/bin/sh
#predict contact map for one sequence using SVM.
if [ $# -ne 3 ]
then
	echo "need 3 parameters:fasta seq_file, output file alignmentfile."
	exit 1
fi
/scratch/arne/contactpreds/chaperonin/bin/predict_map.pl /scratch/arne/PconsC2-extra/svmcon1.0//script /scratch/arne/PconsC2-extra/SCRATCH-1D_1.1/bin/run_SCRATCH-1D_predictors.sh /scratch/arne/PconsC2-extra/svmcon1.0/server/svm_classify /scratch/arne/PconsC2-extra/svmcon1.0/model//model.g3 $1 $2 $3
