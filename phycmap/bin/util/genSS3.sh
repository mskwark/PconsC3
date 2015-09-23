#!/bin/bash

if [ $# -ne 2 ]
then
        echo "Usage: ./PSIPRED <input_file> <out_dir> "
        exit
fi

# $1 the sequence name with suffix .seq
# $2 out directory

RaptorX_HOME=/hd0/phycmap.release/bin
###### BY TINA, May 5, 2003
BLASTPGP=$RaptorX_HOME/util/BLAST/bin/blastpgp
MAKEMAT=$RaptorX_HOME/util/BLAST/bin/makemat
BLASTDB=$RaptorX_HOME/databases/NR_new/nr90
PSIPREDDIR=$RaptorX_HOME/util/PSIPRED
###### BY TINA, May 5, 2003

DESTDIR=$2

###### BY TINA, May 6, 2003
# first make PSP directory if it doesn't exist
if [ ! -d $DESTDIR ] ; then
    mkdir $DESTDIR
fi
###### BY TINA, May 6, 2003

bname=`basename $1 .seq`
rootname=R$bname

echo $rootname 

cat  $1 > $rootname.seq

echo "Running PSI-BLAST with sequence" $1 "...."
#$BLASTPGP -b 0 -j 3 -h 0.0008 -d $BLASTDB -i $1 -C $rootname.chk
#####$BLASTPGP -F T -b 0 -a 22 -j 3 -h 0.001 -d $BLASTDB -i $1 -C $rootname.chk

cp $RaptorX_HOME/tmp/$bname"_sse".chk $rootname.chk
echo $rootname.chk > $rootname.pn
echo $1 > $rootname.sn


#echo $rootname.chk > $rootname.pn
#echo $1 > $rootname.sn
$MAKEMAT -P $rootname

echo Pass1 ....

$PSIPREDDIR/bin/psipred $rootname.mtx $PSIPREDDIR/data/weights.dat $PSIPREDDIR/data/weights.dat2 $PSIPREDDIR/data/weights.dat3 > $rootname.ss

echo Pass2

$PSIPREDDIR/bin/psipass2 $PSIPREDDIR/data/weights_p2.dat 1 0.98 1.09 $rootname.ss2 $rootname.ss > $rootname.horiz


echo "Final output files:" $rootname.ss2 $rootname.horiz $rootname.ss
mv $rootname.ss $DESTDIR/$bname.ss
mv $rootname.ss2 $DESTDIR/$bname.ss2
mv $rootname.horiz $DESTDIR/$bname.horiz

#remove temporary files
echo Cleaning up ....
#-rm -f $rootname.pn $rootname.sn $rootname.aux error.log $rootname.mtx
#-rm -f $rootname.chk 
rm -f $rootname.*
