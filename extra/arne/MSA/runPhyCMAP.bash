#!/bin/bash -x


if [ -n ` hostname | grep  hpc2n ` ] 
then
    install_dir="/pfs/nobackup/home/a/arnee/Software/PconsC2-extra/PhyCmap/phycmap.release/"
elif  [ -n ` hostname | grep  nsc ` ] 
then
    install_dir="/proj/bioinfo/software/PconsC2-extra/PhyCmap/phycmap.release/"
else
    install_dir="/scratch/arne/PconsC2-extra/PhyCmap/phycmap.release/"
fi

#this script can be used for webserver other than standalone, since we do not
#have a slim version of epad and matlab bioinformatics package, and R
error_file_io="Writing or reading files error!"
error_buildfeature="Build feature failed!"
error_epad="epad error";
error_getfeature="getfeature error";
error_rrr="rrr.pl error";
error_ilp="Ilp running error!";

seqfile=$1
MSAFILE=$2
cur
seqbase=`basename $seqfile|sed -e s/\.seq$// |sed -e s/\.fa$// |sed -e s/\.fasta$//  `;
while [ "$1" != "" ]; do
    case $1 in
        -cpu ) shift
	    BLAST_CPU=$1
	    ;;
    esac
    shift
done
if [ $# -gt 2 ]; then 
    pdbid=$3
else
    pdbid=$seqbase
fi



if [ $BLAST_CPU -gt 1 &> /dev/null ] ; then :
    echo "Using $BLAST_CPU cpu in blast searching...";
else
    BLAST_CPU=1
fi

echo "$BLAST_CPU : $pdbid : $seqfile : $seqbase"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$install_dir/bin
export BIOPERL_DIR=$install_dir/bin/BioPerl-1.6.1/
export PDBTOOLS_DIR=$install_dir/bin/pdbtools/


CNFSEARCHDIR=$install_dir/bin/
epadCaDir="$install_dir/bin/epad.bin/"
epadCbDir="$install_dir/bin/epad.cb/"
bindir="$install_dir/bin"
pretagdir=$install_dir/bin/Pretag_To_EPAD/
currdir=`pwd -P`
workdir=`echo '$dir=int(rand(100000));$dir=".phycmaptmp.$dir.$ARGV[0]";if(-d $dir||-e $dir){ }else{print $dir}' | perl - $pdbid` 
workdir=`pwd`/$workdir

mkdir -p $workdir
if [ $? -ne 0 ];then echo "ERR1 $error_file_io" ; exit $? ; fi

cp $seqfile $workdir
seqfile=`basename $seqfile`
sendbackdir=$3
if [[ $# == 3 && $sendbackdir == "" ]]; then
sendbackdir=$currdir;
echo $sendbackdir
else
    sendbackdir=$currdir
fi

cd $workdir

echo "PhyCMAP: Generating features..."

tgtfile="../$pdbid.tgt"
a3mfile="../$pdbid.a3m"

if [[ -f $tgtfile && -f $a3mfile ]] ; then
ln -s $tgtfile $workdir
ln -s $a3mfile $workdir
tgtfile=`basename $tgtfile`;
a3mfile=`basename $a3mfile`;

else 
( cd $CNFSEARCHDIR;
./buildFeature -i $workdir/$seqfile -o $pdbid.tgt -c $BLAST_CPU &> $workdir/buildFeature.log ; 
if [ $? -ne 0 ];then echo "ERR2 $error_buildfeature" ; exit -1 ; fi
ln -s $CNFSEARCHDIR/$pdbid.tgt  $workdir ;
ln -s $CNFSEARCHDIR/tmp/$seqbase.a3m $workdir/$pdbid.a3m
) 
if [ $? -ne 0 ];then echo "ERR3" ; exit -1 ; fi


tgtfile=$pdbid.tgt
a3mfile=$pdbid.a3m

fi


$bindir/reformat.pl -r -noss $a3mfile $a3mfile.fasta &> $workdir/reformat.log
touch $pdbid.rr
#compute the tgt file and a2m file
#copy raptorx2:/home/majianzhu/LRR/CNFsearch and setup it!
a2mfile="$pdbid.a2m"
$bindir/converta3m2a2m.sh $a3mfile $tgtfile $bindir &> $workdir/converta3m.log
#compute the epad file, depend on matlab bioinformatics package,R, [comptue-0-4 ]

(cd $pretagdir;
./TGT_To_Pretag -i $workdir/$tgtfile -o $workdir/$pdbid.pretag.ca -s 1 -m 1 -H 1
./TGT_To_Pretag -i $workdir/$tgtfile -o $workdir/$pdbid.pretag.cb -s 2 -m 1 -H 2 -c 1
) &> $workdir/pretag.log
#exit 0;

(cd $epadCaDir; ./bin/EPADCalc -e EPAD -A $workdir/$pdbid.pretag.ca -F 11 | sed -E "s/\*//g" > $workdir/$pdbid.epadca.prob ; exit $? ) &> $workdir/epadca.log
if [ $? -ne 0 ] ;then echo "ERR4 $error_epad" ; exit -1 ;fi

(cd $epadCbDir; ./bin/EPADCalc -e EPAD -A $workdir/$pdbid.pretag.cb -F 11 | sed -E "s/\*//g" > $workdir/$pdbid.epadcb.prob ; exit $? ) &> $workdir/epadcb.log
if [ $? -ne 0 ] ;then echo "ERR5 $error_epad" ; exit -1 ;fi


bpsfile="$pdbid.bps.csv"
mifile="$pdbid.zydi-0"
moreevfile="$pdbid.moreev.csv"
$bindir/getProteinFeature -seq $seqfile -msa $a3mfile.fasta -mifile $moreevfile -zy0file $pdbid.zydi-0 -bpsfile $pdbid.bps.csv &> $workdir/getproteinfeature.log

if [ $? -ne 0 ] ;then echo "ERR6 $error_getfeature" ; exit -1 ;fi

#compute epadca epadcb, require epad
#for 6125 it has been computed

echo "PhyCMAP: Predicting pairwise probabilities..."

rfpredfile="$pdbid.predcb"
tempoutfile="$pdbid.pred.feature"

a=`which R`
if [ $? -ne 0 ] ;then echo "ERR-R R not installed" ; exit -1 ;fi

$bindir/rrr-org.pl -evfile $moreevfile -lib $PDBTOOLS_DIR   -pdb $pdbid  -act predict  -tpl $tgtfile -epadca $workdir/$pdbid.epadca.prob -epadcb $workdir/$pdbid.epadcb.prob -bps $bpsfile -mi $mifile -out $tempoutfile -outfile $rfpredfile -modelFile $bindir/model_rf379_24up_cb_new  -r_exe `which R`  -methodStr rf379 -featureSetStr 3:379  -Routputfile $pdbid.rout  &> $workdir/r.stdout 

if [ $? -ne 0 ] ;then echo "ERR7 $error_rrr" ; exit -1 ;fi

echo "PhyCMAP: Computing the contact map with constraints..."

$bindir/runsdcp.pipe.sh $pdbid $rfpredfile $tgtfile $bindir &> $workdir/runsdcp.log
if [ $? -ne 0 ] ;then echo "ERR8 $error_ilp" ; exit -1 ;fi


cp $workdir/$pdbid.rr $currdir/$pdbid.rrunsort
(head -n5 $currdir/$pdbid.rrunsort ; tail -n +6 $currdir/$pdbid.rrunsort |sort -n -r -k5 ) > $pdbid.rr
 
$bindir/rr_format.sh $pdbid.rr > $pdbid.rr2 
mv $pdbid.rr2  $currdir/$pdbid.rr
#rm $currdir/$pdbid.rrunsort


if [ "$currdir" != "$install_dir/test" ] ; then
#rm -rf $workdir ;
fi

