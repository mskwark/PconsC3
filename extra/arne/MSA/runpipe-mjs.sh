#!/bin/bash

#this script can be used for webserver other than standalone, since we do not
#have a slim version of epad and matlab bioinformatics package, and R
error_file_io="Writing or reading files error!"
error_buildfeature="Build feature failed!"
error_epad="";
error_getfeature="";
error_rrr="rrr.pl error";
error_ilp="Ilp running error!";

seqfile=$1
seqbase=`basename $seqfile|replace '.seq' '' `;
pdbid=$seqbase;
echo "PDBID: $pdbid";
if [ ! -e $seqfile.a3m ]; then
    echo "The A3M formatted alignment $seqfile.a3m does not exist!";
    exit;
fi

while [ "$1" != "" ]; do
    case $1 in
        -cpu ) shift
	    BLAST_CPU=$1
	    ;;
    esac
    shift
done

if [ $BLAST_CPU -gt 1 &> /dev/null ] ; then :
    echo "Using $BLAST_CPU cpu in blast searching...";
else
    BLAST_CPU=1
fi
install_dir=/net/radon/hd0/phycmap.release
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$install_dir/bin
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
if [ $? -ne 0 ];then echo $error_file_io ; exit $? ; fi

cp $seqfile $workdir
cp $seqfile.a3m $workdir/$seqbase.a3m
$bindir/reformat.pl a3m clu $seqfile.a3m $workdir/$seqbase.clu;
seqfile=`basename $seqfile`
sendbackdir=$3
if [[ $# == 3 && $sendbackdir == "" ]]; then
sendbackdir=$currdir;
echo $sendbackdir
else
    sendbackdir=$currdir
fi

cd $workdir

echo "PhyCMAP: Generating features... $workdir $pdbid"

tgtfile="../$pdbid.tgt"
a3mfile="../$pdbid.a3m"

if [[ -f $tgtfile && -f $a3mfile ]] ; then
ln -s $tgtfile $workdir;
ln -s $a3mfile $workdir;
tgtfile=`basename $tgtfile`;
a3mfile=`basename $a3mfile`;

else 
( cd $CNFSEARCHDIR;
ln -s $workdir/$seqbase.a3m $CNFSEARCHDIR/tmp/$seqbase.a3m;
tail -n +2 $workdir/$seqbase.clu > $CNFSEARCHDIR/tmp/$seqbase.clu;
echo  "./buildFeature -i $workdir/$seqfile -o $pdbid.tgt -c $BLAST_CPU &> $workdir/buildFeature.log";
./buildFeature -i $workdir/$seqfile -o $pdbid.tgt -c $BLAST_CPU &> $workdir/buildFeature.log ; 
if [ $? -ne 0 ];then echo $error_buildfeature ; exit -1 ; fi
ln -s $CNFSEARCHDIR/$pdbid.tgt  $workdir ;
ln -s $CNFSEARCHDIR/tmp/$seqbase.a3m $workdir/$pdbid.a3m
) 
if [ $? -ne 0 ];then exit -1 ; fi


tgtfile=$pdbid.tgt
a3mfile=$pdbid.a3m

fi

echo "TGT: $tgtfile  A3M: $a3mfile";

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
) &> pretag.log
#exit 0;

(cd $epadCaDir; ./bin/EPADCalc -e EPAD -A $workdir/$pdbid.pretag.ca -F 11 | replace '*' '' > $workdir/$pdbid.epadca.prob ; exit $? ) &> epadca.log
if [ $? -ne 0 ] ;then echo $error_epad ; exit -1 ;fi

(cd $epadCbDir; ./bin/EPADCalc -e EPAD -A $workdir/$pdbid.pretag.cb -F 11 | replace '*' '' > $workdir/$pdbid.epadcb.prob ; exit $? ) &> epadcb.log
if [ $? -ne 0 ] ;then echo $error_epad ; exit -1 ;fi


bpsfile="$pdbid.bps.csv"
mifile="$pdbid.zydi-0"
moreevfile="$pdbid.moreev.csv"
echo $bindir/getProteinFeature -seq $seqfile -msa $a3mfile.fasta -mifile $moreevfile -zy0file $pdbid.zydi-0 -bpsfile $pdbid.bps.csv
$bindir/getProteinFeature -seq $seqfile -msa $a3mfile.fasta -mifile $moreevfile -zy0file $pdbid.zydi-0 -bpsfile $pdbid.bps.csv &> getproteinfeature.log

if [ $? -ne 0 ] ;then echo $error_getfeature ; exit -1 ;fi

#compute epadca epadcb, require epad
#for 6125 it has been computed

echo "PhyCMAP: Predicting pairwise probabilities..."

rfpredfile="$pdbid.predcb"
tempoutfile="$pdbid.pred.feature"

$bindir/rrr.pl -lib $PDBTOOLS_DIR   -pdb $pdbid  -act predict  -tpl $tgtfile -epadca $workdir/$pdbid.epadca.prob -epadcb $workdir/$pdbid.epadcb.prob -bps $bpsfile -mi $mifile -out $tempoutfile -outfile $rfpredfile -modelFile $bindir/model_rf379_24up_cb_new  -r_exe `which R`  -methodStr rf379 -featureSetStr 3:379  -Routputfile $pdbid.rout -evfile $moreevfile &> r.stdout 

if [ $? -ne 0 ] ;then echo $error_rrr ; exit -1 ;fi

echo "PhyCMAP: Computing the contact map with constraints..."

$bindir/runsdcp.pipe.sh $pdbid $rfpredfile $tgtfile $bindir &> runsdcp.log
if [ $? -ne 0 ] ;then echo $error_ilp ; exit -1 ;fi


cp $workdir/$pdbid.rr $currdir/$pdbid.rrunsort
(head -n5 $currdir/$pdbid.rrunsort ; tail -n +6 $currdir/$pdbid.rrunsort |sort -n -r -k5 ) > $pdbid.rr
 
$bindir/rr_format.sh $pdbid.rr > $pdbid.rr2 
mv $pdbid.rr2  $currdir/$pdbid.rr
rm $currdir/$pdbid.rrunsort


if [ "$currdir" != "$install_dir/test" ] ; then
rm -rf $workdir ;
fi

