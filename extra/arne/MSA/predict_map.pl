#! /usr/bin/perl -w
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
my $info = Sys::Info->new;
my $cpu  = $info->device( CPU => %options );


##################################################################
# Predict contact map for one single sequence in fasta format
# Input: SVMcon path, predict_ssa.sh in SSpro 4 package, svm classifier of SVM-light, 
# trained SVM model, input sequence file in FASTA format, output file in CASP format
#
# Author: Jianlin Cheng
# January 2, 2006
##################################################################

if (@ARGV != 7)
{
	die "need 7 parameters: script dir (SVMcon dir), ssa predictor (SCRATCH), svm classifier (SVM-light), svm model (trained SVM model), input file in FASTA format, output file align_file.\n"; 
}
$script_dir = shift @ARGV;
$ssa_predictor = shift @ARGV;
$svm_predictor = shift @ARGV; 
$svm_model = shift @ARGV;
$fasta_file = shift @ARGV;
$output_file = shift @ARGV; 
$align_file = shift @ARGV; 

-d $script_dir || die "can't find script dir: $script_dir\n";

if (! -f $ssa_predictor)
{
	die "can't find the ss, sa predictor.\n"; 
}
if (! -f $svm_predictor)
{
	die "can't find svm classifier.\n"; 
}
if (! -f $svm_model)
{
	die "can't find the model definition file.\n"; 
}
if (! -f $fasta_file)
{
	die "can't find the fasta file.\n"; 
}
if (! -f $align_file)
{
	die "can't find the align file.\n"; 
}
open(FASTA, "$fasta_file") || die "can't open fasta file.\n"; 
$target_name = <FASTA>; 
close FASTA;
$target_name = substr($target_name,1); 

#generate alignment, predict ss and sa
$pos = rindex($fasta_file, "/"); 
if ($pos < 0)
{
	$ssa_file = $fasta_file; 
}
else
{
	$ssa_file = substr($fasta_file, $pos+1); 
}
#$ssa_file =~ s/\./_/g; 
$ssa_file .= ".ssa"; 
#$align_file = "${ssa_file}.align"; 
print "predict secondary structure and solvent accessibility...\n";
#print "$ssa_predictor $fasta_file $ssa_file\n"; 
system("$ssa_predictor $fasta_file $ssa_file $cpu->count"); 
#notice: two files are generated from ssa predictor: one is ssa output, one is alignment file. 

# Modifed to be able to use SCRATCH preditors (different format and separate files)
#create tmp file for contact map prediction
open(TMP, ">$output_file.tmp") || die "can't create temporary file.\n"; 
open(INPUT, "$fasta_file") || die "can't open the ssa file.\n"; 
$seq='';
while (<INPUT>) {
    chomp;
    if (!/\>/) {$seq .= $_;}
}
close INPUT;

open(INPUT, "$ssa_file.ss") || die "can't open the ssa file.\n"; 
$ss='';
while (<INPUT>){
    chomp;
    if (!/\>/) {$ss .= $_;}
}
close INPUT;

open(INPUT, "$ssa_file.acc") || die "can't open the ssa file.\n"; 
$sa='';
while (<INPUT>){
    chomp;
    if (!/\>/) {$sa .= $_;}
}
close INPUT;

#print "SEQ: ",$seq,"\n";
#print "SA: ",$sa,"\n";
#print "SS: ",$ss,"\n";

$length = length($seq); 
if ($length != length($ss) || $length != length($sa))
{
	die "sequence length doesn't match.\n"; 
}

for ($i = 0; $i < $length; $i++)
{
	if ($i < $length - 1)
	{
		print TMP substr($seq, $i, 1), " ";
	}
	else
	{
		print TMP substr($seq, $i, 1), "\n";
	}
}

for ($i = 0; $i < $length; $i++)
{
	if ($i < $length - 1)
	{
		print TMP substr($ss, $i, 1), " ";
	}
	else
	{
		print TMP substr($ss, $i, 1), "\n";
	}
}

for ($i = 0; $i < $length; $i++)
{
	if ($i < $length - 1)
	{
		if (substr($sa, $i, 1) eq "e" )
		{
			print TMP 50, " ";
		}
		else
		{
			print TMP 10, " "; 
		}
	}
	else
	{
		if (substr($sa, $i, 1) eq "e" )
		{
			print TMP 50, "\n";
		}
		else
		{
			print TMP 10, "\n"; 
		}
	}
}

for ($i = 0; $i < $length; $i++)
{
	if ($i < $length - 1)
	{
		print TMP "0 0 0 ";
	}
	else
	{
		print TMP "0 0 0\n";
	}
}
close TMP; 



# We also need to genera

#generate svm dataset for separation >=6
print "generate SVM dataset...\n";
#print "$script_dir/generate_input_with_title.pl $script_dir $output_file.tmp $align_file  6 8 > $output_file.svm\n";
system("$script_dir/generate_input_with_title.pl $script_dir $output_file.tmp $align_file  6 8 > $output_file.svm"); 

print "classify data points using SVM...\n";
#make svm predictions
system("$svm_predictor $output_file.svm $svm_model $output_file.res");

open(RES, "$output_file.res") || die "can't open result file.\n"; 
open(OUT, ">$output_file") || die "can't create output file.\n"; 
@scores = <RES>;
close RES;
close OUT; 

open(SET, "$output_file.svm") || die "can't read generate data points.\n";
@set = <SET>;
close SET;

#generate CASP format file
$casp = $output_file; 
open(CASP, ">$casp") || die "can't create casp file.\n";
print CASP "PFRMAT RR\n";
print CASP "TARGET $target_name";
print CASP "AUTHOR SVMcon\n";
print CASP "METHOD SVM contact map predictor (separation >= 6)\n";
print CASP "MODEL  1\n"; 

#pint out sequence
for ($i = 1; $i <= $length; $i++)
{
	print CASP substr($seq, $i-1, 1);
	if ($i % 50 == 0 || $i == $length)
	{
		print CASP "\n";
	}
}

#check consistency
@set == 2 * @scores || die "number of data points doesn't match with number of scores.\n";

while (@scores)
{
	$score = shift @scores;
	$pair = shift @set;
	shift @set;
	chomp $pair;
	$pair = substr($pair, 1);
	($id1, $id2) = split(/\s+/, $pair);

	if ($score > 0)
	{
		$sigmoid = 1 / (1 + exp(-$score) );
		$sigmoid *= 100000;
		$sigmoid = int($sigmoid);
		$sigmoid /= 100000;
		push @pairs, {
			id1 => $id1,
			id2 => $id2,
			score => $sigmoid
		};

	}

	#if ($score > 0)
	#{
	#	printf CASP "%3s%4s 0  8 ", $id1, $id2;
	#	#convert score to probability
	#	$sigmoid = 1 / (1 + exp(-$score) );
	#	printf CASP "%5s\n", $sigmoid;
	#}
}

@sorted_pairs = sort { $a->{"id2"} <=> $b->{"id2"} } @pairs;
@out_pairs = sort { $a->{"id1"} <=> $b->{"id1"} } @sorted_pairs;

for ($i = 0; $i < @out_pairs; $i++)
{
	printf CASP "%3s%4s 0  8 %5s\n", $out_pairs[$i]{"id1"},  $out_pairs[$i]{"id2"}, $out_pairs[$i]{"score"};
}

print CASP "END\n";


#remove temporary files
#`rm $output_file.tmp $output_file.res $output_file.svm`; 
#`rm $ssa_file`; 
#`rm $align_file`; 
