#!/usr/bin/env perl

if ($#ARGV != 3) {
  print "Find expected dist of random sequences generated with\n";
  print "a markov model.\n";
  print "usage: $0 model_file seq_len num_runs out_file [probelen]\n";
  exit(0);
}
$modelFile = shift @ARGV;
$seqLen    = shift @ARGV;
$numRuns   = shift @ARGV;
$outfile   = shift @ARGV;

$probeLen = 15;
if ($#ARGV >= 0) {
  $probeLen = shift @ARGV;
}
$wordLen = 11;
if ($#ARGV >= 0) {
  $wordLen = shift  @ARGV;
}

$seqAFile = "A.$$";
$seqBFile = "B.$$";
$alnFile = "AB.aln.$$";

$genA = "/home/mchaisso/UTILS/model_dna/mmgenseq.pl $modelFile $seqLen $seqAFile";
$genB = "/home/mchaisso/UTILS/model_dna/mmgenseq.pl $modelFile $seqLen $seqBFile";


$countCommon = "/home/mchaisso/UNICHRPRO/par_unicotup/LINUX/tupal -i $seqAFile -i $seqBFile -p $probeLen -w $wordLen -l $probeLen -o $alnFile";

open(OUT, ">$outfile") or die "cannot open $outfile\n";
for $i (0 .. $numRuns-1) {
  `$genA`;
  `$genB`;
  if (! -e $seqAFile or ! -e $seqBFile) {
    print "error generating sequence file \n";
    cleanup();
    exit (0);
  }
  `$countCommon`;
  $res = $?;
  if ($res == 0 and -e $alnFile) {
    $res = `wc -l $alnFile`;
    chomp $res;
    $res =~/\s+(\d+)\s+.*/;
    $numMatch = $1;
    print "got $numMatch matches $res\n";
    print OUT "$numMatch\n";
  }
  cleanup();
}
close OUT;


sub cleanup {
  if (-e $seqAFile ) { `rm $seqAFile`;}
  if (-e $seqBFile ) { `rm $seqBFile`;}
  if (-e $alnFile )  { `rm $alnFile`;}
}  
