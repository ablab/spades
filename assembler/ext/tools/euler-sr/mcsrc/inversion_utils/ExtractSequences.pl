#!/usr/bin/env perl

if ($#ARGV != 3) {
  print "usage: $0 humanbed chimpbed humandir chimpdir\n";
  exit(0);
}

$humanBed = shift @ARGV;
$chimpBed = shift @ARGV;
$humanDir  = shift @ARGV;
$chimpDir  = shift @ARGV;

open(HF, $humanBed) or die "cannot open $humanBed\n";
open(CF, $chimpBed) or die "cannot open $chimpBed\n";

while (<HF>) {
  $hline = $_;
  $cline = <CF>;

  @hcoords = split(/\s+/,$hline);
  @ccoords = split(/\s+/,$cline);

  $humanOut = "human.fasta";
  $chimpOut = "chimp.fasta";

  $humanSeq = "$humanDir/$hcoords[0].fa";
  $chimpSeq = "$chimpDir/$ccoords[0].fa";
  print "extracting @hcoords @ccoords\n";
  system("~/projects/mcsrc/sequtils/powerpc/extractseq $humanSeq $hcoords[1] $hcoords[2] $humanOut");
  system("~/projects/mcsrc/sequtils/powerpc/extractseq $chimpSeq $ccoords[1] $ccoords[2] $chimpOut");

  @all = @hcoords;
  push @all, @ccords;
  $combinedSeqName = join("", @all);

  system("FastaToMatlabSeq.pl $humanOut > $combinedSeqName.txt");
  system("FastaToMatlabSeq.pl $chimpOut >> $combinedSeqName.txt");
}
