#!/usr/bin/env perl
if ($#ARGV < 0) {
  print "\n usage: $0 bedfile [ratio]\n";
  print " prints information to a file to be read by matlab PrintDotPlots.m\n\n"; 
  exit(0);
}
use POSIX;

open(IN, @ARGV[0]) or die "cannot open @ARGV[0]\n";
%hCoords = ();
%cCoords = ();
shift @ARGV;
$ratio = 3;
if ($#ARGV >= 0) {
  $ratio = shift @ARGV;
} 
while(<IN>) {
  $line = $_;
  ($human, $chimp) = split(/\s+/, $line);
  $human =~ /(chr\S+):(\d+)\-(\d+)/;
  $hChr = $1; $hS = $2; $hE = $3;
  $hLen = $hE- $hS + 1;
  $hSL = $hS - POSIX::floor($hLen*$ratio);
  $hEL = $hE + POSIX::floor($hLen*$ratio);

  $chimp =~ /(chr\S+):(\d+)\-(\d+)/;
  $cChr = $1; $cS = $2; $cE = $3;
  $cLen = $cE- $cS + 1;
  $cSL = $cS - POSIX::floor($cLen*$ratio);
  $cEL = $cE + POSIX::floor($cLen*$ratio);
  print "human.$hChr.$hSL.$hEL chimp.$cChr.$cSL.$cEL human.$hChr.$hS.$hE.chimp.$cChr.$cS.$cE.jpg #human $hChr $hS $hE chimpanzee $cChr $cS $cE#\n";
}


