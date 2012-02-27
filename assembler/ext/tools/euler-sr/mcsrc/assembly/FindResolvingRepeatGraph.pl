#!/usr/bin/env perl

if ($#ARGV != 2) {
  print "usage: $0 fasta_file tuple_size tuple_incrment\n";
  print "builds de bruijn graphs on fasta_file until there is a tuple\n";
  print "size that creates a simple graph (two edges)\n";
  exit(0);
}
$fasta = shift @ARGV;
$tupleSize = shift @ARGV;
$step      = shift @ARGV;

$resolved = 0;

while ($resolved == 0) {
  $command = "~/projects/mcsrc/assembly/i686/debruijn $fasta $fasta.$tupleSize.gvz -tupleSize $tupleSize -edgeFile $fasta.$tupleSize.$edge";
  
  `$command`;
  $nEdge = `grep -c ">" $fasta.$tupleSize.$edge`;
  chomp $nEdge;
  print "$tupleSize $nEdge\n";
  if ($nEdge == 2) {
    $resolved = 1;
  }
  else {
    $tupleSize += $step;
  }
}
    
