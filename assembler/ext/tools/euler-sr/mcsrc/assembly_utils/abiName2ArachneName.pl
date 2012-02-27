#!/usr/bin/env perl

use ReadLibrary;


if ($#ARGV != 1) {
  print "\n usar: $0 reads_file.fasta output.fasta\n";
  exit(0);
}


$in = shift @ARGV;
$fastaOut = shift @ARGV;

open (FASTAIN, $in) or die "no puedo abrir $in\n";
open (FASTAOUT, ">$fastaOut") or die "no puedo escribir $xmlOut\n";
# Print xml header information


$numParsed = 0;
$numFailed = 0;
while (<FASTAIN>) {
  $line = $_;
  chomp $line;
  if ($line =~ />/) {
    ($baseName, $type, $dir, $name) = ReadLibrary::ParseABITitle($line);
    @data = [];
    print FASTAOUT ">$name\n";
  }
  else {
    print FASTAOUT "$line\n";
  }
}

close FASTAOUT;
