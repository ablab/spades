#!/usr/bin/env perl


use common;

if ($#ARGV != 2) {
  print "usage: $0 binFile seqLendir outfile\n";
  exit(0);
}
print "Deprecated.  Use BinInversions.pl\n";
exit(0);


$binFile = shift @ARGV;
$seqLenDir = shift @ARGV;
$outFile = shift @ARGV;


%seqLengths = ();

common::GetSequenceLengths("$seqLenDir/*.size", \%seqLengths);


open (BINFILE, "$binFile") or die "cannot open $binFile\n";
open (OUTFILE, ">$outFile") or die "cannot open $outFile\n";
@bins = ();
%coords = ();
$binNumber = -1;
while (<BINFILE>) {
  $line = $_;
  if ($line =~ /bin/) {
    # process an old bin
    @binSpecList = keys %coords;
    if ($#binSpecList >= 0 ) {
      print OUTFILE "bin: $binNumber\n";
      foreach $binSpec (@binSpecList) {
	# find the coordinates of the inversion for species
	# $binSpec
	$start = -1;
	$end   = -1;
	for ($c = 0; $c <= $#{$coords{$binSpec}}; $c++ ) {
	  if ($start == -1 or $start > $coords{$binSpec}[$c][0]) {
	    $start = $coords{$binSpec}[$c][0];
	  }
	  if ($end == -1 or $end < $coords{$binSPec}[$c][1]) {
	    $end = $coords{$binSpec}[$c][1];
	  }
	}
	print OUTFILE  "  $binSpec $start $end\n";
      }
    }
    # start a new bin
    ++$binNumber;
    %coords = ();
  }
  else {
    $line =~/\s*(\S+)\s+(\-?\d+)\s+(\-?\d+)/;
    $spec = $1;
    $start = $2;
    $end = $3;
    if ($start < 0 or $end < 0) {
      # quick sanity check
      if ($end >=0 or $start >= 0 ){
	print "error, both coords must be either above or below 0. $line\n";
	exit(0);
      }
      $revStart = $start; $revEnd = $end;
      if (! exists $seqLengths{$spec}) {
	print "needed length for $spec \n";
	exit(0);
      }
      $start = $seqLengths{$spec} + $revEnd + 1;
      $end   = $seqLengths{$spec} + $revStart + 1;
    }
    push @{$coords{$spec}}, [ $start, $end];
  }
}
