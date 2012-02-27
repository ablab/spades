#!/usr/bin/env perl

use ReadLibrary;

if ($#ARGV < 2) {
  print "usage: $0 mapfile pairfile namesfile [-t match threshold] [-n readsNamefile ]\n";
  print "Find matches in the mapfile that are of length at least minmatch\n";
  print "and print the corrresponding coordinates.\n";
  print "namesfile is a file of all read names \"brl001.titles\"\n";
  print "minmatch (100) is the minimum length match to bother printing\n";
  exit(0);
}


$in = shift @ARGV;
$pairfile = shift @ARGV;
$namefile = shift @ARGV;
$threshold = 100;
$readsNameFile = "";
while ($#ARGV >= 1) {
  $option = shift @ARGV;
  if ($option eq "-t") {
    $threshold = shift @ARGV;
  }
  elsif ($option eq "-n")  {
    $readNameFile = shift @ARGV;
  }
}

open(IN, "$in") or die "cannot open $in\n";
open(PAIRS, "$pairfile") or die "cannot open $pairfile\n";


# read in the names of EVERY read
@names = ();
open(NAM, "$namefile") or die "cannot open $namefile\n";
while (<NAM>) {
  ($a,$b,$c,$name) = ReadLibrary::ParseABITitle($_);
  push @names, $name;
}

if ($readNameFile ne "") {
  open(RN, ">$readNameFile") or die "cannot write to $readNameFile\n";
}

@pairs = ();
ReadLibrary::ReadMatePairIndexFile(\*PAIRS, \@pairs);
@locations = ();
@types     = ();
@directions = ();
$first = 1;
while(<IN>) {
  $line = $_;
  chomp $line;
  if ($line =~ /^>(.*)/) {
    # reset the hits list
    $title = $1;
    ($base, $dir, $type) = ReadLibrary::ParseABITitle($title);
    if ($first == 0) {
      push @locations, [@hits];
      push @directions, [@dirs];
      push @types, $type;
      $end = $#locations;
    }
    else {
      $first = 0;
    }
    @hits = ();
    @dirs = ();
  }
  else {
    @vals = split (/\s+/, $line);
    if ($vals[0] > $threshold) {
      push @hits, $vals[3];
      push @dirs, $vals[6];
    }
  }
}

for $l (0 .. $#locations) {
  if ($pairs[$l] != -1) {
    $mp = $pairs[$l];
    for $p (0 .. $#{$locations[$l]}) {
      for $q ( 0 .. $#{$locations[$mp]}) {
	if ($locations[$l][$p] < $locations[$mp][$q]) {
	  $first = $locations[$l][$p];
	  $second = $locations[$mp][$q];
	  $firstDir = $directions[$l][$p];
	  $secondDir = $directions[$mp][$q];
	  $firstName = $names[$l];
	  $secondName = $names[$mp];
	}
	else {
	  $first = $locations[$mp][$q];
	  $second =$locations[$l][$p];
	  $firstDir = $directions[$mp][$q];
	  $secondDir = $directions[$l][$p];
	  $firstName = $names[$mp];
	  $secondName = $names[$l];
	}
	print "$first $firstDir $second $secondDir $types[$l] $l $mp \n";
	if ($readNameFile ne "") {
	  print RN "$firstName $secondName\n";
	}
      }
    }
    # don't duplicate the data
    $pairs[$mp] = -1;
  }
}

if ($readNameFile ne "") {
  close RN;
}

