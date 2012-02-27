#!/usr/bin/env perl

use common;

if ($#ARGV < 2) {
  print "usage: $0 infile outfile spec1 [spec2 ...]\n";
  exit(0);
}
$infile = shift @ARGV;
$outfile = shift @ARGV;
$toRemove = ();
while ($#ARGV >= 0) {
  $val = shift @ARGV;
  push @delspec, $val;
}
print "removing @delspec\n";

@refSpecies = ();
@bins = ();
@boundaries = ();
@species = ();
@deletedSpecies = ();
common::ReadCharFile($infile, \@bins, \@boundaries, \@species);


@toRemove = ();
while ($#delspec >= 0 ) {
  $foundone = 0;
  for ($s = 0; $s <= $#species; $s++ ) {
    if ($species[$s] eq $delspec[0]) {
      push @toRemove, $s;
      shift @delspec;
      $foundone = 1;
    }
  }
  if ($foundone == 0) {
    print "Did not find @delspec[0], bailing out \n";
    exit(0);
  }
}
print "deleting: @toRemove\n";

# find out what species are gone
$i = 0;
$numDeleted = 0;
for $i (0 .. $#species) {
  # check to see if species i is removed
  $j = 0;
  while ($j <= $#toRemove) {
    if ($toRemove[$j] == $i) {
      push @deletedSpecies, $species[$i-$numDeleted];
      splice @species, $i - $numDeleted, 1;
      $numDeleted++;
      $j = $#toRemove+2;
    }
    else {
      $j++;
    }
  }
}
print "species: @species\n";

$b = 0;
print "nb: $#bins\n";
while ($b <= $#bins) {
  $s = 0;
  # remove all rows corresponding to a removed character
  while ($s <= $#{$bins[$b]}) {
    $j = 0;
    while ($j <= $#toRemove) {
      if ($deletedSpecies[$j] eq $bins[$b][$s][0]) {
	splice @{$bins[$b]}, $s, 1;
	# leave s constant (by decrementing before incrementing)
	$s--;
	# break out of the loop early
	$j = $#toRemove + 1;
      }
      $j++;
    }
    $s++;
  }
  # if all rows are removed, remove this character
  if ($#{$bins[$b]} < 0) {
    splice @bins, $b, 1;
    splice @boundaries, $b, 1;
  }
  else {
    # For each remaining row remove the characters corresponding 
    # to removed columns
    $s = 0;
    for $s (0 .. $#{$bins[$b]}) {
      # delete columns
      $j = $#toRemove;
#      print "matrix before: @{$bins[$b][$s][4]}\n";
      while ($j >= 0) {
#	print "splicing $toRemove[$j] @{$bins[$b][$s][4]}[$toRemove[$j]] ";
	splice @{$bins[$b][$s][4]}, $toRemove[$j], 1;
	$j--;
      }
#      print "after: @{$bins[$b][$s][4]}\n";
    }
    $b++;
  }
}

common::PrintCharFile($outfile, \@species, \@bins, \@boundaries);

