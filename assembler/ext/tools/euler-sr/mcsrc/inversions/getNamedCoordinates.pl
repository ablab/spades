#!/usr/bin/env perl

use common;

$charFileName = shift @ARGV;
$target = shift @ARGV;
$region = shift @ARGV;
@refSpecies = ();
@bins = ();
@boundaries = ();
common::ReadCharFile($charFileName, \@bins, \@boundaries, \@localSpecies);


foreach $i (0 .. $#bins) {
  if ($boundaries[$i][0] != 0 and 
      $boundaries[$i][1] != 0) {
    # found a bin that has been mapped to the human sequence
    print "$target\thuman\tnone\t$boundaries[$i][0]\t$region\t1\n";
    print "$target\thuman\tnone\t$boundaries[$i][1]\t$region\t1\n";
  }
  else {
    # didn't find a bin that has been mapped to human.
    # just use the first entry that has coordinates
    $j = 0;
    while ($j < $#{$bins[$i]}) {
      if ($bins[$i][$j][1] != 0 and
	  $bins[$i][$j][2] != 0) {
	print "$target\thuman\tnone\t$bins[$i][$j][1]\t$region\t1\n";
	print "$target\thuman\tnone\t$bins[$i][$j][2]\t$region\t1\n";
	$j = $#{$bins[$i]};
      }
      $j++;
    }
  }
}

