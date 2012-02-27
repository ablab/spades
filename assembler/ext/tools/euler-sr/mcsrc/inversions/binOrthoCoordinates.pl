#!/usr/bin/env perl

use common;

if ($#ARGV != 3) {
  print "usage: $0 orthoFile dataFile outFile distance\n";
  exit(0);
}
$orthoFile = shift @ARGV;
$dataFile  = shift @ARGV;
$outFile   = shift @ARGV;
$dist      = shift @ARGV;

%orthoPos = ();
common::ReadOrthoCoordinates($orthoFile, \%orthoPos);
%data = ();
common::ReadDataFile($dataFile, \%data);


@bins = ();

$counter = 0;
%binnedSequences = ();

foreach $seqName (keys %orthoPos) {
  if (! exists $binnedSequences{$seqName}) {
    $orthStartPos = @{$orthoPos{$seqName}}[0];
    $orthEndPos   = @{$orthoPos{$seqName}}[1];
    $cognate = $data{$seqName}[-1];
    $cognateOrthStart = @{$orthoPos{$cognate}}[0];
    $cognateOrthEnd   = @{$orthoPos{$cognate}}[1];

    $binFound = 0;
    #  print "$seqName, $#bins\n";
    $counter++;
    if ($counter % 100 == 0) {
      print "$counter $#bins\n";
    }
    for ($b = 0; $b <= $#bins and $binFound == 0; $b++ ) {
      $closeEnough = 0;
      $max = $#{$bins[$b]};
      if ($#{$bins[$b]} > 25) {
	$max = 25;
      }
      for ($bi = 0; $bi <= $max and $closeEnough == 0; $bi++) {
	$startDist = abs($orthStartPos - $bins[$b][$bi][1]);
	$endDist   = abs($orthEndPos   - $bins[$b][$bi][2]);
	#      print "sd: $orthStartPos $bins[$b][$bi][1]  $bins[$b][$bi][2] $startDist $endDist\n";
	if ($startDist < $dist and $endDist < $dist) {
	  $closeEnough = 1;
	}
      }
      if ($closeEnough) {
	push  @{$bins[$b]}, [$seqName, $orthStartPos, $orthEndPos];
	# add the cognate sequence as well
	push  @{$bins[$b]}, [$cognate, $cognateOrthStart, $cognateOrthEnd];
	$binnedSequences{$seqName} = 1;
	$binnedSequences{$cognate} = 1;
	$binFound = 1;
      }
    }
    if ($binFound == 0) {
      $last = $#bins;
      push @{$bins[$last +1]}, [$seqName, $orthStartPos, $orthEndPos];
      push  @{$bins[$last +1]}, [$cognate, $cognateOrthStart, $cognateOrthEnd];
    }
  }
}
       
open (OUT, ">$outFile") or die "cannot open $outfile\n";
for ($b = 0; $b <= $#bins; $b++ ) {
  print OUT "bin: $b\n";
  for ($bi = 0; $bi <= $#{$bins[$b]}; $bi++) {
    print OUT "@{$bins[$b][$bi]}\n";
  }
  print OUT "\n\n";
}
