#!/usr/bin/env perl

use common;
if ($#ARGV != 3) {
  print "usage: $0 speciesFile charFile localSpeciesFile outFile\n";
  exit(0);
}
$speciesFileName      = shift @ARGV;
$charFileName         = shift @ARGV;
$localSpeciesFileName = shift @ARGV;
$outFileName          = shift @ARGV;
%pairs = ();
@refSpecies = ();
@bins = ();
@boundaries = ();
common::ReadCharFile($charFileName, \@bins, \@boundaries, \@localSpecies);

@species = ();
common::ReadSpeciesFile($speciesFileName, \@species);


print "got species: @species\n";

@localSpecies = ();
common::ReadSpeciesFile($localSpeciesFileName, \@localSpecies);


# the species list appear in order of appearance of each 

for $binNum (0 .. $#bins) {
  for $binIdx (0 .. $#{$bins[$binNum]}) {
    @newCharList = ();
    $refSpecIndex = 0;
    $localSpecIndex = 0;
    $pos = 0;

    while ( $refSpecIndex <= $#species or
	    $localSpecIndex <= $#localSpecies) {
      if ($refSpecIndex <= $#species and
	  $localSpecIndex <= $#localSpecies and
	  $species[$refSpecIndex] eq $localSpecies[$localSpecIndex]) {
	# the two sequences have the same corresponding char
	$refSpecIndex++;
	$localSpecIndex++;
	push @newCharList, $bins[$binNum][$binIdx][4][$pos];
	$pos++;
      } elsif ($localSpecIndex > $#localSpecies or
	       ($refSpecIndex <= $#species and
		$species[$refSpecIndex] lt $localSpecies[$localSpecIndex])) {
	# something exists in ref that's not in qry
	# add an unknown

	push @newCharList, 2;
	$refSpecIndex++;
      } elsif ($refSpecIndex >  $#species or
	       ($localSpecIndex < $#localSpecies and
		$localSpecies[$localSpecIndex] lt $species[$refSpecIndex])) {
	$localSpecIndex++;
	$pos++;
      }
    }
    $bins[$binNum][$binIdx][4] = [ @newCharList ];
  }
}

common::PrintCharFile($outFileName, \@species, \@bins, \@boundaries);

print "got local species: @localSpecies\n";
print "got species @species\n";
