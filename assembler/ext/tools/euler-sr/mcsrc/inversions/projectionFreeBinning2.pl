#!/usr/bin/env perl


use common;

if ($#ARGV != 3) {
  print "usage: $0 invFile speciesFile size_dir outfile\n";
  exit(0);
}

$invFile = shift @ARGV;
$speciesFile = shift @ARGV;
$sizeDir = shift @ARGV;
$outFile = shift @ARGV;

open (OUT, ">$outFile") or die "cannot open $outFile\n";

# Grab the sequences
@refSpecies;
%pairs;
common::ReadInvFile($invFile, \@refSpecies, \%pairs);
print "done\n";
# Read what species to care about
@species = ();
common::ReadSpeciesFile($speciesFile, \@species);

# Get the sequence lengths
%sequenceLengths = ();

common::GetSequenceLengths("$sizeDir/*.size", \%sequenceLengths);

@bins = ();

# clean up the pairs hash by removing empty elements
@refKeys = common::OrderedKeys(\%pairs);
foreach $ref (@refKeys) {
  @qryKeys = common::OrderedKeys(\%{$pairs{$ref}});
  foreach $qry (@qryKeys) {
    if ($#{$pairs{$ref}{$qry}} < 0) {
      print "removing empty: $ref $qry\n";
      delete $pairs{$ref}{$qry};
    }
  }
  @qryKeys = common::OrderedKeys(\%{$pairs{$ref}});
  if ($#qryKeys < 0) {
    print "removing empty $ref\n";
    delete $pairs{$ref};
  }
}


@bins = ();

@refKeys = common::OrderedKeys(\%pairs);
foreach $ref (@refKeys) {
  print "$ref\n";
  @qryKeys = common::OrderedKeys(\%{$pairs{$ref}});
  foreach $qry (@qryKeys) {
    for($inv = 0; $inv <= $#{$pairs{$ref}{$qry}}; $inv++ ) {
      $binFound = 0;
      $qryStart = $sequenceLengths{$qry} - $pairs{$ref}{$qry}[$inv][5];
      $qryEnd   = $sequenceLengths{$qry} - $pairs{$ref}{$qry}[$inv][4];
      $refStart = $pairs{$ref}{$qry}[$inv][1];
      $refEnd = $pairs{$ref}{$qry}[$inv][2];
#      print "$ref $qry $refStart $refEnd $qryStart $qryEnd\n";
      # look through bins for one that fits this inversion
      for ($bin = 0; $bin < scalar @bins and $binFound == 0; $bin++ ) {
#	print "  bin: $bin\n";
	for ($binit = 0; $binit < @{$bins[$bin]} and $binFound == 0; $binit++) {
	  $refOvpLen = 0;
	  $qryOvpLen = 0;
	  if ($bins[$bin][$binit][0] eq $ref) {
	    # check to see if this overlaps the ref seq
	    $refOvpLen = common::ComputeOverlap($refStart, $refEnd,
						$bins[$bin][$binit][1], $bins[$bin][$binit][2]);
	  }
	  if ($bins[$bin][$binit][0] eq $qry) {
	    # check to see if this overlaps the ref seq
	    $qryOvpLen = common::ComputeOverlap($qryStart, $qryEnd,
						$bins[$bin][$binit][1], $bins[$bin][$binit][2]);
	  }
#	  print "    item : $binit $bins[$bin][$binit][0] $bins[$bin][$binit][1], $bins[$bin][$binit][2] $refOvpLen $qryOvpLen\n";
	  if ($refOvpLen > 0 or $qryOvpLen > 0) {
	    # found a bin for this inversion. store it
#	    print "adding $ref $qry $inv to $bin\n";
	    push @{$bins[$bin]}, [$ref, $refStart, $refEnd];
	    push @{$bins[$bin]}, [$qry, $qryStart, $qryEnd];
	    $binFound = 1;
	  }
	}
      }
      if ($binFound == 0) {
	$last = $#bins;
	$last++;
#	print "creating bin: $last\n";
	@bins[$last] = ();
	push @{$bins[$last]}, [$ref, $pairs{$ref}{$qry}[$inv][1], $pairs{$ref}{$qry}[$inv][2]];
	push @{$bins[$last]}, [$qry, $qryStart, $qryEnd];
      }
    }
  }
}

# the bins themselves may overlap, condense the bins

#$b1 = 0;
#while ($b1 < $#bins) {
#  $b2 = $b1 + 1;
#  while ($b2 < $#bins) {
#    $foundOverlap = 0;
#    for ($bi1 = 0; $bi1 <= $#{$bins[$b1]} and $foundOverlap==0; $bi1++) {
#      for ($bi2 = 0; $bi2 <= $#{$bins[$b2]} and $foundOverlap == 0; $bi2++ ) {
#	if ($bins[$b1][$bi1][0] eq $bins[$b2][$bi2][0]) {
#	  $overlap = common::ComputeOverlap($bins[$b1][$bi1][1], $bins[$b1][$bi1][2],
#					    $bins[$b2][$bi2][1], $bins[$b2][$bi2][2]);
#	  if ($overlap > 0) {
#	    $foundMatch = 1;
#	    for ($mb = 0; $mb <= $#{$bins[$b2]}; $mb++ ) {
#	      push @{$bins[$b1]}, @{$bins[$b2]}[$mb];
#	    }
#	    splice @bins, $b2, 1;
#	  }
#	} # end checking overlap
#      }
#    } # end looking at all pairs of coordinates
#    if ($foundOverlap == 0) {
#      $b2++;
#    }
#  }
#  $b1++;
#}

for ($b = 0; $b <= $#bins; $b++ ) {
  print OUT "bin: $b\n";
  for ($bi = 0; $bi <= $#{$bins[$b]}; $bi++) {
    print OUT "  $bins[$b][$bi][0] $bins[$b][$bi][1] $bins[$b][$bi][2]\n";
  }
}

sub max {
  my ($a, $b) = @_;
  if ($a > $b) { return $a; }
  return $b;
}

sub min {
  my ($a, $b) = @_;
    if ($a < $b) { return $a;}
  return $b;
}
