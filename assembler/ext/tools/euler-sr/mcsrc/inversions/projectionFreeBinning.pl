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

# Read what species to care about
@species = ();
common::ReadSpeciesFile($speciesFile, \@species);
%validSpecies = ();
for ($s = 0; $s <= $#species; $s++ ) {
  $validSpecies{$species[$s]} = 1;
}
# Get the sequence lengths
%sequenceLengths = ();

common::GetSequenceLengths("$sizeDir/*.size", \%sequenceLengths);

@bins = ();

# clean up the pairs hash by removing empty elements
@refKeys = common::OrderedKeys(\%pairs);
foreach $ref (@refKeys) {
  @qryKeys = common::OrderedKeys(\%{$pairs{$ref}});
  if (exists $validSpecies{$ref} ) {
    foreach $qry (@qryKeys) {
      if ($#{$pairs{$ref}{$qry}} < 0 or (! exists $validSpecies{$qry})) {
        print "removing empty: $ref $qry\n";
        delete $pairs{$ref}{$qry};
      }
    }
  }
  if (exists $validSpecies{$ref}) {
    @qryKeys = common::OrderedKeys(\%{$pairs{$ref}});
  }
  if ($#qryKeys < 0 or ! exists $validSpecies{$ref}) {
    print "removing empty $ref\n";
    delete $pairs{$ref};
  }
}

$b = 0;
$pairsLeft = 1;
while ($pairsLeft == 1) {
  # Get the first element in the pairs hash
  @refKeys = (); @qryKeys = ();
  @refKeys = common::OrderedKeys(\%pairs);
  if ($#refKeys < 0) {
    last;
  }
  $pos = 0;
  if ($#refKeys < 0) {
    print "invalid pair structure\n";
    exit(0);
  }
#  print "ref keys: @refKeys\n";
  $ref = $refKeys[0];
  @qryKeys = common::OrderedKeys(\%{$pairs{$ref}});
  if ($#qryKeys < 0) {
    print "invalid pair structure qry for $ref $#qryKeys\n";
    exit(0);
  }
  $qry = @qryKeys[0];
  $refSpec = $pairs{$ref}{$qry}[0][0];
  $refStart = $pairs{$ref}{$qry}[0][1];
  $refEnd   = $pairs{$ref}{$qry}[0][2];
  $qrySpec = $pairs{$ref}{$qry}[0][3];
  if (! exists $sequenceLengths{$qrySpec} ) {
    print "could not length for $qrySpec\n";
    exit(0);
  }
  $qryStart = $sequenceLengths{$qrySpec} - $pairs{$ref}{$qry}[0][5];
  $qryEnd   = $sequenceLengths{$qrySpec} - $pairs{$ref}{$qry}[0][4];
  @toCheck = ();
  push @toCheck, [$refSpec, $refStart, $refEnd];
  push @toCheck, [$qrySpec, $qryStart, $qryEnd];
  $sRefSpec = $refSpec;
  $sQrySpec = $qrySpec;
#  print "finding overlaps for $ref $qry $start $end\n";
  %bin = ();
  print "new bin: $b $ref $qry $refStart $refEnd\n";
  @ref = (); @start = (); @end = ();
  $firstIter = 1;
  while ($#toCheck >= 0) {
    @overlapping = ();
    print "seed start:     $toCheck[0][0], $toCheck[0][1], $toCheck[0][2]\n";
    FindOverlappingIntervals(\%pairs,
			     $toCheck[0][0], $toCheck[0][1], $toCheck[0][2],
			     \%sequenceLengths,
			     0.0001,
			     \@overlapping, \@toCheck);
    shift @toCheck;
    if ($#overlapping >= 0) {
      for ($ovp = 0; $ovp <= $#overlapping; $ovp++ ) {
	$ref = $overlapping[$ovp][0];
	print "\tadding overlap $ref $start $end got overlap: @{$overlapping[$ovp]} \n";
	push @{$bin{$ref}}, [ @{$overlapping[$ovp]} ];
      }
    }
    elsif ($firstIter == 1) {
      print "Singleton bin, ref now: $sRefSpec, qry now: $sQrySpec\n";
      push @{$bin{$sRefSpec}}, [$sRefSpec, $refStart, $refEnd];
      push @{$bin{$sQrySpec}}, [$sQrySpec, -$qryStart, -$qryEnd];
      print "splicing @{$$pairs{$sRefSpec}{$sQrySpec}[0]}\n";
      splice @{$pairs{$sRefSpec}{$sQrySpec}}, 0, 1;
      if ($#{$pairs{$sRefSpec}{$sQrySpec}} < 0 ) {
	print "DELETING pair $sRefSpec $sQrySpec\n";
	delete $pairs{$sRefSpec}{$sQrySpec};
	$deletedQry = 1;
      }
      my @sQrySpecKeys = common::OrderedKeys(\%{$pairs{$sRefSpec}});
      if ($#sQrySpecKeys < 0) {
	print "DELETING REFERENCE!!! $sRefSpec\n";
	delete $pairs{$sRefSpec};
      }
    }
    $firstIter = 0;
    print "--------------------------------------\n\n\n";
  }
  push @bins, {%bin};
  $b++;

  print "Done with bin \n\n\n\n";
}

for ($b = 0; $b < $#bins; $b++ ) {
  print OUT "bin: $b\n";
  @refKeys = common::OrderedKeys(\%{$bins[$b]});
  foreach $ref (@refKeys) {
    for ($l = 0; $l <= $#{$bins[$b]{$ref}}; $l++) {
      print OUT "  $bins[$b]{$ref}[$l][0] $bins[$b]{$ref}[$l][1] $bins[$b]{$ref}[$l][2]\n";
    }
  }
}

sub FindOverlappingIntervals {
  my ($pairs,
      $refSpecies, $refStart, $refEnd,
      $seqLengths,
      $ovpRatio, $overlapping, $toCheck) = @_;
  # Look for overlaps from the reference species.  The coordinates
  # in the ref overlap with other ref alignment coordinates against
  # other species.
  my $refRevStart = $$seqLengths{$refSpecies} - $refEnd;
  my $refRevEnd   = $$seqLengths{$refSpecies} - $refStart;

  print "checking ref: $refSpecies\n";
  if (exists $$pairs{$refSpecies} ) {
    @qrySpecies = common::OrderedKeys(\%{$$pairs{$refSpecies}});
    foreach $qry (@qrySpecies) {
      print "$refSpecies vs $qry \n";
      $deletedQry = 0;
      # look to see if refSpecies overlaps query with the window refStart ... refEnd
      for ( my $inv = 0; $deletedQry == 0 and
	    $inv <= $#{$$pairs{$refSpecies}{$qry}}; ) {
	my $invStart = $$pairs{$refSpecies}{$qry}[$inv][1];
	my $invEnd   = $$pairs{$refSpecies}{$qry}[$inv][2];
	my $qryInvSpec  = $$pairs{$refSpecies}{$qry}[$inv][3];
	my $qryInvStart = $$pairs{$refSpecies}{$qry}[$inv][4];
	my $qryInvEnd   = $$pairs{$refSpecies}{$qry}[$inv][5];
	my $ovpLen = 0;
	if ($refStart <= $invStart) {
	  if ($refEnd <= $invStart) {
	    $ovpLen = 0;
	  } else {
	    $ovpLen = min($refEnd, $invEnd) - $invStart;
	  }
	}
	else {
	  if ($refStart >= $invEnd ) {
	    $ovpLen = 0;
	  } else {
	    $ovpLen = min($refEnd, $invEnd) - $refStart;
	  }
	}
	if ($ovpLen > 0 and 
	    (($ovpLen / abs($refEnd - $refStart + 1) > $ovpRatio) and
	     ($ovpLen / abs($invEnd - $invStart + 1) > $ovpRatio))) {
	  print "accept ovplen: $ovpLen $invStart $invEnd\n";
	  push @{$overlapping}, [$refSpecies, $invStart, $invEnd];
	  push @{$overlapping}, [$qryInvSpec, -1*$qryInvStart, -1*$qryInvEnd];
	  if (! exists $$seqLengths{$qryInvSpec}) {
	    print "could not find length for $qryInvSpec\n";
	    exit(0);
	  }
	  push @{$toCheck}, [$qryInvSpec, 
			     $$seqLengths{$qryInvSpec} - $qryInvEnd,
			     $$seqLengths{$qryInvSpec} - $qryInvStart ];
	    my $ovpLen = 0;
	  print "splicing @{$$pairs{$refSpecies}{$qry}[$inv]}\n";
	  splice @{$$pairs{$refSpecies}{$qry}}, $inv, 1;
	  if ($#{$$pairs{$refSpecies}{$qry}} < 0 ) {
	      print "DELETING pair $refSpecies $qry\n";
	      delete $$pairs{$refSpecies}{$qry};
	      $deletedQry = 1;
	    }
	}
	else {
	  if ($ovpLen > 0 ) {
	    my $refLen = abs($refEnd - $refStart + 1);
	    my $qryLen = abs($qryEnd - $qryStart + 1);
	    my $invLen = abs($invEnd - $invStart + 1);
	    my $refRatio = $ovpLen / $refLen;
	    my $qryRatio = $ovpLen / $qryLen;
	    my $invRatio = $ovpLen / $invLen;
	    print "reject ovplen: $ovpLen / $refLen = $refRatio $ovpLen / ";
	    print "$qryLen = $qryRatio  $ovpLen / $invLen = $invRatio\n";
	    print "$invStart $invEnd $qryStart $qryEnd\n";
	  }
	  ++$inv;
	}
      }
    }
    @qrySpecies = common::OrderedKeys(\%{$$pairs{$refSpecies}});
    if ($#qrySpecies < 0) {
      print "DELETING REFERENCE!!! $refKey\n";
      delete $$pairs{$refSpecies};
    }
  }

  # try to find a query overlap.  Remember, start is after end in 
  # the query coordinates.

  @refKeys = common::OrderedKeys(\%{$pairs});
  $qryStart = $$seqLengths{$refSpecies} - $refEnd + 1;
  $qryEnd   = $$seqLengths{$refSpecies} - $refStart + 1;
  print "checking query overlap. $qryStart $qryEnd \n";
  foreach $refKey (@refKeys) {
    if (exists $$pairs{$refKey}{$refSpecies} ) {
      $deletedQry = 0;
      for ( my $inv = 0; $deletedQry == 0 and
	    $inv <= $#{$$pairs{$refKey}{$refSpecies}}; ) {
	my $invStart = $$pairs{$refKey}{$refSpecies}[$inv][4];
	my $invEnd   = $$pairs{$refKey}{$refSpecies}[$inv][5];
	my $refInvSpec  = $$pairs{$refKey}{$refSpecies}[$inv][0];
	my $refInvStart = $$pairs{$refKey}{$refSpecies}[$inv][1];
	my $refInvEnd   = $$pairs{$refKey}{$refSpecies}[$inv][2];
	my $ovpLen = 0;
	if ($qryStart <= $invStart) {
	  if ($qryEnd <= $invStart) {
	    $ovpLen = 0;
	  } else {
	    $ovpLen = min($qryEnd, $invEnd) - $invStart;
	  }
	} else {
	  if ($qryStart >= $invEnd ) {
	    $ovpLen = 0;
	  } else {
	    $ovpLen = min($qryEnd, $invEnd) - $qryStart;
	  }
	}
	if ($ovpLen > 0 and
	    (($ovpLen / abs($qryEnd - $qryStart + 1) > $ovpRatio) and
	     ($ovpLen / abs($invEnd - $invStart + 1) > $ovpRatio))) {
	  print "found query overlap: $refKey $refSpecies $qryStart $qryEnd $invStart $invEnd\n";
	  push @{$overlapping}, [$refSpecies, -1*$invStart, -1*$invEnd];
	  push @{$overlapping}, [$refInvSpec, $refInvStart, $refInvEnd];
	  my $refStartRev = $$seqLengths{$refInvSpec} - $refInvEnd;
	  my $refEndRev   = $$seqLengths{$refInvSpec} - $refInvStart;
	  push @{$toCheck}, [$refInvSpec, $refInvStart, $refInvEnd];
	  print "splicing @{$$pairs{$refKey}{$refSpecies}[$inv]}\n";
	  splice @{$$pairs{$refKey}{$refSpecies}}, $inv, 1;
	  if ($#{$$pairs{$refKey}{$refSpecies}} < 0 ) {
	    print "DELETING pair $refKey $refSpecies\n";
	    delete $$pairs{$refKey}{$refSpecies};
	    $deletedQry = 1;
	  }
	}
	else {
	  if ($ovpLen > 0) {
	    my $qryLen = abs($qryEnd - $qryStart + 1);
	    my $invLen = abs($invEnd - $invStart + 1);
	    my $qryRatio = $ovpLen / $qryLen;
	    my $invRatio = $ovpLen / $invLen;
	    print "ovplen: $ovpLen / $qryLen = $qryRatio  $ovpLen / $invLen = $invRatio\n";
	  }
	  ++$inv;
	}
      }
    }
    @qrySpeciesList = common::OrderedKeys(\%{$$pairs{$refKey}});
    if ($#qrySpeciesList < 0) {
      print "DELETING REFERENCE!!! $refKey\n";
      delete $$pairs{$refKey};
    }
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
