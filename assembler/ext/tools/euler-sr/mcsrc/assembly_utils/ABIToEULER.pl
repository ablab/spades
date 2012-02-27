#!/usr/bin/env perl


use ReadLibrary;

if ($#ARGV != 2) {
  print "\n";
  print "usage: $0 readsFile matesFile eulerFormatFile\n";
  print "\n";
  print "transform the titles in readsfile so that they adhere \n";
  print "to the strict format imposed by name.rul, that describes \n";
  print "mate pair information for euler.\n";
  print "for paired-end reads, name.rul contains a line: \n";
  print "Double-barreled reads:          x y     ALL     201-400   3000 10000\n";
  print "the corresponding title is: LIBRARY_PLATE_WELL.primer\n";
  print "LIBRARY_PLATE_WELL for paired-ends should be the same, and primer \n";
  print "matesfile is in the format :\n";
  print "readname1 readname2\n";
  print "should be x for front, and y for reverse\n\n";

  exit(0);
}

$readsFile = shift @ARGV;
$matesFile = shift @ARGV;
$outFile   = shift @ARGV;

open (RF, "$readsFile") or die "cannot open $readsFile\n";
open (MF, "$matesFile") or die "cannot open $matesFile\n";
open (OF, ">$outFile") or die "cannot open $outFile\n";

# read in mate-pair information
%mates = ();
while (<MF>) {
  $line = $_;
  chomp $line;
  @names = split(/\s+/, $line);
  $mates{$names[0]} = $names[1];
  $mates{$names[1]} = $names[0];
}

%mateIndices = ();
$mateIndex = 0;
$unpairedIndex = 0;
# read in the read informatoin
while (<RF>) {
  $line = $_;
  if ($line =~ />/) {
    ($base, $dir, $type, $name ) = ReadLibrary::ParseABITitle($line);
    if (not exists $mates{$name}) {
      # this is an unpaired read. No matter what the type is, 
      # treat it as unpaired
      $library = "unpaired";
      $plate   = "1";
      $well    = $unpairedIndex; $unpairedIndex++;
      $primer  = "s";
    }
    else {
      if ($dir  eq "F") {
	$primer = "x";
      }
      elsif ($dir eq "R") {
	$primer = "y";
      }
      else {
	print "error, direction '$dir' is not defined\n";
      }
      if ($type == 0) {
	# fosmid
	$library = "fosmid";
	$plate   = "0";
      }
      elsif ($type == 1) {
	$library = "cosmid";
	$plate = "1";
      }
      elsif ($type == 3) {
	$library = "fosmid";
	$plate = "3";
      }

      if (exists $mateIndices{$base}) {
	$well = $mateIndices{$base};
      }
      else {
	$well = $mateIndex;
	$mateIndices{$base} = $mateIndex;
	$mateIndex++;
      }
    }
    print OF ">$library\_$plate\_$well.$primer\n";
  }
  else {
    print OF "$line";
  }
}

