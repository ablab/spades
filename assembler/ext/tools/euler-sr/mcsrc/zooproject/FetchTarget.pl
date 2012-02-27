#!/usr/bin/env perl

use zooproject;
if ($#ARGV < 1) {
  print "usage: FetchTarget.pl targetFileName targetName [-organism organism]\n";
  exit(0);
}

$targetFileName = shift @ARGV;
$targetName  = shift @ARGV;
$src = $ENV{"EUSRC"};
$nameTableFileName = "$src/zooproject/nametable.txt";
$organism = "";
while ($#ARGV >= 0) {
  $option = shift @ARGV;
  if ($option eq "-organism") {
    $organism = shift @ARGV;
  }
}
# get map from nisc species names to species unique ids
my %niscToName;
my %niscToId;
my %nameToId;
my @idToName;

# Read in name ids
zooproject::ParseNameTable($nameTableFileName,
			   \%niscToName, \%niscToId, \%nameToId, \@idToName);


# Read in clone
my @cloneInfo;
zooproject::ParseTargetTable($targetFileName, \@cloneInfo);


# Iterate over target lines and output 

$ci = scalar @cloneInfo;
print "got $ci clone info\n";

while ($#cloneInfo >= 0) {
  # pick the first species, and find all others that match it
  $clone = $cloneInfo[0];
  @cloneIds = ();
  $cloneUID = $niscToID[$$clone{"Organism"}];
  $cloneOrganism = $niscToName{$$clone{"Organism"}}; 
  print "comparing $organism $cloneOrganism\n";
  if ($organism eq "" or $organism eq $cloneOrganism) {
    while ($#cloneInfo >= 0 and 
	   $cloneInfo[0]{"Organism"} eq $$clone{"Organism"}) {
      if ($cloneInfo[0]{"Status"} eq "Sequenced") {
	push @cloneIds, $cloneInfo[0]{"GenbankId"};
      }
      shift @cloneInfo;
    }
    if ($#cloneIds >= 0) {
      $specName = $$clone{"Organism"};
      #    print "species $specName has ids @cloneIds\n";
      $species = $niscToName{$$clone{"Organism"}};
      $id      = $niscToId{$$clone{"Organism"}};
      ($numFetched, $numTotal) = zooproject::FetchSequences($targetName, $id, \@cloneIds);
      print "summary: $species $id, $numFetched / $numTotal\n";
    }
  }
  else {
    shift @cloneInfo;
  }
}
