#!/usr/bin/env perl

use zooproject;

if ($#ARGV < 0) {
  print "usage: PrintSpeciesIds.pl targettablename nametablename [-warn]\n";
  exit(0);
}

if ($#ARGV == 0) {
  `cat @ARGV[0]`;
  exit(0);
}

$targetTableName = shift @ARGV;
if ($#ARGV >= 0) {
  $nameTableName = shift @ARGV;
}
$nameTableName = $ENV{"EUSRC"} . "/zooproject/zoonames.txt";
if ($#ARGV >= 0) {
  $nameTableName = shift @ARGV;
}
$warn = 0;
if ($#ARGV >= 0) {
  $opt = shift @ARGV;
  if ($opt eq "-warn") {
    $warn = 1;
  }
}
my @targetTable;
zooproject::ParseTargetTable($targetTableName, \@targetTable);


# Read in name ids
my %niscToName;
my %niscToId;
my %nameToId;
my @idToName;
zooproject::ParseNameTable($nameTableName,
			   \%niscToName, \%niscToId, \%nameToId, \@idToName);

$clone = 0;
for ($clone = 0; $clone <= $#targetTable - 1; $clone++ ) {
  if ($targetTable[$clone]{"Organism"} ne $targetTable[$clone+1]{"Organism"}) {
    if (exists $niscToId{$targetTable[$clone]{"Organism"}}) {
      $id = $niscToId{$targetTable[$clone]{"Organism"}};
      print "$id ";
    }
    else {
      if ($warn) {
	$org = $targetTable[$clone]{"Organism"};
	print "did not find $org\n";
      }
    }
  }
}
