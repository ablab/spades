#!/usr/bin/env perl

use ReadLibrary;
use ReadMatch;

if ($#ARGV <  5) {
  print "\nusage: $0 reads locations mate-names start end outfile [-b blacklist] [-nomate]\n\n";
  print "prints the titles of the reads that are mapped between (start, end) of the \n";
  print "reference genome.\n";
  print "reads      The fasta file with reads\n";
  print "locations  The locations file fo format \n";
  print "  'mate1pos mate1dir mate2pos mate2dir type mate1index mate2index'\n";
  print "mate-names The names of the mate-pairs, two names per line\n";
  print "start,     end are the start and ending positions in the reference\n";
  print "           genome to extract reads\n";
  print "blacklist  The blacklist file.  Exclude any read mentioned in this\n";
  print "           file, and its mate\n\n";
  print "-nomate    Only extract sequences that match, not the sequence and the mate\n";
  exit(0);
}

$readsFile = shift @ARGV;
$locationsFile = shift @ARGV;
$namesFile = shift @ARGV;
$start = shift @ARGV;
$end   = shift @ARGV;
$outfile= shift @ARGV;
$extractMate =  1;
while ($#ARGV >= 0) {
  $opt = shift @ARGV;
  if ($opt eq "-b" ) {
    $blacklistFile = shift @ARGV;
  }
  if ($opt eq "-nomate") {
    $extractMate = 0;
  }
}


open(READS, "$readsFile") or die "cannot open $infile\n";
open(NAM, "$namesFile") or die "cannot open $namesFile\n";
open(OUT, ">$outfile") or die "cannot open $outfile\n";


# parse the reads to ignore
%blacklist = ();
if ($blacklistFile != "") {
  open(BF, "$blacklistFile") or die "cannot open $blacklistFile\n";
  while (<BF>) {
    $line = $_;
    chomp $line;
    if ($line =~ />/) {
      ($title, $nam, $type, $name) = ReadLibrary::ParseABITitle($line);
    }
    else {
      $name = $line;
    }
    $blacklist{$name} = 1;
  }
  print "parsed $blacklistFile\n";
}


# parse the names file. this has to be done first
# to use the name as a key for the locations
@names = ();
%mates  = ();
while (<NAM>) {
  $line = $_;
  chomp $line;
  @namePair = split(/\s+/, $line);
  push @names, [@namePair];
  $mates{$namePair[0]} = $namePair[1];
  $mates{$namePair[1]} = $namePair[0];
}

@nmates = values %mates;
print "read: $#nmates mates\n";

# parse the locations of the reads
%locations = ();
ReadLibrary::ReadLocationsFile($locationsFile, \%locations);
@loc = values %locations;
print "parsed $locationsFile $#loc \n";


# now read in the reads, and decide if each
# is 
#   1. Either end of the mate-pair is in the interval 
#      specified.
#   2. Neither end of the mate-pair is blacklisted.
#

$print = 0;

print "indexing iwth $ReadLibrary::MQryStart $ReadLibrary::MQryEnd\n";
%keep = ();
while(<READS>) {
  $line = $_;
  if ($line =~ /^>(.*)/) {
    $title = $1;
    ($base,$dir,$type,$name) = ReadLibrary::ParseABITitle($line);
    $mate = $mates{$name};
    $nmatches = scalar @{$locations{$name}};
    $m = 0;
    if ($extractMate) {
      if (exists $mates{$name} and $nmatches > 0) {
	if (not exists $blacklist{$name} and
	    not exists $blacklist{$mate}) {
	  #print "checking $locations{$name}[0][$ReadLibrary::MQryStart]
	  # $locations{$name}[0][$ReadLibrary::MQryEnd] \n";
	  if (($locations{$name}[$m]->{queryStart} >= $start and
	       $locations{$name}[$m]->{queryStart} <= $end) or
	      ($locations{$name}[$m]->{queryEnd} >= $start and
	       $locations{$name}[$m]->{queryEnd} <= $end)) {
	    $keep{$name} = 1;
	    $keep{$mate} = 1;
	  }
	}
      }
    }
    else {
      if (not exists $blacklist{$name}) {
	#print "checking $locations{$name}[0][$ReadLibrary::MQryStart]
	# $locations{$name}[0][$ReadLibrary::MQryEnd] \n";
	if (($locations{$name}[$m]->{queryStart} >= $start and
	     $locations{$name}[$m]->{queryStart} <= $end) or
	    ($locations{$name}[$m]->{queryEnd} >= $start and
	     $locations{$name}[$m]->{queryEnd} <= $end)) {
	  $keep{$name} = 1;
	}
      }
    }
  }
}

close READS;
open(READS, "$readsFile") or die "cannot open $infile\n";
while(<READS>) {
  $line = $_;
  if ($line =~ /^>(.*)/) {
    $print = 0;
    $title = $1;
    ($base,$dir,$type,$name) = ReadLibrary::ParseABITitle($line);
    if ($keep{$name} == 1) {
      $print = 1;
      print OUT "$line";
    }
  }
  else {
    if ($print == 1) {
      print OUT "$line";
    }
  }
}

close OUT;
