#!/usr/bin/env perl
use strict;

if ($#ARGV != 1) { print "usage: $0 reads overlaplen \n"; exit(-1);}
open(READSFILE, "@ARGV[0]") or die "cannot open @ARGV[0]\n";
my $overlaplen = @ARGV[1];

my ($title, $line, $read, $title2);

#Get the reads file along with all the lengths.
my %lengths = ();
my $len;
$title = <READSFILE>;
$title =~ /> (.*)/;
$title = $1;
$read = "";
while (<READSFILE>) {
  chomp;
  $line = $_;
  if ($line =~ />/) {
    $len = length($read);
    $lengths{$title} = $len;
    $line =~ /> (.*)/;
    $title = $1;
    $read = "";
  }
  else {
    $read = $read . $line;
  }
}
close(READSFILE);

<STDIN>;
<STDIN>;
<STDIN>;
my @fields = ();
my ($lena, $lenb);
my $numOverlaps = 0;
my $isreverse = 0;
while (<STDIN>) {
  @fields = split /\s+/;
  if (@fields[2] != "100.00") {next;}
  if (@fields[4] != 0) { next;}
  if (@fields[0] eq @fields[1]) { next;}
  if (@fields[3] < $overlaplen) { #print "skipping @fields[0] @fields[1] too short @fields[3]\n"; 
    next;}

  if (@fields[8] < @fields[9]) { $isreverse = 0;} else { $isreverse = 1; }
  $lena = $lengths{@fields[0]};
  $lenb = $lengths{@fields[1]};
  
  if (@fields[6] == 1) {
    # look to see if overlap is completely encompassing.
    if (@fields[6] + @fields[3] - 1 == @fields[7]) { $numOverlaps = $numOverlaps + 1; next; }
    # forward overlaps wit beginning
    if (!$isreverse) {
      if (@fields[9] == $lenb) { $numOverlaps = $numOverlaps + 1; next; }
      else {
	# otherwise problem
#	 print "problem with forward begin align $lena $lenb\n";
#	 print $_;
	next;
      }
    }
    else {
      # reverse overlaps with beginning
      if (@fields[9] == 1) { $numOverlaps = $numOverlaps + 1; }
      else {
#	 print "problem with reverse begin align $lena $lenb\n";
#	 print $_;
	next;
      }
    }
  }
  elsif(@fields[7] == $lena) {
    if (!$isreverse) {
      if (@fields[8] == 1) { $numOverlaps = $numOverlaps + 1; }
      else {
#	 print "problem with forward end align $lena $lenb\n";
#	 print $_;
      }
    }
    else {
      if (@fields[8] == $lenb) { $numOverlaps = $numOverlaps + 1; }
      else {
#	 print "problem with reverse end align $lena $lenb\n";
#	 print $_;
      }
    }
  }
  elsif (!$isreverse && (@fields[8] + @fields[3] == @fields[9])) { $numOverlaps = $numOverlaps + 1; }
  elsif ($isreverse && (@fields[9] + @fields[3] == @fields[8])) { $numOverlaps = $numOverlaps + 1; }
  else {
#    print "unclassifiable problem $lena $lenb\n";
#    print "$_";
  }
}

print "got $numOverlaps overlaps\n";
