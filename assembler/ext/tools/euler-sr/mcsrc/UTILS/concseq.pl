#!/usr/bin/perl
use strict;

my $destbase = shift(@ARGV);
my $destName     = "$destbase.corr";
my $rulefileName = "$destbase.rule";
my $matesName    = "$destbase.mates";
my $locsName     = "$destbase.loc";

my @seqfiles  = ();
my %clonelens = {};
my %clonevars = {};
my %readlocs  = {};

my $arg = "";
my $readLoc = 0;
while ($#ARGV >= 0) {
  $arg = shift(@ARGV);
  if ($arg =~ /\-/) {
    if ($arg eq "-d") {
      $readLoc = 1;
    }
  }
  else {
    push @seqfiles, $arg;
    $clonelens{@seqfiles[$#seqfiles]} = shift(@ARGV);
    $clonevars{@seqfiles[$#seqfiles]} = shift(@ARGV);
  }
}

print "opening mane file for writing\n";
open (DESTSEQ, ">$destName") or die "cannot open $destName for writing\n";
open (RULE, ">$rulefileName") or die "cannot open $rulefileName for writing\n";
open (MATES, ">$matesName") or die "cannot open $matesName for writing\n";
if ($readLoc) {
  open(LOCS, ">$locsName")or die "cannot open $locsName for writing\n";
}
print "done openiing files for writing\n";
print RULE <<'_ruleprologue';
/*			Names are: Library_Plate_well.primer	
				primer  library plate   length range	*/

Single reads:                   s       ALL     ALL
Single reads:                   p       ALL     ALL
_ruleprologue


my $seqindex  = 0;
my $mateindex = 0;
my $set   = 0;
my $seqoffset  = 0;
my $mateoffset = 0;
my $index = 0;
my ($low, $high);
my $matenum;
my $matea;
my $mateb;  
my $line;
my $seq;
foreach $seq (@seqfiles) {
  print "working on seq $seq\n";
  open (SEQ, "$seq.corr") || die "cannot open seq file $seq.corr\n";
  $low = $clonelens{$seq} - $clonevars{$seq};
  $high = $clonelens{$seq} + $clonevars{$seq};
  print RULE "Double-barreled reads:          f r     ALL	2-2     $low    $high\n";
  while (<SEQ>) {
    $line = $_;
    chomp($line);
    if ($line =~ /^>(\D+)\d+/) {
      print DESTSEQ ">$1$index\n";
      $index = $index + 1;
    }
    else 
      { print DESTSEQ "$line\n"; }
  }
  close(SEQ);
  open (MATE, "$seq.mates") || die "cannot open mates file $seq.mates\n";
  while (<MATE>) {
    /(\D+)(\d+) (\d+) (\d+)/;
    $matenum = $2 + $mateoffset;
    $matea   = $3 + $seqoffset;
    $mateb   = $4 + $seqoffset;
    print MATES "$1$matenum $matea $mateb $set\n";
    $mateindex = $mateindex + 1;
  }

  if ($readLoc) {
    open(LOC, "$seq.loc") || die "cannot open locations file $seq.len\n";
    while (<LOC>) {
      print LOCS $_;
    }
    close(LOC);
  }

  $seqoffset  = $index;
  $mateoffset = $mateindex;
  $set        = $set + 1;

}
close RULE;
close MATES;
close DESTSEQ;
if ($readLoc) {
  close(LOCS);
}


