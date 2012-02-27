#!/usr/bin/env perl

if ($#ARGV < 1) {
  print "usage: $0 mapfile outbase [lengthThreshold]\n";
  print "    splits mapfile into 3 files 0 matches, 1 high-scoring match, and > 1 high scoring match\n";
  print "    lengthThreshold (300) The minimal length HSP to consider a hit\n";
  exit(0);
}
$mapfile = shift @ARGV;
$outbase = shift @ARGV;
$lengthThreshold = 300;
if ($#ARGV >=0) {
  $lengthThreshold = shift @ARGV;
}

open(MAP, $mapfile) or die "cannot open $mapfile\n";
open(NOMAP, ">$outbase.none") or die "cannot write to $outbase.none\n";
open(ONEMAP, ">$outbase.one") or die "cannot write to $outbase.one\n";
open(MULTIMAP, ">$outbase.multi") or die "cannot write to $outbase.multi\n";

$first = 1;
while(<MAP>) {
  $line = $_;
  if ($line =~ /^>(.*)/) {
    # reset the hits list
#    print "creating new hits\n";
    $prevTitle = $title;
    $title = $1;
    if ($first == 0) {
      if ($#hits < 0) {
	print NOMAP ">$prevTitle\n";
      }
      elsif ($#hits == 0) {
	print ONEMAP ">$prevTitle\n";
	print ONEMAP "@hits";
      }
      else {
	print MULTIMAP ">$prevTitle\n";
	print MULTIMAM @hits;
      }
    }
    $first = 0;
    @hits = ();
  }
  else {
    push @hits, $line;
  }
}
if ($first == 0) {
  if ($#hits < 0) {
    print NOMAP ">$prevTitle\n";
  } elsif ($#hits == 0) {
    print ONEMAP ">$prevTitle\n";
    print ONEMAP "@hits";
  } else {
    print MULTIMAP ">$prevTitle\n";
    print MULTIMAM @hits;
  }
}
