#!/usr/bin/env perl

if ($#ARGV != 1) {
  print "usage: $0 extracts_title_file strips_title_file\n";
  exit(1);
}
$etf = shift @ARGV;
$stf = shift @ARGV;

open (ET, "$etf") or die "cannot open $etf\n";
open (ST, "$stf") or die "cannot open $stf\n";

@etc  = ();
while (<ET>) {
  $etl = $_;
  $pos = index $etl, "|";
  @contigs = get_contigs($etl, $pos);
  push @etc, [ @contigs ];
}

@stc = ();
while (<ST>) {
  $stl = $_;
  $pos = 0;
  @contigs = ();
  @contigs = get_contigs($stl, $pos);
  $contigStr = join(' ', @contigs);
  push @stc, $contigStr;
#  print "just added $contigStr\n";
}

for $i (0 .. $#etc) {
  @idxs = ();
  print "$i: ";
  for $j (0 .. $#{$etc[$i]}) {
    $continue = 1;
    $k = 0;
    while ($k < $#stc && $continue) {
#      print "looking at $etc[$i][$j] \n";
#      print "comparing @stc[$k]\n";
      if (index(@stc[$k], $etc[$i][$j]) >= 0) {
	push @idxs, $k;
	$continue = 0;
      }
      $k++;
    }
  }
  if ($#idxs == -1) {
    print "NOT FOUND: @{$etc}[$i]\n";
  }
  elsif ($#idxs > 0) {
    $start = @idxs[0];
    $gap = 0;
    for $l (1 .. $#idxs) {
      if (@idxs[$l] != $start) {
	print "GAP in @{$etc[$i]} indices: @idxs\n";
	$gap = 1;
      }
    }
    if ($gap == 0) {
      print " passed\n";
    }
  }
  else {
    print " passed\n";
  }
}


sub get_contigs {
  my($str, $pos) = @_;
  $pos = index $str, 'contig', $pos;
  if($pos < 0) {
    return [];
  }
  my($sub) = substr $str, $pos;
  my @contigs = ();
  while ($sub =~ /(contig \d+).*/) {
    $c = $1;
    chomp $c;
    push @contigs, $c;
    $pos = index($str, $c) + length($c);
    $sub = substr $str, $pos;
  }
  return @contigs;
}
