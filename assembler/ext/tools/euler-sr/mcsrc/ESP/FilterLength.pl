#!/usr/bin/env perl

if ($#ARGV < 1) {
  print "\nusage: $0 mapfile filteredMapFile [lengthThreshold]\n";
  print "    removes matches with length less than lengthThreshold\n";
  print "    lengthThreshold (300)\n\n";
  exit(0);
}
$mapfile = shift @ARGV;
$outfile = shift @ARGV;
$lengthThreshold = 300;
if ($#ARGV >=0) {
  $lengthThreshold = shift @ARGV;
}

open(MAP, $mapfile) or die "cannot open $mapfile\n";
open(OUT, ">$outfile") or die "cannot write to $outfile\n";


$first = 1;
while(<MAP>) {
  $line = $_;
  if ($line =~ /^>(.*)/) {
    # reset the hits list
    $prevTitle = $title;
    $title = $1;
    if ($first == 0) {
      # only print this if it is a quality match
      if ($#hits >= 0) {
	print OUT ">$prevTitle\n";
	print OUT "@hits";
      }
    }
    $first = 0;
    @hits = ();
  }
  else {
    @vals = split (/\s+/, $line);
    $length = abs($vals[2] - $vals[1]);
    if ( $length > $lengthThreshold) {
      push @hits, $line;
    }
  }
}


# process the last
if ($first == 0) {
  if ($#hits >= 0) {
    print ONEMAP ">$prevTitle\n";
    print ONEMAP "@hits\n";
  }
}
