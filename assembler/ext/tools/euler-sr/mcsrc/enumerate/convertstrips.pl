#!/usr/bin/env perl

if ($#ARGV != 1) {
  print "usage: $0 infile outfile\n";
  exit(0);
}

$infile = shift @ARGV;
$outfile = shift @ARGV;
open (IN, "$infile") or die "cannot open $infile\n";
open (OUT, ">$outfile") or die "cannot open $outfile\n";

$total = 0;
@lines = <IN>;
close IN;
$count = $#lines + 1;
print OUT "$count\n";
$i = 0;
for $i (0 .. $#lines) {
  @lines[$i] =~ /\d+\s([\S].*)/;
  print OUT "$1\n";
}

close OUT;

  
