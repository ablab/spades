#!/usr/bin/env perl

if ($#ARGV != 1) {
  print " usage: $0 pattern outfile\n";
  exit(0);
}
$filepattern = shift @ARGV;
$outfile = shift @ARGV;


@files = glob($filepattern);
%fh = ();
foreach $file (@files) {
  $file =~ /.*\.(\d+)\..*/;
  $number = $1;
  print "got #: $number\n";
  $fh{$number} = $file;
}

@ks = sort {$a <=> $b} keys %fh;

open(OUT, ">$outfile") or die "cannot open $outfile\n";

foreach $k (@ks) {
  open (IN, "$fh{$k}");
  @l = <IN>;
  shift @l;
  print OUT "bin: $k\n";
  print OUT @l;
}
