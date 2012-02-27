#!/usr/bin/env perl
if ($#< = 2) {
  print "usage: TrimReads.pl in out length\n";
  exit(0);
}
$in = shift @ARGV;
$out = shift @ARGV;
$length = shift @ARGV;

open(IN, "$in") or die "cannot open $in\n";
open(OUT, ">$out") or die "cannot open $out\n";
while (<IN>) {
	#print title
	print OUT $_;
	$read = <IN>;
	$truncRead = substr($read, 0, $length);
	print OUT "$read\n";
}
	
