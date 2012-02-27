#!/usr/bin/env perl
if ($#ARGV < 1) {
		print "usage: eland2fasta.pl in.eland out.fasta\n";
		exit(0);
}

$in = shift @ARGV;
$out = shift @ARGV;

open(IN, "$in") or die "cannot open eland file $in\n";
open(OUT, ">$out") or die "cannot open fasta file $out\n";

while(<IN>) {
		if ($_ =~ />(\S+)\s+(\w+)/) {
				print OUT ">$1\n$2\n";
		}
		else {
				print "ERROR parsing line $_";
				exit(1);
		}
}
		
