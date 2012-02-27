#!/usr/bin/env perl

if ($#ARGV != 1) {
		print "usage: JoinReads.pl inFile outFile\n";
		print "       Concatenates all reads in 'inFile' into a single\n";
		print "       lined outfile\n";
		exit(1);
}
$in = shift @ARGV;
$out = shift @ARGV;


open(IN, $in) or die "cannot open $in\n";
open(OUT, ">$out") or die "cannot open $out\n";

print OUT ">joined\n";
# get the first title
<IN>;
while(<IN>) {
		if ($_ =~ /^>/) {

				print OUT "NNN";
				print OUT $readSeqn;
				$readSeqn = "";
		}
		else {
				chomp $_;
				$readSeqn .= $_;
		}
}
