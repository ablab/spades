#!/usr/bin/env perl
if ($#ARGV < 1) {
		print "usage: Trim454Reads.pl in out\n";
		exit(0);
}
$in = shift @ARGV;
$out = shift @ARGV;
open(IN, $in) or die "cannot open $in\n";
open(OUT, ">$out") or die "cannot open $out\n";
$trimFront = 0;
while ($#ARGV >= 0) {
		$opt = shift @ARGV;
		if ($opt eq "-trimFront") {
				$trimFront = shift @ARGV;
		}
}
print "trim: $trimFront\n";
$curTitle = "";
while(<IN>) {
		if ($_ =~/^>/) {
				$curTitle = $_;

				# process the previous read.
				$readSeq = substr($readSeq, $trimFront);
				if ($readSeq =~ /[^N]N$/) {
						$readSeq =~ /(.*[^N])([N]+)$/;
						print "$readSeq\n";
						$readSeq = $1;
						print "$readSeq\n";
				}
        print OUT "$prevTitle";
				print OUT "$readSeq\n";
				
				$readSeq = "";
				$prevTitle = $curTitle;
		}
		else {
				chomp $_;
				$readSeq .= $_;
		}
}
if ($curTitle ne "") {
		if ($readSeq =~ /[^N]N$/) {
				$readSeq =~ /(.*[^N])[N]+$/;
				$readSeq = $1;
		}
		print OUT "$curTitle";
		print OUT "$readSeq\n";
}
