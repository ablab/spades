#!/usr/bin/env perl


$in = shift @ARGV;
$len = shift @ARGV;

open(IN, $in) or die "cannot open $in\n";
while(<IN>) {
		if (/^>/) {
				print;
		}
		else {
				$read = substr($_, 0, $len);
				print "$read\n";
		}
}

