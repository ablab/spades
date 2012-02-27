#!/usr/bin/env perl

$hist = {};
while (<>) {
		$val = $_;
		chomp $val;
		$hist{$val}++;
}


foreach $val (sort {$a <=> $b } keys %hist) {
		print "$val $hist{$val}\n";
}
