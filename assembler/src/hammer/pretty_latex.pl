#!/usr/bin/perl -w

use strict;

my $ecolinum = 4565344;

my $dirname = shift;
opendir DIR, $dirname or die $!;

my @solid; my @solidwrong;
my @bad; my @badwrong;
my @reads; my @readswrong;

while (my $file = readdir(DIR)) {
	if ( $file eq "." || $file eq "..") next;
	if ( $file =~ /([0-9]+)\.results\.solidkmers/ ) {
		my $num = $1;
		open FILE, $dirname . "/" . $file or die $!;
		my $l = <FILE>;
		if ( $l =~ /trusted: ([0-9]+)/ ) {
			$solid[$num] = $1;
		} else { die "Wrong format at line 1 of $dirname/$file"; }
		close FILE;
	}
}

die "OK";

while (my $kmer = <FILE>) {
	my $qual = <FILE>;
	$kmer =~ s/\s+$//;
	print $kmer . " " . $qual;
}

