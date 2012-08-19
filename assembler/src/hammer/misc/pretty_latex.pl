#!/usr/bin/perl -w

use strict;

my $ecolinum = 4565344.0;

my $dirname = shift;
opendir DIR, $dirname or die $!;

my @solid; my @solidwrong;
my @bad; my @badwrong;
my @reads; my @readswrong;

my $maxnum = -1;

while (my $file = readdir(DIR)) {
	if ( $file eq "." || $file eq "..") { next; }
	if ( $file =~ /([0-9]+)\.results\.solidkmers/ ) {
		my $num = $1; if ( $num > $maxnum ) $maxnum = $num;
		open FILE, $dirname . "/" . $file or die $!;
		my @l = <FILE>;
		if ( $l[0] =~ /trusted: ([0-9]+)/ ) {
			$solid[$num] = $1;
			print $l[0];
		} else { die "Wrong format at line 1 of $dirname/$file"; }
		if ( $l[1] =~ /erroneous: ([0-9]+)/ ) {
			$solidwrong[$num] = $1;
			print $l[1];
		} else { die "Wrong format at line 2 of $dirname/$file"; }
		if ( $l[2] =~ /bad: ([0-9]+)/ ) {
			$bad[$num] = $1;
			print $l[2];
		} else { die "Wrong format at line 3 of $dirname/$file"; }
		if ( $l[3] =~ /erroneous: ([0-9]+)/ ) {
			$badwrong[$num] = $1;
			print $l[3];
		} else { die "Wrong format at line 4 of $dirname/$file"; }
		close FILE;
	}
	if ( $file =~ /([0-9]+)\.results\.newreads/ ) {
		my $num = $1; if ( $num > $maxnum ) $maxnum = $num;
		open FILE, $dirname . "/" . $file or die $!;
		my @l = <FILE>;
		if ( $l[0] =~ /trusted: ([0-9]+)/ ) {
			$reads[$num] = $1;
			print $l[0];
		} else { die "Wrong format at line 1 of $dirname/$file"; }
		if ( $l[1] =~ /erroneous: ([0-9]+)/ ) {
			$readswrong[$num] = $1;
			print $l[1];
		} else { die "Wrong format at line 2 of $dirname/$file"; }
		close FILE;
	}
}

for ( my $i=0; $i < $maxnum; ++$i ) {
print ($i+1) . " & " . $solid[$num] . " & " . ($solid[$num] - $solidwrong[$num]) . " & " . $bad[$num] . " & " . $badwrong[$num];
print " & " . $reads[$num] . " & " . ($reads[$num] - $readswrong[$num]) . " & " . ($reads[$num] - $readswrong[$num])/$ecolinum . " \\ \n";
}


