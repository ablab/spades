#!/usr/bin/env perl
if ($#ARGV < 1) {
		print "usage: PrintEdgeIntvs.pl intvFile edgeNumber\n";
		exit(1);
}

$intvFile = shift @ARGV;
$edgeNumber = shift @ARGV;

open(IF, $intvFile) or die "cannot open $intvFile\n";

$printIntv = 0;
while(<IF>) {
		$il = $_;
		if ($il =~ /EDGE (\d+) .*/) {
				if ($edgeNumber == $1) {
						$printIntv = 1;
				}
				else {
						if ($printIntv == 1) {
								exit(0);
						}
				}
		}
		else {
				if ($printIntv) {
						$il =~ /INTV (\d+) .*/;
						print "$1\n";
				}
		}
}
