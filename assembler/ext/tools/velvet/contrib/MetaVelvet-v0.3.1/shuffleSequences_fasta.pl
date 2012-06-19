#!/usr/bin/perl

if (!@ARGV) {
	print "Usage: $0 forward_reads.fa reverse_reaads.fa outfile.fa\n";
	print "\tforward_reads.fa / reverse_reads.fa : paired reads to be merged\n";
	print "\toutfile.fa : outfile to be created\n";
	system.exit(0);	
}

$filenameA = $ARGV[0];
$filenameB = $ARGV[1];
$filenameOut = $ARGV[2];

die "Could not open $filenameA" unless (-e $filenameA);
die "Could not open $filenameB" unless (-e $filenameB);

open FILEA, "< $filenameA";
open FILEB, "< $filenameB";

open OUTFILE, "> $filenameOut";

my ($lineA, $lineB);

$lineA = <FILEA>;
$lineB = <FILEB>;

while(defined $lineA) {
	print OUTFILE $lineA;
	$lineA = <FILEA>;
	while (defined $lineA && $lineA !~ m/>/) { 
		print OUTFILE $lineA;
		$lineA = <FILEA>;
	}

	print OUTFILE $lineB;
	$lineB = <FILEB>;
	while (defined $lineB && $lineB !~ m/>/) { 
		print OUTFILE $lineB;
		$lineB = <FILEB>;
	}
}
