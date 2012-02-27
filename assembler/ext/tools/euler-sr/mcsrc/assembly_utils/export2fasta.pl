#!/usr/bin/env perl
if ($#ARGV < 1) {
		print "usage: export2fasta.pl in.eland out.fasta\n";
		exit(0);
}

$in = shift @ARGV;
$out = shift @ARGV;

open(IN, "$in") or die "cannot open export file $in\n";
open(OUT, ">$out") or die "cannot open fasta file $out\n";
#HWI-EAS258_2_FC204VU		1	1	425	708		1	GTAAGTGATCTTGCTGAGCTCAAAGACAACAAACG	32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 30 32 28 32 32 32 32 	chrX.fa		6814121	R	35	94	0			0N	Y
#HWI-EAS258_2_FC204VU		1	1	371	659		2	TTCGATCTAACATCTTTAAAATTTGTAGAACTATG	4 4 4 4 4 4 4 4 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 	chrIII.fa		8870323	R	10A24	7	107			-196	F	Y

while(<IN>) {
		if ($_ =~ /(\S+)\s+(\d+)\s+(\d+)\s+([\-]?\d+)\s+([\-]?\d+)\s+(\d+)\s+(\w+)\s.*/) {
				print  OUT ">$1:$2:$3:$4:$5/$6\n$7\n";
		}
		else {
				print "ERROR parsing line $_";
				exit(1);
		}
}
		
