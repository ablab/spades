#!/usr/bin/env perl
if ($#ARGV < 0) {
		print "usage: $0 fasta_file\n";
		print "fasta_file has the read coordinates stored as pos=NNN dir=N\n";
		print "dir is 0 (forward) or 1 (reverse)\n";
}
open(IN, "$ARGV[0]") or die "cannot open $ARGV[0]\n";


while(<IN>) {
		if ($_ =~ /pos=\s*(\d+)/) {
				$pos = $1;
				if ($_ =~ /strand=\s*(\d)/) {
						$strand = $1;
						print "$strand $pos\n";
				}
		}
}
