#!/usr/bin/env perl
if ($#ARGV < 1) {
  print "usage: AlignmentToFasta.pl in.aln out.fasta\n";
  exit(1);
}
$in  = shift @ARGV;
$out = shift @ARGV;
open(IN, "$in") or die "cannot open $in\n";
open(OUT, ">$out") or die "cannot open $out\n";
#66386   0       000000000000000000000000000000 0        GTTTTCCCTCAGCTACTATATTTCTGTTCT

@profile = ();
$readIndex = 0;
while (<IN>) {
		$_ =~ /(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\S+)/;
		$pos = $1;
		$strandNum = $2;
		if ($strandNum == 0) {
				$strand = "FOR";
    }
		else {
				$strand = "REV";
		}
		$title = ">". $readIndex. "_" . $pos  . "_" . $strand . " pos=$pos strand=$strandNum";
		print OUT "$title\n";
		print OUT "$5\n";
		$readIndex++;
}

#$rate = $nErrors/$total;
#print "$nErrors $total $rate\n";
close OUT;
