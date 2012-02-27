#!/usr/bin/env perl

$in  = shift @ARGV;
open(IN, "$in") or die "cannot open $in\n";
#66386   0       000000000000000000000000000000 0        GTTTTCCCTCAGCTACTATATTTCTGTTCT

@profile = ();

while (<IN>) {
		$_ =~ /(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\S+)/;
	  $nErrors += $4;
		$total += length($3);
		if ($4 > 0) {
				print $_;
		}
#		@errorProf = split(//, $3);
#		for ($i = 0; $i <= $#errorProf; $i++) {
#				$profile[$i] += $errorProf[$i];
#		}
}

$rate = $nErrors/$total;
print "$nErrors $total $rate\n";
