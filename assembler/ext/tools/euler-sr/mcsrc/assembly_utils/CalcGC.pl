#!/usr/bin/env perl

$ncg = 0;
$total = 0;
while(<>) {
		if (!/^>/) {
				$ncg += $_ =~ tr/CG/cg/;
				$total += length($_);
		}
}

$pctGC = $ncg / $total;
print "$pctGC $ncg $total\n";
