#!/usr/bin/env perl


$spectFile= shift @ARGV;
open(SP, "$spectFile") or die "cannot open $spectFile\n";
<SP>; # discard N reads;
$n1 = 0; $n2 = $0; $n3 = 0;
while(<SP>) {
		$_ =~ /(\S+)\s+(\d+)/;
		$mult = $2;
		if ($mult == 1) { $n1++;}
		elsif ($mult == 2) {$n2++;}
		else {$n3++;}
}
print "$n1 $n2 $n3\n";
