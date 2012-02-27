#!/usr/bin/env perl
if ($#ARGV < 1) {
		print "usage: QFastaToMakedFasta.pl in_sequence.txt qualLimit\n";
		exit(0);
}
$in = shift @ARGV;
$qualLimit = shift @ARGV;

open(IN, $in) or die "cannot open $in\n";

while(<IN>) {
		$title = $_;
		$title = substr($title, 1);
		$read = <IN>;
		chomp $read;
		$orig = $read;
		$qtitle = <IN>;
		$qual = <IN>;
		chomp $qual;
		@qualVals = split(//, $qual);
#		print "qual: $qual, qv: @qualVals\n";
		@readNucs = split(//, $read);
		for ($i = 0; $i <= $#qualVals; $i++) {
				$q = $qualVals[$i];
				$n = ord($q) - 64;
#				print "$q $n\n";
				if ($n < $qualLimit) {
						$readNucs[$i] = lc($readNucs[$i]);
				}
		}
		$read = join("", @readNucs);
		print ">$title";
		print "$read\n";
}



