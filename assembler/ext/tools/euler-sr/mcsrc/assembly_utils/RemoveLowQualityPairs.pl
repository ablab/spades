#!/usr/bin/env perl

if ($#ARGV < 2) {
		print "usage: RemoveLoqQualityPairs.pl reads readpairs maxerr\n";
		exit(1);
}

$file1 = shift @ARGV;
$file2 = shift @ARGV;
$maxErr = shift @ARGV;

open(F1, $file1) or die "cannot open $file1\n";
open(F2, $file2) or die "cannot open $file2\n";

open(O1, ">$file1.filt") or die "cannot open $file1.filt\n";
open(O2, ">$file2.filt") or die "cannot open $file2.filt\n";

while(1) {
		$title1 = <F1>;
		$read1 = <F1>;
		$title2 = <F2>;
		$read2 = <F2>;
		if ($title1 eq "") {
				# done
				exit(0);
		}

		$dis = 0;
		if ($read1 =~ /[nN]/|| $read2 =~ /[nN]/) {
				$dis = 1;
		}
		else {
				$nErr1 = ($read1 =~ tr/actg/ACTG/);
				$nErr2 = ($read2 =~ tr/actg/ACTG/);
				if ($nErr1 <= $maxErr && $nErr2 <= $maxErr) {
						print O1 "$title1";
						print O1 "$read1";
						print O2 "$title2";
						print O2 "$read2";
				}
				else {
						$dis = 1;
				}
		}
		if ($dis) {
				$numDis++;
		}
}

print "discarded: $numDis\n";
