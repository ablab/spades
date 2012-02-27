#!/usr/bin/env perl
if ($#ARGV < 1) {
		print "usage: contigs length\n";
		exit(0);
}
$in = shift @ARGV;
$length = shift @ARGV;

open(IN, "$in") or die "cannot open $in\n";
$contigStart = 0;
$contigEnd   = 0;
$contig = "";
while(<IN>) {
		$line = $_;
		if ($line =~/^>(.*)/) {
				$curTitle = $1;
				$contigStart = 0;
				if (length($contig) > $length) {
						$pre = substr $contig, 0, $length;
						$suf = substr $contig, length($contig) - $length, $length;
						print ">pre_$prevTitle\n";
						print "$pre\n";
						print ">suf_$prevTitle\n";
						print "$suf\n";
				}
				$contig = "";
				$prevTitle = $curTitle;
		}
		else {
				chomp $line;
				$contig .= $line;
		}
}
