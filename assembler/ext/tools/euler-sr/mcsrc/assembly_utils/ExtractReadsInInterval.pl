#!/usr/bin/env perl
if ($#ARGV < 2) {
				print "usage: ExtractReadsInInterval.pl readsFile intvStart intvEnd\n";
				print "         The reads must have the mapstart and mapend keywords in the fasta title.\n";
				exit(0);
}
$readsFile = shift @ARGV;
$intvStart = shift @ARGV;
$intvEnd = shift @ARGV;

open(RF, $readsFile) or die "cannot open $readsFile\n";
$printRead = 0;
while(<RF>) {
		$line = $_;
		
		if ($line =~ /^>.*mapstart=(\d+).*mapend=(\d+)/) {
				$start = $1;
				$end = $2;
				if ($start > $invStart and $end < $intvEnd) {
						print $line;
						$printRead = 1;
				}
				else {
						$printRead = 0;
				}
		}
		elsif ($line =~ /^>.*pos=(\d+)/) {
				$start = $1;
				$end = $start+1;
				if ($start > $invStart and $end < $intvEnd) {
						print $line;
						$printRead = 1;
				}
				else {
						$printRead = 0;
				}
		}
		else {
				if ($printRead) {
						print $line;
				}
		}
}
