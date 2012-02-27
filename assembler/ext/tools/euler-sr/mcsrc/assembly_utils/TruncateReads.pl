#!/usr/bin/env perl
if ($#ARGV != 1) {
 print "usage: TruncateReads.pl readsFile readLength\n";
 exit(1);
} 
$rf = shift @ARGV;
$readLen = shift @ARGV;
open(RF, "$rf") or die "cannot open $rf\n";
$title = "";
while(<RF>) {
		if ($_ =~ /^>/) {
				if ($title ne "") {
						$subseq = substr($seq, 0, $readLen);
						print $title;
						print "$subseq\n";
				}
				$title = $_;
				$seq = "";
		}
		else {
				chomp($_);
				$seq .= $_;
		}
}
