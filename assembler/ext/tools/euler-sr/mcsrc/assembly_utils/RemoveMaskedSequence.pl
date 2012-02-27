#!/usr/bin/env perl
#if ($#ARGV < 1) {
#		print "usage: RemoveMaskedSequence.pl [-maxMasked m (3)]\n";
#		exit(0);
#}

#$in = shift @ARGV;
#$out = shift @ARGV;
$maxMasked = 3;
while ($#ARGV >= 0 ) {
		$opt = shift @ARGV;
		if ($opt eq "-maxMasked") {
				$maxMasked = shift @ARGV;
		}
}
#open(IN, "$in") or die "cannot open $in\n";
#open(OUT, ">$out") or die "cannot open $out\n";
$seq = "";
$title = "";
$lineLength = 60;
$readNumber = 0;
while(<>) {
		if ($_ =~ /^>/) {
				$newTitle = $_;
				$readNumber++;
				if ($title ne "") {
						# process the sequence for X's.
						# first count the vector seq
						# First remove all N's from the beginning and end.
						$seq =~ /[N]*([^N].*[^N])[N]*/;
						$seq = $1;
						$seqCopy = $seq;
						$nN = $seqCopy =~ tr/N/n/;
						$seqLen = length($seq);
						if (($seqLen == 0) || ($nN / $seqLen > 0.5)) {
								$discard = 1;
						}
						else {
								# find the longest non-vector sequence, and print that.
								if ($nN > 3 ) {
										@seqA = ();
										@seqA = split(//, $seq);
										if ($#seqA > 1) {
												# charge the first boundary
												$s = 0;
												@starts = ();
												@ends   = ();
												@lengths = ();

												while ($s <= $#seqA) {
														$sprev = $s;
														while ($s <= $#seqA and
																	 $seqA[$s] ne 'A' and
																	 $seqA[$s] ne 'T' and
																	 $seqA[$s] ne 'G' and
																	 $seqA[$s] ne 'C') {
																$s++;
														}
														$start = $s;
														while ($s <= $#seqA and
																	 ($seqA[$s] eq 'A' or
																		$seqA[$s] eq 'T' or
																		$seqA[$s] eq 'G' or
																		$seqA[$s] eq 'C')) {
																$s++;
														}
														$end = $s;
														push @starts, $start;
														push @ends, $end;
														push @lengths, $end - $start;
														if ($sprev == $s) {
																print "something bad with the loop\n";
																exit(0);
														}
												}
										}
										$cur = 0;
										$longest = 0;
										for ($cur = 0; $cur <= $#lengths; $cur++ ) {
												if ($lengths[$cur] > $lengths[$longest]) {
														$longest = $cur;
												}
										}
										chomp $title;
										$trimEnd = length($seq) - $ends[$longest];
										$title .= " trimFront=@starts[$longest] trimEnd=$trimEnd\n";
										$truncated = substr $seq, $starts[$longest], $ends[$longest] - $starts[$longest];
								}
								else {
										$truncated = $seq;
								}
								print  $title;
								while (length($truncated) > 0) {
										$line = substr $truncated, 0, $lineLength;
										print "$line\n";
										$truncated = substr $truncated, $lineLength;
								}
						}
				}
				$title = $newTitle;
				$seq = "";
		}
		else {
				chomp;
				$seq .= $_;
		}
}
