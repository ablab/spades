#!/usr/bin/env perl

#$in = shift @ARGV;
#$out = shift @ARGV;
#open(IN, "$in") or die "cannot open $in\n";
#open(OUT, ">$out") or die "cannot open $out\n";
$seq = "";
$title = "";
$lineLength = 60;
while(<>) {
		if ($_ =~ /^>/) {
				$newTitle = $_;
				if ($title ne "") {
						# process the sequence for X's.
						# first count the vector seq
						$seqCopy = $seq;
						$nX = $seqCopy =~ tr/X/x/;
						$seqLen = length($seq);
						if ($nX / $seqLen > 0.5) {
								$discard = 1;
						}
						else {
								# find the longest non-vector sequence, and print that.
								if ($nX > 0 ) {
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
