#!/usr/bin/env perl


# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
#  A      B                               C      D       E       F       G       H        I       J      K        L
#pre_22  gi|48994873|gb|U00096.2|        100.00  200     0       0       1       200     4312314 4312513 3e-111   396


open(IN, $ARGV[0]) or die "cannot open $ARGV[0]\n";
$prevQuery = "";
while(<IN>) {
		#       A        B       C       D       E       F       G       H      I        J
#		$_ =~ ;
		$line = $_;
		if ($line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
				$curQuery = $1;  
				$start    = $9;
				$end      = $10;
				if ($curQuery ne $prevQuery) {
						if ($curQuery =~ /^pre/) {
								$preStart = $start;
								$preEnd   = $end;
						}
						if ($curQuery =~ /^suf/) {
								$sufStart = $start;
								$sufEnd   = $end;
								if ($sufStart < $sufEnd) {
										$length = $sufEnd - $preStart;
										print "$preStart $sufEnd $length\n";
								}
								else {
										$length = $preStart - $sufEnd;
										print "$sufEnd $preStart $length\n";
								}
						}
				}
				$prevQuery = $curQuery;
		}
		else {
				print "no match\n";
		}
}
		
	 
