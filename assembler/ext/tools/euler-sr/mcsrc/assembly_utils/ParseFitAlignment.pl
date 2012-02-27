#!/usr/bin/env perl

if ($#ARGV < 0) {
		print "usage: ParseFitAlignment.pl file.aln [file2.aln ...]\n";
		exit(1);
}
$readLength = -1;

while ($#ARGV >= 0) {
		push @alignFiles, shift @ARGV;
}

@errorProfile = ();
@errorCount = ();


$nReads = 0;

foreach $alignFile (@alignFiles) {
	open(FIT, "$alignFile") or die "cannot open $alignFile\n";
	
	# parse alignment line
	while(<FIT>) {
		$alignment= $_;
		if ($alignment =~ /(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)/) {
			$start = $1; $strand = $2;
			$alignStr = $3; $nErr = $4;
			$seqStr = $5;
			$printSeq = 1;
		}
		elsif ($alignment =~ /(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/) {
			$start = $1; $strand = $2;
			$alignStr = $3; $nErr = $4;
		}
		
		
		
		if ($readLength == -1) {
			# init everything
			$readLength = length $alignStr;
			for ($i = 0; $i < $readLength; $i++) {
				$errorProfile[$i] = ();
				for ($j = 0; $j < $readLength; $j++) {
					push @{$errorProfile[$i]}, 0;
				}
				push @errorCount, 0;
			}
		}

		++$nReads;
		@align = split(//, $alignStr);
		for ($i = 0; $i < $readLength; $i++) {
			$errorProfile[$nErr][$i] += @align[$i];
		}
		@errorCount[$nErr]++;
	}
}

@posTotal = ();
for ($p = 0; $p < $readLength; $p++) {
		push @posTotal, $0;
}

#for ($p = 0; $p < $readLength; $p++) {
#		print "$p,";
#}
for ($nErr = 0; $nErr < $readLength; $nErr++) {
#		print "$nErr,";
		for ($p = 0; $p < $readLength; $p++) {
#				print "$errorProfile[$nErr][$p],";
				$posTotal[$p] += $errorProfile[$nErr][$p];
		}
#		print "\n";
}
#
#for ($p = 0; $p < $readLength; $p++) {
#		print "$posTotal[$p],";
#}
#print "\n";
for ($p = 0; $p < $readLength; $p++) {
		$frac = $posTotal[$p] / $nReads;
		print "$frac ";
}
print "\n";
#for ($p = 0; $p < $readLength; $p++) {
#		print "$errorCount[$p],";
#}
#print "\n";
#for ($p = 0; $p < $readLength; $p++) {
#		$frac = $errorCount[$p] / $nReads;
#		print "$frac,";
#}
#
#print "\n";

#print "$nReads\n";
