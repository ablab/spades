#!/usr/bin/env perl
if ($#ARGV < 0) {
		print "usage: RetainHighestBlastHit.pl blasttab\n";
		exit(1);
}

$blastInFile = shift @ARGV;

open(BI, $blastInFile) or die "cannot open $blastInFile\n";

$blastLine = <BI>;
$blastLine =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/;
$prevTitle = $1;
$prevScore = $12;
$prevLine  = $blastLine;

$bestLine  = $prevLine;
$bestScore = $prevScore;

$done = 0;


while($done == 0) {
		# Parse the current score
		# Loop invariant. At the beginning of the loop, $prevTitle is not equal to 
		# prev prev title.
		do {
				if (($blastLine = <BI>)) {
						$blastLine =~/(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/;
						$curTitle = $1;
						$curScore = $12;
						if ($curTitle eq $prevTitle && 
								$curScore > $bestScore) {
								$bestLine = $blastLine;
								$bestScore = $curScore;
						}
				}
				else {
						$done = 1;
				}
		} while ($done == 0 && $curTitle eq $prevTitle);
		print $bestLine;
		$prevTitle = $curTitle;
		$prevScore = $curScore;
		$bestLine  = $blastLine;
}
# print the last line
if ($prevTitle ne $curTitle) {
		print $bestLine;
}


