#!/usr/bin/env perl
$pttFile = shift @ARGV;
$coordsFile = shift @ARGV;

open(PTT, $pttFile) or die "cannot open $pttFile\n";
open(CF, $coordsFile) or die "cannot open $coordsFile\n";

@origCoords = ();
$cn = 0;
while(<CF>) {
		$cl = $_;
		$cl =~ /(\d+)\s+(\d+).*/;
		$start = $1; $end = $2;
		push @origCoords, [($start, $end)];
#		print "$coords[$cn], @{$coords[$cn]}\n";
		$cn++;
}

@coords = sort {$a->[0] <=> $b->[0]} @origCoords;

# first 2 lines of ptt file are documentation
<PTT>;
<PTT>;

$numMatch = 0;
$numNotMatch = 0; 
while(<PTT>) {
		$_ =~ /(\d+)\.\.(\d+).*/;
		$start = $1; $end = $2;
		$matched = 0;
#		print "checking $start $end\n";
		for ($i = 0; $i <= $#coords && $matched == 0; $i++) {
				if ($coords[$i][0] <= $start &&
						$coords[$i][1] >= $end) {
							$numMatch++;
						$matched = 1;
				}
		}
		if ($i == $#coords+1) {
				print $_;
				$numNotMatch++;
		}
}
print "$numMatch $numNotMatch\n";
