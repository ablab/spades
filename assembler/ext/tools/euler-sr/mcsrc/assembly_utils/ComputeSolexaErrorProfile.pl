#!/usr/bin/env perl

$in = shift @ARGV;
open(IN, $in) or die "cannot open $in\n";

$line = <IN>;
while ($line =~ /^\#/) {
		$line = <IN>;
}
$totReads = 0;
@errorProf = ();
@netErrorProf = ();
@mutProf = ();
for ($i = 0; $i < 100; $i++) {
		$profInit[$i] = 0;
}
do {
		$line =~ /(\w+) (\d+) (\d+) (\d+) ([FR]) (\w+)/;
		$read = $1;
		if ($read  !~ /\./) {
				$seq = $6;
				$dir = $5;
				if ($dir eq "R") {
						$seq =~ tr/ACTG/TGAC/;
						$seq = reverse($seq);
				}
				if (length($read) == length($seq)) {
						$len = length($read);
						$readErrors = 0;
						@errors = ();
						for ($i = 0; $i < $len; $i++) {
								$r = substr($read, $i, 1);
								$s = substr($seq, $i, 1);
								if ($r ne $s) {
										$readErrors++;
								push @errors, $i;
								}
						}
						if ($#errors >= 0  and 
								$profInit[$#errors] == 0) {
								$profInit[$#errors] = 1;
								for ($e = 0; $e < $len; $e++) {
										$errorProf[$#errors][$e] = 0;
								}
						}
						for ($e = 0; $e <= $#errors; $e++) {
								$errorProf[$#errors][$errors[$e]]++;
								$netErrorProf[$errors[$e]]++;
								if ($#errors < 7) {
										if ( $#mutProf < 0) {
												for ($m = 0; $m < $len; $m++) {
														$mutProf[$m] = 0;
												}
										}
										$mutProf[$errors[$e]]++;
								}
						}
						
				}
		}
		$totReads++;
} while ($line = <IN>);

#for ($ne = 0; $ne <= $#errorProf; $ne++ ) {
#		$numErrors = $ne + 1;
#		print "$numErrors @{$errorProf[$ne]}\n";
#}
#print "net: @netErrorProf\n";
print "$totReads @mutProf\n";
