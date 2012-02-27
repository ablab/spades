#!/usr/bin/env perl
#!/usr/bin/env perl


if ($#ARGV < 1) {
		print "usage: RandomlyAssignN.pl in out\n";
		exit(0);
}
$in = shift @ARGV;
$out = shift @ARGV;
@nucs = ('G', 'A', 'C', 'T');
$index = 0;
open(IN, $in) or die "cannot open $in\n";
open(OUT, ">$out") or die "cannot open $out\n";
$trimFront = 0;
while ($#ARGV >= 0) {
		$opt = shift @ARGV;
}

$curTitle = "";
$prevTitle = "";
$last = 0;
while(!$last) {
		$line = <IN>;

		if ($line eq "") {
				$last = 1;
		}

		if ($last || ($line =~/^>/)) {
				$curTitle = $line;

				# process the previous read.
				@readChars = split(//, $readSeq);
				for ($i = 0; $i <= $#readChars; $i++) { 
						if ($readChars[$i] eq 'N') {
								$readChars[$i] = $nucs[$index];
								$index++;
								$index = $index % 4;
						}
				}
				$readSeq = join("", @readChars);
				if ($prevTitle ne "") {
						print OUT "$prevTitle";
						print OUT "$readSeq\n";
				}
				$readSeq = "";
				$prevTitle = $curTitle;
				if ($last) {
						exit 0;
				}
		}
		else {
				chomp $line;
				$readSeq .= $line;
		}
}
