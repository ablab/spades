#!/usr/bin/env perl
if ($#ARGV < 2) {
  print "usage: FilterReads.pl reads.fasta nameParseFile outputName\n";
	print "     -maxC pct - filter reads if they are pct C or more\n";
	print "     -maxA pct - same for A\n";
  exit(0);
}

$readsName = shift @ARGV;
$mateTableName = shift @ARGV;
$outputName = shift @ARGV;

$maxA = 0.95;
$maxC = 0.95;

while ($#ARGV >= 0) {
		$opt = shift @ARGV;
		if ($opt eq "-maxA") {
				$maxA = shift @ARGV;
		}
		elsif ($opt eq "-maxC") {
				$maxC = shift @ARGV;
		}
		else {
				print "bad option: $opt\n";
				exit(0);
		}
}

open(RI, $readsName) or die "cannot open $readsName\n";
open(RO, ">$outputName") or die "cannot write to $outputName\n";
open(MT, $mateTableName) or die "cannot open $mateTableName\n";
@regexps = ();
@types = ();
while(<MT>) {
		$_ =~ /"([^"]+)".*Type=(\d+)/; 
		$reg = $1;
    $type= $2;
#    print "$type, got regexp: '$reg' from $_";
		$reg =~ s/\\\(/\(/g;
		$reg =~ s/\\\)/\)/g;
		push @regexps, $reg;
		
    push @types, $type;
}

$title = <RI>;
$type  = -1;
for ($r = 0; $r <= $#regexps; $r++) {
		if ($title =~ /$regexps[$r]/) {
				$types= $types[$r];
				last;
#				print "type: $type   of $title using $regexps[$r]\n";
		}
}
$read = "";
while(<RI>) {
		$line = $_;
		if ($line =~ /^>/) {
				# proces the current read.
				if ($type == 0) {
						# 0 is 454 regexp

				} 
				elsif ($type == 1) {
						# 1 is Illumina regexp
						$numA = ($read =~ tr/A/A/ );
						$numC = ($read =~ tr/C/C/ );
						$numG = ($read =~ tr/G/G/ );
						$numT = ($read =~ tr/T/T/ );
						
						$length = length($read);
						$readIsGood = 1;
						if ($length == 0) {
								$readIsGood = 0;
						}
						elsif ($numC / $length > $maxC) {
								$readIsGood = 0;
						}
						elsif ($numA / $length > $maxA) {
								$readIsGood = 0;
						}
						elsif ($numC / $length > 0.85 && $numA / $length > 0.15) {
								$readIsGood = 0;
						}
						elsif ($read =~ /AAAAAAAAAAAAAAAAAAAAAAAAA/) {
								$readIsGood = 0;
						}
						elsif ($read =~ /CCCCCCCCCCCCCCCCCCCCCCCCCCC/) {
								$readIsGood = 0;
						}
						elsif ($read =~ /GGGGGGGGGGGGGGGGGGGGGGGGG$/) {
								$readIsGood = 0;
						}
						if ($readIsGood) {
								print RO $title;
								print RO "$read\n";
						}
						else {
#								print "filtering\n";
#								print "$title";
#								print "$read\n";
						}
				}


				# if the next read has been read, bail
				$title = $line;
				$type = -1;
				for ($r = 0; $r <= $#regexps; $r++) {
						if ($title =~ /$regexps[$r]/) {
								$type = $types[$r];
								last;
#								print "$title using $regexps[$r]\n";
						}
				}
				$read = "";
		}
		else {
				chomp $line;
				$read .= $line;
		}
}
