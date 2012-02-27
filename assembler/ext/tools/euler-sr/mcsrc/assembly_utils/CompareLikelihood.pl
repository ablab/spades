#!/usr/bin/env perl
use POSIX;

if ($#ARGV < 2) {
		print "usage: CompareLikelihood.pl error_profile.txt nErrors file.aln [file2.aln...] [options]\n";
		print "   error_profile is the per-position error profile of reads\n";
		print "   nerrors is the number of errors to consider when checkign for\n";
		print "   catastrophic errors\n";
		print "   file.aln ... are the alignment files\n";
		print " -printIndependent 'file' print independent error profiles to 'file'\n";
		print " -pcf 'file' print catastrophies to fasta 'file'\n";
		exit(1);
}

$epFile = shift @ARGV;
$nErrCheck = shift @ARGV;

open(EP, "$epFile") or die "cannot open $epFile\n";
$epLine = <EP>;
@ep = split(",", $epLine);

$errRegion = POSIX::floor($nErrCheck * 4 / 3);
@nIndep = ();
@nCat   = ();
for ($p = 0; $p <= $#ep; $p++) {
		$nIndep[$p] = 0;
		$nCat[$p] = 0;
}

$printIndependent = 0;
while ($#ARGV >= 0) {
		$opt = shift @ARGV;
		if ($opt eq "-printIndependent") {
				$printIndependent = 1;
				$split  = shift @ARGV;
		}
		elsif ($opt eq "-pcf") {
				$catasFastaFile = shift @ARGV;
		}
		else {
				push @alignFiles, $opt;
		}
}

$readIndex = 0;
if ($printIndependent) {
		$indepFile = "$split.indep";
		$cataFile = "$split.catas";
		open(IND, ">$indepFile") or die "cannot open $indepFile\n";
		open(CAT, ">$cataFile") or die "cannot open $cataFile\n";
}

print "cata fasta file: $catasFastaFile\n";

print "nalign: $#alignFiles\n";
$fastaOutOpened = 0;
while($#alignFiles >= 0) {
		$alnFile = shift @alignFiles;
		open(ALN, "$alnFile") or die "cannot open $alnFile\n";
		while(<ALN>) {
				$alignment = $_;
				$printSeq = 0;
				if ($alignment =~ /(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\S+)/) {
						$start = $1; $strand = $2;
						$alignStr = $3; $nErr = $4;
						$seqStr = $5;
						$printSeq = 1;
						if ($fastaFileOpened == 0) {
								$fastaFileOpened = 1;
#								open(INDF, ">$indepFastaFile") or die "cannot open $indepFastaFile\n";
								open(CATF, ">$catasFastaFile")  or die "cannot open $cataFastaFile\n";
						}
#						print "$start $strand $alignStr $nErr\n";
				}
				elsif ($alignment =~ /(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/) {
						$start = $1; $strand = $2;
						$alignStr = $3; $nErr = $4;
				}						

				@align = split(//, $alignStr);
				$lenAlign = scalar @align;
				if ($nErr >= $nErrCheck) {
						$logIndExp = 0;
						$logCatExp = 0;
						for ($i = 0; $i < $nErr; $i++) {
								$pos = $i + $lenAlign - $nErr;
								$lenEp = scalar @ep;
								$errVal = @align[$i + ($lenAlign - $nErr)];
								$logIndExp += (1-$errVal)*log(1-$ep[$pos]) + ($errVal * log($ep[$pos]));
								$logCatExp += (1-$errVal)*log(0.25) + ($errVal)*log(0.75);
						}

						if ($logIndExp  > $logCatExp) {
								$nIndep[$nErr]++;
								print IND $alignment;
						}
						elsif ($logCatExp > $logIndExp ) {
								$nCat[$nErr]++;
								print CAT $alignment;
								if ($printSeq) {
										print CATF ">$readIndex\n";
										print CATF "$seqStr\n";
								}
						}

#						print "logInd: $logIndExp logCat: $logCatExp align: $alignStr\n";
				}
				else {
						print IND $alignment;
				}
				$readIndex++;
		}
		close ALN;
}

for ($i = 0; $i <= $#ep; $i++) {
		print "$i $nIndep[$i] $nCat[$i]\n";
}



