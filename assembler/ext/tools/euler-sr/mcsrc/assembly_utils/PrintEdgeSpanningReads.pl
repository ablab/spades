#!/usr/bin/env perl

$readsFile = shift @ARGV;
$pathFile  = shift @ARGV;
$outFile   = shift @ARGV;

open(RF, $readsFile) or die "cannot open $readsFile\n";
open(PF, $pathFile) or die "cannot open $pathFile\n";
open(OF, ">$outFile") or die "cannot open $outFile\n";


$printRead = 0;
# discards the # reads line
<PF>;
$pathLine = <PF>;
$pathLine =~ /(\d+)/;
$numEdges = $1;
# discard the RC line.
<PF>;

$read = "";
$title = "";
$curTitle = "";
while(<RF>) {
		
		$readLine = $_;
		if ($readLine =~ /^>/) {
				$curTitle = $readLine;
				if ($title ne "") {
						if ($numEdges > 1) {
								print OF $title;
								print OF "$read\n";
						}
						$pathLine = <PF>;
						$pathLine =~ /(\d+)/;
						$numEdges = $1;
# discard the RC line.
						<PF>;
				}
				$title = $curTitle;
				$read  = "";
		}
		else {
				chomp $readLine;
				$read .= $readLine;
		}
}
