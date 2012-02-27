#!/usr/bin/env perl

$readsFile = shift @ARGV;
$pathFile  = shift @ARGV;

open(RF, "$readsFile") or die "cannot open $readsFile\n";
open(PF, "$pathFile") or die "cannot open $pathFile\n";

<PF>;
$printSeq = 0;
$readLine = <RF>;
while(<PF>) {
		$path = $_;
		<PF>;
		chomp $path;
		if ($path =~ /^0/) {
				$printSeq = 1;
		}
		else {
				$printSeq = 0;
		}
		if ($printSeq) {
				print $readLine;
		}
		$cont = 1;
		while (<RF>) {
				$readLine = $_;
				if (/>/) {
						# break from this loop, done printing read
						last;
				}
				if ($readLine =~ />/) {
						last;
				}
				elsif ($printSeq) {
						print $readLine;
				}
		}
}
				
				
				
		
