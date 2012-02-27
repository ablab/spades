#!/usr/bin/env perl

if ($#ARGV < 1) {
		print "usage: PrintRetainedReads.pl readFile pathFile\n";
		exit(0);
}

$readsFile = shift @ARGV;
$pathFile  = shift @ARGV;


open(RF, "$readsFile") or die "cannot open $readsFile\n";
open(PF, "$pathFile") or die "cannot open $pathFile\n";


while(<RF>) {
		$readTitle = $_;
		$read      = <RF>;
		$pathFor = <PF>;
		$pathRev = <PF>;
		if ($pathFor !~ /^0/) {
				print $readTitle;
				print $read;
		}
}
