#!/usr/bin/env perl
if ($#ARGV < 1) {
		print "usage: MaskTupleList.pl list mask\n";
		exit(0);
}

$list = shift @ARGV;
$mask = shift @ARGV;

open(LI, $list) or die "cannot open $list\n";
open(MA, $mask) or die "cannot open $mask\n";

<LI>;
<MA>;
if (!<LI>) {
		exit 1;
}
$l = $_;

if (!<MA>) {
		exit 1;
}
$m = $_;

while(1) {
		$l =~ /(\S+) (\d+)/;
		$li = $1;
		$mult = $2;
		$m =~ /(\S+)/;
		$mi = $1;
		if ($li eq $mi) {
				print "$mult\n";
				if (!($l = <LI>)) { last; }
				if (!($m = <MA>)) { last; }
		}
		elsif ($li lt $mi) {
				if (!($l = <LI>)) { last; }
		}
		else {
				if (!($m = <MA>)) { last; }
		}
}				
		

		
	  


