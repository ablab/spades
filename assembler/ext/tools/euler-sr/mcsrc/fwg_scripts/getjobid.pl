#!/usr/bin/env perl

$base = shift @ARGV;
@resList = `qstat -u $ENV{"USER"} | grep $base`;
$nres = $#resList + 1;
if ($nres > 0) {
	foreach $res (@resList) {
		$res =~ /(\d+) .*/;
		$pids .= "$1 ";
	}
	print "$pids\n";
}
