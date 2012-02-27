#!/usr/bin/env perl

if ($#ARGV < 1) {
		print "usage: WatchMem.pl pid outFile\n";
		exit(1);
}

$pid = shift @ARGV;
$outputFile = shift @ARGV;
open(OF, ">$outputFile") or die "cannot open output $outputFile\n";

$sampleTime = 10;
if ($#ARGV >= 0) {
		$sampleTime = shift @ARGV;
}
$statFile = "/proc/$pid/status";
$vmPeak = 0;
$vmUnits = "N/A";
while (-e $statFile ) {
		open(F, $statFile) or die "cannot open $statFile\n";
		@stats = <F>;
		grep(/VmPeak:\s+(\d+) (\w+)/ && ($vmPeak = $1; $vmUnits = $2), @stats);

		sleep($sampleTime);
}

print OF "$vmPeak $vmUnits\n";
