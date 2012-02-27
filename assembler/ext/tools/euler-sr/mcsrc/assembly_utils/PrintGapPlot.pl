#!/usr/bin/env perl

$file  = shift @ARGV;
$minLength = shift @ARGV;
$start = shift @ARGV;
$end   = shift @ARGV;

for ($kmer = $start; $kmer <= $end; $kmer++) {
		$n50 = `~/projects/mcsrc/assembly_utils/FindGapsInCoverage.pl $file $kmer $minLength | N50`;
		print "$kmer $n50";
}
