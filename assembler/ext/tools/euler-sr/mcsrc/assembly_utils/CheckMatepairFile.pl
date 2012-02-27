#!/usr/bin/env perl

use POSIX;

$matepairName = shift @ARGV;
$edge1ReadName = shift @ARGV;
$edge2ReadName = shift @ARGV;


open(MP, $matepairName) or die "cannot open $matepairName\n";

open(E1, $edge1ReadName) or die "cannot open $edge1ReadName\n";
open(E2, $edge2ReadName) or die "cannot open $edge2ReadName\n";


@e1paths = <E1>;
@e2paths = <E2>;

# convert paths into read indices
chomp @e1paths;
chomp @e2paths;

@e1read = ();
@e2read = ();
%e1readH = {};
%e2readH = {};

for ($e = 0; $e <= $#e1paths; $e++ ){
		$r = POSIX::floor($e1paths[$e]/2);
		push @e1read, $r;
		$e1readH[$r] = 1;
}

for ($e = 0; $e <= $#e2paths; $e++ ){
		$r = POSIX::floor($e2paths[$e]/2);
		push @e2read, $r;
		$e2readH[$r] = 1;
}
$read = 0;
while(<MP>) {
		$line = $_;
		$line =~ /(\d+) (\d+)/;
		$mateIndex = $1;
		if ((exists $e1readH[$read]) && (exists $e2readH[$mateIndex])) {
				print "$read $mateIndex\n";
		}
		$read++;
}

