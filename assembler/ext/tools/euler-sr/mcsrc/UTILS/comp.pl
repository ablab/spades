#!/usr/bin/env perl

$rc{'a'} = 't';
$rc{'t'} = 'a';
$rc{'g'} = 'c';
$rc{'c'} = 'g';
$rc{'A'} = 'T';
$rc{'T'} = 'A';
$rc{'C'} = 'G';
$rc{'G'} = 'C';

$str = shift @ARGV;
for $i (0.. (length($str)-1)) {
  $comp = $rc{chop($str)} . $comp;
}

print "$comp\n";
