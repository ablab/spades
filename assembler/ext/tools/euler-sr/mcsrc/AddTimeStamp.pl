#!/usr/bin/env perl

$file = shift @ARGV;

$filels = `ls -l $file`;
chomp $filels;
@stamp = split(/\s+/, $filels);
$datestr = @stamp[5];

@date = split(/\-/, $datestr);
$datestr = "$date[1]\\\/$date[2]\\\/$date[0]"; 
`perl -pi -e "s/ \\* Last modified:  .*/ \* Last modified:  $datestr/" $file`;

