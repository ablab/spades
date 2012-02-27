#!/usr/bin/env perl

$in = shift @ARGV;
open(IN, $in) or die "cannot open $in\n";

while(<IN>) {
 $title = $_;
 $seq = <IN>;
 $title2 = <IN>;
 $qual = <IN>;
 $title =~ /^@(.*)/;
 print ">$1\n";
 print $seq;
}
