#!/usr/bin/env perl

open STDIN, '/dev/null'   or die "Can't read /dev/null: $!";
open STDERR, '>/dev/null' or die "Can't write to /dev/null: $!";
$top = `top -n 1`;

print "got top: $top\n";
