#!/usr/bin/env perl

$s = "hello joe";
@ss = split(//, $s);
$js = join(":", @ss);
print "js: $js\n";

