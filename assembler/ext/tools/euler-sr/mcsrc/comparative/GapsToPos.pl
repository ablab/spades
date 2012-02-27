#!/usr/bin/env perl 

open(IN, @ARGV[0]) or die "cannot open @ARGV[0]\n";
while (<IN>) {
#  /(\w+) (\w+) ([\w\.]+) (\w+) (\d+) (\d+) (\d+) (\d+) (\d+) (\w)/ && (print("$6 $3\n"));
  /(\w+) (\w+) (\w+)\.fa (\w+) (\d+) (\d+)/ && (print("$6 $3\n"));
}

