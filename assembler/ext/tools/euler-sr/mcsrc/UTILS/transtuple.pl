#!/usr/bin/env perl
use POSIX;
$tuple = shift @ARGV;
$length = shift @ARGV;

@nucs = ('g', 'a', 'c', 't');
@seq = ();
for ($i =0 ; $i < $length; $i++) {
  $n = $tuple % 4;
  @seq[$length - $i - 1] = $nucs[$n];
  $tuple = POSIX::floor($tuple / 4);
}

print @seq;
print "\n";

