#!/usr/bin/env perl
use POSIX;

$st = "";
push  @nucs, 'g', 'a', 'c', 't';

$len = shift @ARGV;

while ($#ARGV >= 0) {
  $templen = $len;
  $num = shift @ARGV;
  $st = "";
  while ($templen > 0) {
    $rem = $num % 4;
    $st = $st . @nucs[$rem] ;
    $num = POSIX::floor($num / 4);
    $templen--;
  }
  
  print "$st\n";
}

