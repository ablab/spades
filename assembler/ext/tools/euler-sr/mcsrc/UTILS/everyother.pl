#!/usr/bin/env perl
$line = 0;
while(<>) {
 if ($line == 0) {
   print $_;
   $line = 1;
 }
 else {
   $line = 0;
 }
}
   
