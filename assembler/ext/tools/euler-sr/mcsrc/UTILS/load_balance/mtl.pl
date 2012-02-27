#!/usr/bin/env perl


open (IN, shift @ARGV);
$machines = "";
while (<IN>) {
  $line = $_;
  $line =~ /([\S]+)/;
  $machines = $machines . $1 . " ";
}
print $machines;
