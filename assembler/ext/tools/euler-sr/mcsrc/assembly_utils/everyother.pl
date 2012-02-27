#!/usr/bin/env perl
$even = 0;
while(<>) {
  if ($even == 0) {
    print $_;
    $even = 1;
  }
  else {
    $even = 0;
  }
}
