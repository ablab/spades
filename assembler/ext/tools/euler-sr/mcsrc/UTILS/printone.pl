#!/usr/bin/env perl

while (<>) {
  / (\d+)/ && print "$1\n";
}
