#!/usr/bin/env perl

open(FA, ">fileA.txt");
open(FB, ">>fileA.txt") or die "cannot append to fileA.txt\n";


print FA "line1\n";
print FB "line2\n";
print FA "line3\n";
print FB "line4\n";
