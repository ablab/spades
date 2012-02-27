#!/usr/bin/env perl

use findReadInfo;


if ($#ARGV != 0) {
  print "usar: $0 reads_file.fasta\n";
  exit(0);
}


$in = shift @ARGV;

open (FASTAIN, $in) or die "no puedo abrir $in\n";

while (<FASTAIN>) {
  $line = $_;
  if ($line =~ /^>(.*)/) {
    $title = $1;
    ($organismo, $largoClon, $read, $clon, $ext) = &findReadInfo($title);
    print "$title\torg: $organismo\tlargo: $largoClon\tread: $read\tclon: $clon\text: $ext\n";
  }
}
