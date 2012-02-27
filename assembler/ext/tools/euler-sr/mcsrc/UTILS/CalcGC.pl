#!/usr/bin/env perl
if ($#ARGV < 0) {
  print "usage: $0 fasta_file\n";
  print "output: title<TAB>length<TAB>%GC\n";
  exit(0);
}
$infile = shift @ARGV;

open (IN, "$infile\n");

$numgc = 0;
$length = 0;
$title = "_NONE_";
while(<IN>) {
  $line = $_;
  chomp $line;
  if ($line =~ /^>(.*)/) {
    if ($title ne "_NONE_") {
      if ($length == 0) {
	print "$title\tNAN\tNAN\n";
      }
      else {
	$pctgc = $numgc / $length;
	print "$title\t$length\t$numgc\n";
      }
    }
    $title = $1;
    $numgc = 0;
    $length = 0;
  }
  else {
    $numgc = $numgc + tr/GgCc//;
    $length = $length + length($line);
  }
}

if ($title ne "_NONE_") {
  if ($length == 0) {
    print "$title\tNAN\tNAN\n";
  }
  else {
    $pctgc = $numgc / $length;
    print "$title\t$length\t$numgc\n";
  }
}
