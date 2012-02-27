#!/usr/bin/env perl
use strict;

my %rc = ();

$rc{'a'} = 't';
$rc{'t'} = 'a';
$rc{'g'} = 'c';
$rc{'c'} = 'g';
$rc{'A'} = 'T';
$rc{'T'} = 'A';
$rc{'C'} = 'G';
$rc{'G'} = 'C';


while ($#ARGV >= 0) {
  my $str = shift @ARGV;
  
  my $revseq = "";
  my $i;
  for $i (0.. (length($str)-1)) {
    $revseq = $revseq . $rc{chop($str)};
  }
  print "$revseq\n";
  exit
}

my $seq = "";
my $sub;
my $i;
my ($str, $revseq);
while (<>) {
  $str = $_;
  chomp $str;
  $revseq = "";
  if ($str =~/^>/) {
    print "$str.rc\n";
    if ($seq != "") {
      for $i (0.. (length($str)-1)) {
	$revseq = $revseq . chop($seq);
      }
      while (length($revseq) > 0) {
	$sub = substr($revseq, 0, 60);
	print "$sub\n";
	$revseq = substr($revseq, 60);
      }
    }
    $seq = "";
  }
  else {
    $seq = $seq . $str;
  }
}
if ($seq ne "") {
  my $l = length($seq);
  for $i (0.. (length($seq)-1)) {
    $revseq = $revseq . chop($seq);
  }
  while (length($revseq) > 0) {
    $sub = substr($revseq, 0, 60);
    print "$sub\n";
    $revseq = substr($revseq, 60);
  }
}
