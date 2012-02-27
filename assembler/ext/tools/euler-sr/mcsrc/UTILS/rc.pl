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
$rc{'N'} = 'N';
$rc{'n'} = 'n';

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
my $lines = 0;
my $oldTitle = "";
my $newTitle = "";
while (<>) {
  $str = $_;
  chomp $str;
  if ($str =~/^>/ ) {
    $newTitle = "$str.rc";
    if ($seq ne "") {
      print "$oldTitle\n";
      $revseq = "";
      my $l = length($seq) - 1;
      for $i (0.. $l) {
	$revseq = $revseq . $rc{chop($seq)};
      }
      while (length($revseq) > 0) {
	$sub = substr($revseq, 0, 60);
	print "$sub\n";
	$revseq = substr($revseq, 60);
      }
    }
    $oldTitle = $newTitle;
    $seq = "";
  }
  else {
    $seq = $seq . $str;
    $lines = $lines + 1;
  }
}
# print the sequence
if ($seq ne "") {
  print "$newTitle\n";
  my $l = length($seq);
  $revseq = "";
  for $i (0.. (length($seq)-1)) {
    $revseq = $revseq . $rc{chop($seq)};
  }
  while (length($revseq) > 0) {
    $sub = substr($revseq, 0, 60);
    print "$sub\n";
    $revseq = substr($revseq, 60);
  }
}
