#!/usr/bin/env perl

if ($#ARGV < 1) {
  print "usage: extractSequences seq_title_list fasta_file\n\n";
  print "the titles in 'seq_title_list must be in the order of their\n";
  print "appearance in 'fasta_file'\n";
  print "They will be printed to stdout\n";
  exit(0);
}

$tf = shift @ARGV;
$sf = shift @ARGV;
open(TITLES, "$tf") or die "cannot open $tf\n";
open(SEQS, "$sf") or die"cannot open $sf\n";

@titles = <TITLES>;
chomp @titles;
$titleMap = {};
for ($t = 0; $t <= $#titles; $t++) {
#		print "storing title: $titles[$t]\n";
		$titleMap{$titles[$t]} = 1;
}

chomp @titles;
$curTitle = 0;
$printSeq = 1;
while(<SEQS>) {
  $line = $_;
#  print "looking for $titles[$curTitle] $printSeq $curTitle\n";
  $ct = $titles[$curTitle];
  if ($line =~ />(\S+)/) {
    $title = $1;
#		print "using title: $title\n";
  }
  if (exists $titleMap{$title}) {
    $lin = $line;
    chomp $lin; 
    print "$line";
    $printSeq = 1;
    $curTitle++;
  }
  elsif ($line !~ />/ && $printSeq == 1) {
    print "$line";
  }
  elsif ($line =~ />/) {
    $printSeq = 0;
  }
}
