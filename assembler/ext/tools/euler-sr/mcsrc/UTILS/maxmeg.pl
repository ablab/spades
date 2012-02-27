#!/usr/bin/env perl
use strict;

if ($#ARGV ne 1) {
  printf "usage: maxmeg infile outfile\n";
  exit(-1);
}

my $infile = shift @ARGV;
my $outfile = shift @ARGV;

open (IN, "$infile") or die "cannot open $infile\n";

my %alignscore = ();
my %alignstring = ();
<IN>;<IN>;<IN>;

while (<IN>) {
  if (/WARNING/) { next;}
  /\w*\.con(\d+)\t\w+\t+(\d+)/;
  if (exists $alignscore{$1}) {
    if ($alignscore{$1} < $2) {
      $alignscore{$1} = int($2);
      $alignstring{$1} = $_;
    }
  }
  else {
    $alignscore{$1} = $2;
    $alignstring{$1} = $_;
  }
}
close(IN);
my $i;

sub sortint {
  my($a, $b) = @_;
  return int($a) < int($b);
}

open(OUT, ">$outfile") or die "cannot open $outfile\n";

for $i (sort{ $a <=> $b }  keys(%alignstring)) {
  print OUT $alignstring{$i};
}

close(OUT);


