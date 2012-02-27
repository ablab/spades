#!/usr/bin/env perl
use strict;

if ($#ARGV != 1) {
  print "usage: $0 refrence query\n";
  exit(-1);
}


my $k = 5;

my $ref = shift @ARGV;
my $query = shift @ARGV;

open (REF, "$ref") or die "cannot open $ref\n";

my $refseq = "";
my $line;

my $label = "";
my %references = ();
$line = <REF>;
chomp($line);
$line =~ /\W+(\w+)/;
$label = $1;
while (<REF>) {
  $line= $_;
  chomp($line);
  if ($line =~ />/) {
    $references{$label} = $refseq;
    $line =~ />(\w+) .*/;
    $label = $1;
    $refseq = "";
  }
  else {
    $refseq = $refseq . $line;
  }
}
$references{$label} = $refseq;
close(REF);

open (QUERY, "$query") or die "cannot open $query\n";

my $queryseq = "";
my %queries = ();
my $label = "";
$label = <QUERY>;
chomp($label);
$label =~ /\D+(\d+) .*/;
$label = int($1);

while (<QUERY>) {
  $line= $_;
  chomp($line);
  if ($line =~ />/) {
    $queries{$label} = $queryseq;
    $line =~ /\D+(\d+)/;
    $label = int($1);
    $queryseq = "";
  }
  else {
    $queryseq = $queryseq . $line;
  }
}
$queries{$label} = $queryseq;
close(QUERY);


my $word;
my $refpos = 0;

my $index;

#print "refseq: $refseq\n";
#print "queryseq: $queryseq\n";
#printf("length ref: %d  query: %d\n", length($refseq), length($queryseq));

my ($key, $refkey);
my $revcomp;
for $key (sort { $a <=> $b } keys(%queries)) {
  for $refkey (keys (%references)) {
    $index = index(lc($references{$refkey}), lc($queries{$key}));
#    if ($index < 0) { 
#      print "$refkey doesn't have $key\n";
#    }
    while ($index >= 0) {
      printf("$refkey (%d) $key (%d)  $index\n", length($references{$refkey}), length($queries{$key}));
      $index = index(lc($references{$refkey}), lc($queries{$key}), $index+1); 
    }
    $revcomp = rc(lc($queries{$key}));
    $index = index(lc($references{$refkey}), $revcomp);
    while ($index >= 0) {
      printf("$refkey (%d) $key (%d) $index  *\n", length($references{$refkey}), length($queries{$key}));
      $index = index(lc($references{$refkey}), $revcomp, $index+1); 
    }
  }
}      
    
sub rc {

my %rcs = ();
$rcs{'a'} = 't';
$rcs{'t'} = 'a';
$rcs{'g'} = 'c';
$rcs{'c'} = 'g';
my $str = @_[0];
my $revseq = "";
my $i;
for $i (0.. (length($str)-1)) {
  $revseq = $revseq . $rcs{chop($str)};
}
return $revseq;
}
