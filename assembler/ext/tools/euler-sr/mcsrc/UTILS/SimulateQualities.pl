#!/usr/bin/env perl
use strict;
if ($#ARGV != 1) {
  print "usage: $0 readsfile qualitiesfile\n";
}
my ($readsfile, $qualitiesfile) = @ARGV;

open(READS, "$readsfile") or die "cannot open $readsfile\n"; 
open(QUALITIES, ">$qualitiesfile") or die "cannot open $qualitiesfile\n";  

my ($name, $read);
my $line;
$name = <READS>;
chomp($name);
my $i;
while (<READS>) {
  $line = $_;
  chomp($line);
  if ($line =~ />/) {
    print QUALITIES "$name\n";
    for $i (0 .. length($read)-1) {
      print QUALITIES "30 ";
      if ($i % 60 == 0 && $i > 0) {
	print QUALITIES "\n";
      }
    }
    if (length($read)-1 % 60 != 0 ) {
      print QUALITIES "\n";
    }
    $read = "";
    $name = $line;
  }
  else {
    $read = $read . $line;
  }
}
  
