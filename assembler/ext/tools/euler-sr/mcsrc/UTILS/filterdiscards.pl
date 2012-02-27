#!/usr/bin/perl

if ($#ARGV != 2) {
  print "usage: $0 readsfile discardsfile newfile\n";

$reads = @ARGV[0];
$discardsfile = @ARGV[1];
$newfile = @ARGV[2];


open(READS, "$reads") or die "cannot open $reads\n";
open(DISCARDS, "$discardsfile") or die "cannot open $discardsfile\n";
open(NEW, ">$newfile") or die "cannot open $newfile for write\n";
@discards = <DISCARDS>;
$indices  = {};
grep ((/\D+(\d+)/) and ($indices{$1} = 1)), @discards;

while (<MATES>) {
  $line = $_;
  $line =~ /\w+ (\d+) (\d+)/;
  $a = $1 - 1;
  $b = $2 - 1;
  if ( (not $indices{$a}) and (not $indices{$b})) {
    print NEW $line;
  }
}


    
  
