#!/usr/bin/env perl

use POSIX;
if ($#ARGV < 0 ) {
  print "\nusage: $0 bedfile\n";
  print " prints the coordinates in the bedfile in one file per chromosome\n\n";
  exit(0);
}

open(IN, @ARGV[0]) or die "cannot open @ARGV[0]\n";
shift @ARGV;
$ratio = 1;
if ($#ARGV >= 0) {
  $ratio = shift @ARGV;
}
print "ratio: $ratio\n";
%hCoords = ();
%cCoords = ();

while(<IN>) {
  $line = $_;
  ($human, $chimp) = split(/\s+/, $line);
  $human =~ /(chr\S+):(\d+)\-(\d+)/;
  $hChr = $1; $hS = $2; $hE = $3;
  if (not exists $hCoords{$hChr}) {
    @{$hCoords{$hChr}} = ();
  }
  push @{$hCoords{$hChr}}, ($hS, $hE);

  $chimp =~ /(chr\S+):(\d+)\-(\d+)/;
  $cChr = $1; $cS = $2; $cE = $3;
  if (not exists $cCoords{$cChr}) {
    @{$cCoords{$cChr}} = ();
  }
  push @{$cCoords{$cChr}}, ($cS, $cE);
}

foreach $chr (keys %hCoords) {
  $hFile = "human.$chr";
  open(H, ">$hFile") or die "cannot open $hFile\n";
  for ($h = 0; $h < $#{$hCoords{$chr}}; $h += 2) {
    $length = $hCoords{$chr}[$h+1] - $hCoords{$chr}[$h] + 1;
    $start = $hCoords{$chr}[$h] - POSIX::floor($length*$ratio);
    $end   = $hCoords{$chr}[$h+1] + POSIX::floor($length*$ratio);
    print H "$start $end\n";
  }
  close H;
}

foreach $chr (keys %cCoords) {
  $cFile = "chimp.$chr";
  open(C, ">$cFile") or die "cannot open $cFile\n";
  for ($c = 0; $c < $#{$cCoords{$chr}}; $c += 2) {
    $length = $cCoords{$chr}[$c+1] - $cCoords{$chr}[$c] + 1;
    $start = $cCoords{$chr}[$c] - POSIX::floor($length*$ratio);
    $end   = $cCoords{$chr}[$c+1] + POSIX::floor($length*$ratio);
    print C "$start $end\n";
  }
  close C ;
}
