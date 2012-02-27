#!/usr/bin/env perl

if ($#ARGV < 2 ) {
  print "usage: $0 mapfile suffix alignment_files \n";
  exit(0);
}

$mapFile = shift @ARGV;
open (MAP_FILE, "$mapFile") or die "cannot open map: $mapFile \n";

%names = {};

while (<MAP_FILE>) {
  $line = $_;
  $line =/(\w+)\s+(\d+)/;
  $names{$1} = $2;
}

$suffix = shift @ARGV;


foreach $file (@ARGV) {
  $file =~ /(\w+)\.(\w+)\W*/;
  $spec1 = $1;
  $spec2 = $2;
  open (GLUE_FILE, "$file") or die "cannot open $file\n";
  open (TRANS_GLUE_FILE, ">$file.$suffix") or die "cannot open $file.$sufix\n";
  while (<GLUE_FILE>) {
    $line = $_;
    if ($line =~ /(\w+)\s*(@);\s*(\w+)\s*(\!)\s*;/) {
      print TRANS_GLUE_FILE "$1 $names{$spec1}; $3 $names{$spec2};\n";
    }
    else {
      print "error parsing glue file\n";
      exit(0);
    }
  }
  close (TRANS_GLUE_FILE);
  close (GLUE_FILE);
}
