#!/usr/bin/env perl

use ReadLibrary;

# program to read through a fasta file and determine which reads are mate-pairs.
# This assumes only 2 naming conventions:
#
# ABI long:  000000532223A01.F.ab1
#
# ABI short: S401149.PSM-578-M13R.abl
#
# Read pairs from the same clone have the same name, except for forward and reverse.
# 

if ($#ARGV != 0) {
  print "usage: $0 reads_file\n";
  print "searches through reads_file for mate-pair information encoded in the \n";
  print "read names. \n";
  print "Read names that are parsed right now are: \n";
  print "/S401149\.(.{6,16})([FR])/\n";
  print "example ( add example ) \n";
  print "/0+(.{8,12})\.([FR])\.ab1/\n";
  print "example ( add example ) \n";
  print "/000000(.{4,10})@.{4,10}([FR]).\.[bg]uf/\n";
  print "example ( add example ) \n";
  exit(0);
}

$infile = shift @ARGV;

open(IN, $infile) or die "cannot open $infile\n";
%pairs = ();
$index = 0;
@titles = ();
@dirs   = ();
while (<IN>) {
  $line = $_;
  if ($line =~/^>(.*)/) {
    $title = $1;
    $dir = -1;
    $base = "";
#    print "matching $title";
    ($base, $dir, $type) = ReadLibrary::ParseABITitle($title);
#   print "got bdt: $base, $dir, $type\n";
    $side = GetSide($dir);
    if (base ne "" and ! exists $pairs{$base} ) {
      @{$pairs{$base}} = (-1, -1, $type);
    }
    else {
      # sanity check, both mate pairs should not be assigned.
      if ($pairs{$base}[0] != -1 and
	  $pairs{$base}[1] != -1) {
	print "error, template $base exists more than twice at $index\n";
	print "$pairs{$base}[0] != -1 $pairs{$base}[1] != -1 \n";
	exit(0);
      }
    }
    if ($base ne "" and $side != -1) {
      $pairs{$base}[$side] = $index;
    }
    $index++;
    push @titles, $base;
    push @dirs, $side; 
  }
}

for $i (0 .. $#titles) {
#  print "$titles[$i] $dirs[$i]\n";
  print "$pairs{$titles[$i]}[Other($dirs[$i])]\n";
}

sub Other {
  my ($val) = @_;
  if ($val == 0) { return 1; }
  if ($val == 1) { return 0; }
  return -1;
}

sub GetSide {
  my ($side) = @_;
  if ($side eq "F") { return 0; }
  if ($side eq "R") { return 1; }
  return -1;
}
