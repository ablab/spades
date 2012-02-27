#!/usr/bin/env perl
use strict;


if ($#ARGV ne 1) {
  print "usage: fixreads.pl readsfile paramfile\n";
  exit(0);
}

my $readsfile = shift @ARGV;
my $paramfile = shift @ARGV;


#  Get the parameters from the params file
open(PARAMFILE, "$paramfile") or die "cannot open $paramfile\n";
my @paramfile  = <PARAMFILE>;
close(PARAMFILE);

my @widths     = ();
my @thresholds = ();
my @ratios     = ();
my @deltas     = ();

my $exec = shift @paramfile;
chomp($exec);
if ($#paramfile < 0) {
  print "error, no instructions!\n";
  exit(0);
}
#     width thres. ratio (not quite perfect rexp)
grep (/(\d+) +(\d+) +(\d*\.{0,1}\d+) +(\d+)/ 
                 && push(@widths, $1) 
                 && push(@thresholds, $2) 
                 && push(@ratios, $3)
                 && push(@deltas, $4), @paramfile);

if ($#widths ne $#paramfile) {
  print "There is an error somewhere in the param file. Bailing out.\n";
#  exit(-1);
}

my $fixNo = "";
my $readsname;
my $outname;
my $i;
my $runstr;
my $txtout = "$readsfile." . $$;
if (-e "$txtout") {
`rm $txtout`;
}
print "creating $txtout\n";
`touch $txtout`;

foreach $i (0 .. ($#widths - 1)) {
  $readsname = $readsfile . "$fixNo";
  $fixNo = $i + 1;
  $outname   = $readsfile  . "$fixNo";
  $runstr = "csh -c \"$exec -c -i $readsname -k @widths[$i] -v @thresholds[$i] -r @ratios[$i] -D @deltas[$i] -o $outname  >> $txtout \" ";
  print "$runstr\n";
  `$runstr`;
}
$readsname = $readsfile . "$fixNo";
$fixNo = $#widths + 1;
$outname   = $readsfile  . "$fixNo";
$runstr = "csh -c \"$exec -c -i $readsname -k @widths[$i] -v @thresholds[$i] -r @ratios[$i] -o $outname -d \" ";
print "$runstr >> $txtout \n";
`$runstr`;
