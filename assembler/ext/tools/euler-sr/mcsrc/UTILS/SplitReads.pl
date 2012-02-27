#!/usr/bin/env perl

if ($#ARGV != 3) {
  print "usage: $0 reads reads-per-file outbase\n";
  print "       Divide the file 'reads' into n files, each file has 'reads-per-file' in it\n";
  exit(0);
}

$readsFile = shift @ARGV;
$readsPerFile = shift @ARGV;
$outBase  = shift  @ARGV;
$ext      = shift @ARGV;

$r = 0;
$fileNumber = 0;
open(RF, "$readsFile") or die "cannot open $readsFile\n";
while(<RF>) {
  $line = $_;
  if ($line =~ />/) {
    if ($r % $readsPerFile == 0) {
      $outName = "$outBase.$fileNumber.$ext";
      if ($r > 0) {
        close OUT;
      }
      open(OUT, ">$outName") or die "cannot open $outName\n";    
      $fileNumber++;
    }
    ++$r;
  }
  print OUT $line;
}
