#!/usr/bin/env perl
if ($#ARGV != 2) {
  print "usage: $0 filepattern sequencefile outfile\n";
  exit(0);
}

$filepattern = shift @ARGV;
$sequenceFile = shift @ARGV;
$outfile = shift @ARGV;

@files = glob("$filepattern");
%fh = ();
foreach $file (@files) {
  $file =~ /.*\.(\d+)\..*/;
  $number = $1;
  $fh{$number} = $file;
}

@ks = sort {$a <=> $b} keys %fh;
print "got keys: @ks\n";
open (SF, "$sequenceFile") or die "cannot open $sequenceFile\n";
$s = <SF>;

open (OUT, ">$outfile") or die "cannot open $outfile\n";
print OUT "$s";

for $k (@ks) {
  print "opening $fh{$k}\n";
  open (IN, $fh{$k}) or die "cannot open $fh{$k}\n";
  @l = <IN>;
  print OUT "bin: $k 0 to 0\n";
  foreach $line (@l) {
    if ($line =~ /(\w+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S.*)/) {
      $spec = $1; $start = $2; $end = $3; $nspec = $4; $rest = $5;
      @dirs = split /\s+/, $rest;
      $nspec = scalar @dirs;
      print OUT "$spec $start $end $nspec $rest\n";
    }
  }
  close IN;
}

close OUT;

