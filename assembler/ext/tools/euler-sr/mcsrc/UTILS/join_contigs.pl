#!/usr/bin/env perl
if ($#ARGV != 3) {
  print "\nusage: $0 infile gaplen outfile line_length\n";
  exit(1);
}

$infile = shift @ARGV;
$gaplen = shift @ARGV;
$outfile = shift @ARGV;
$linelen = shift @ARGV;

open(IN, "$infile") or die "cannot open $infile\n";
open(OUT, ">$outfile") or die "cannot open $outfile\n";

$linenumber = 0;
$gap = '';
for $i (0 .. ($gaplen-1)) {
  $gap = $gap . 'N';
}
$curline = "";
while (<IN>) {
  $line = $_;
  chomp $line;
  if ($line =~ />/ and $linenumber != 0) {
    $curline = $curline . $gap;
  }
  else {
    if ($linenumber == 0) {
      print OUT "$line\n";
    }
    else {
      $curline = $curline . $line;
      if (length($curline) > $linelen) {
	$pline = substr($curline, 0, $linelen);
	print  OUT "$pline\n";
	$curline = substr($curline, $linelen);
      }
    }
  }
  $linenumber = $linenumber + 1;
}

close OUT;
close IN;  
