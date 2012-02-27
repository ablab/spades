#!/usr/bin/env perl
$prefix = shift @ARGV;
@files = glob("$prefix.*");
%dim1 = ();
%dim2 = ();
$min1 = 99999999;
$max1 = 0;
$min2 = 9999999999;
$maxx = 0;
$nfiles = scalar @files;
foreach $f (@files) {
  $f =~ /$prefix\.(\d+)\.(\d+).*/;
  $dim1{$1} = 1;
  $dim2{$2} = 1;
  if ($1 < $min1) {
	   $min1 = $1;
  }
  if ($1 > $max1) {
	  $max1 = $1;
  }
  if ($2 < $min2) {
    $min2 = $2;
  }
  if ($2 > $max2) {
    $max2 = $2;
  }
}

@d1 = sort {$a <=> $b} keys %dim1;
@d2 = sort {$a <=> $b} keys %dim2;
$nd1 = scalar @d1;
$nd2 = scalar @d2;
$min1 = shift @d1;
$max1 = pop @d1;
print "$nd1 $min1 $max1 $nd2 $min2 $max2\n";
