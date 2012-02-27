#!/usr/bin/env perl

use grimm;
if ($#ARGV < 0) {
  print "usage: $0 mgr_micro.txt\n";
  exit(1);
}

$mgrMicroIn =shift @ARGV;
open(MGRMICRO, $mgrMicroIn) or die "cannot open $mgrMicroIn\n";

@mgrFile = <MGRMICRO>;
@blocks = ();
@order  = ();

grimm::ParseMicroGrimmBlocks(\@mgrFile, \@blocks, \@order);

$nb = scalar @blocks;
$no = scalar @order;
for ($b = 0; $b < $nb ; $b++ ){
  $ns = scalar @{$blocks[$b]};
  @inv = ();
  grimm::FindInversions(\@{$order[$b][0]}, \@{$order[$b][1]}, \@inv);
  @anchors = ();
  grimm::GetInvertedAnchors(\@{$blocks[$b]}, \@inv, \@anchors);
  for ($a = 0; $a <= $#anchors; $a++ ) {
    $s1 = $anchors[$a]{"start1"};
    $l1 = $anchors[$a]{"len1"};
    $e1 = $s1 + $l1;
    $s2 = $anchors[$a]{"start2"};
    $l2 = $anchors[$a]{"len2"};
    $e2 = $s2 + $l2;
    print "$s1\t$e1\t$s2\t$e2\n";
  }
}


