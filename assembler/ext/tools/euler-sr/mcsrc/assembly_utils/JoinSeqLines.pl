#!/usr/bin/env perl

$in = shift @ARGV;
$out = shift @ARGV;
open(IN, "$in") or die "cannot open $in\n";
open(OUT, ">$out") or die "cannot open $out\n";
$seq ="";
while(<IN>) {
 if (/^>/) {
				if ($seq ne "") {
			  	print OUT $title;
			  	print OUT "$seq\n";
        }
				$seq = "";
				$title = $_;
  }
  else {
				chomp;
				$seq .= $_;
				}
  }

print OUT "$title";
print OUT "$seq\n";
   
       
