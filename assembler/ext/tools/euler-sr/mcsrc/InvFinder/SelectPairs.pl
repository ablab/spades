#!/usr/bin/env perl

#!/usr/bin/env perl
use POSIX;

if (@ARGV <  2) {
  print "usage: InvToPosFile.pl inffile ref1 [ref2 ...] \n";
  exit(0);
}

$in = shift @ARGV;
%relSpecies = ();
$index = 0;
while ($#ARGV >= 0) {
  $relSpecies{shift @ARGV} = $index;
  $index++;
}


open(IN, "$in") or die "cannot open $in\n";
@rlk = keys %relSpecies;


$ref = "";
$qry = "";
$doPrint = 0;
while (<IN>) { 
  $line = $_;
  if ($line =~ /EN/) {
    # line is a lav title
    if ($line =~ /(\w+)\.\w+\.fa\.(\w+)\..*/) {
      $ref = $1;
      $qry = $2;
#      print "ref: $ref qry: $qry\n";
      $doPrint = 0;
      if (exists($relSpecies{$ref}) && 
	  exists($relSpecies{$qry})) {
	if ($relSpecies{$ref} < $relSpecies{$qry}) {
	  $doPrint = 1;
	  print $line;
	}
      }
      else {
	$doPrint = 0;
      }
    }
    else {
      print "malformed title line $line\n";
      exit(0);
    }
  }
  else {
    if ($line =~ /(\d+) (\d+) (\d+) (\d+) (\d)/) { 
      if ($doPrint) {
	print $line;
      }
    }
    else {
      print "malformed coordinates line $line\n";
      exit(1);lformed
    }
  }
}
