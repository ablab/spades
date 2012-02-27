#!/usr/bin/env perl

$in = shift @ARGV;

open (IN, "$in") or die "cannot open $in\n";


%names = ();
$numSpecies = 0;
@locations = ();
%indices = ();
while (<IN>) {
  $line = $_;
  if ($line !~/didn.*/) {
    $line =~ /(\w+) (\w+) (\w+) (\d+) (\d+)/;
    $orig = $1;
    $spec = $3;
    $loc  = $4;
    if (not exists $indices{$spec}) {
      ++$numSpecies;
      $indices{$spec} = $numSpecies;
    }
    if (not exists $indices{$orig}) {
      ++$numSpecies;
      $indices{$orig} = $numSpecies;
    }
    print "s $indices{$orig} $indices{$spec} $orig $spec\n";
    print "$loc\n"
  }
}

#for $i (keys %names) {
#  print "$i @{$names{$i}} \n";
#}
