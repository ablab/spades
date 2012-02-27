#!/usr/bin/env perl


if ($#ARGV < 0) {
  printf("usage: $0 : filterfile\n");
  exit(-1);
}

my $fileName= @ARGV[0];

open(FITFILE, "$fileName") or die "cannot open $fileName\n";


# ditch first 3 lines.
<FITFILE>;
<FITFILE>;
<FITFILE>;

my ($fstart, $fend);
my ($rstart, $rend);
my $prev;
while (<FITFILE>) {
  $line = $_;
  if ($line =~ /rev pos (\d+) (\d+)/) {
    $rstart = $1;
    $rend   = $2;
    printf("%d\n", $rstart - $fend);
    if ($rstart - $fend < 0) {
      print "$line";
      print "$prev"
    }
  }
  else {
    $line =~/pos (\d+) (\d+)/;
    $fstart = $1;
    $fend   = $2;
    $prev = $line;
  }
  
}
