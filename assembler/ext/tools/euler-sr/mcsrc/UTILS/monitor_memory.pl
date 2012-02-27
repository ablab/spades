#!/usr/bin/env perl

if ($#ARGV < 1) {
  print "usage: $0 progname delay outfile\n";
  print "repcat will output the amount of memory one instance of prog 'progname'\n";
  print "is using.\n";
  exit(0);
}

  
$progname = shift @ARGV;
$delay    = shift @ARGV;
if ($#ARGV >= 0) {
  $outfile  = shift @ARGV;
  open (OUT, ">$outfile") or die "cannot open $outfile\n";
}
else {
  *OUT = *stdout;
}

if ($delay == 0) {
  print "cannot have 0 delay.  A value of $delay was specified\n";
  exit(0);
}



$cont = 1;
while ($cont == 1) {
  $topstr = `top`;
  @lines = split /'\n'/, $topstr;
  $line = "";
  grep /(.*$progname.*)/ && ($line = $1), @lines;
  if ($line eq "") {
    $cont = 0;
  }
  @stats = split /\s+/, $line;
  if ($cont != 0) {
    print OUT "using memory @stats[4] @stats[5] \n";
  }
  sleep($delay);
}

close(OUT)
