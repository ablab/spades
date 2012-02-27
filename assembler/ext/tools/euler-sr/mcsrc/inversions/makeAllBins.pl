#!/usr/bin/env perl

if ($#ARGV != 4) {
  print "usage: $0 filenamePattern sequencefile orthfile seqdir outfile\n";
  exit 0;
}

$pattern = shift @ARGV;
$sequences = shift @ARGV;
$orth = shift @ARGV;
$seqDir = shift @ARGV;
$out  = shift @ARGV;

if (-e $out) {
  `mv $out $out.bak`;
}
print "pattern: $pattern\n";
@files = glob($pattern);

%ku = ();
foreach $f (@files) {
  $f =~ /bin\.(\d+)\..*/;
  $ku{$1} = $f;
}

@ks = sort {$a <=> $b } keys %ku;


$nk = scalar @ks;
print "running on $nk files\n";
$k = shift @ks;
$f = $ku{$k};
if (-e $out) {
  `mv $out $out.bak`;
}
$cmd = "~/projects/mcsrc/inversions/make_grid.pl $f  $sequences $orth $seqDir $out";

print "running $cmd\n";

system("$cmd");
$nSeq = CountSequences($f);

$totSeq = $nSeq;
foreach $k (@ks) {
  $f  = $ku{$k};
  $cmd = "~/projects/mcsrc/inversions/make_grid.pl $f $sequences $orth $seqDir $out";
  print "running $cmd\n";
  system("$cmd ");
}

sub CountSequences {
  my ($fileName) = @_;
  $nSeq = 0;
  open (IN, "$fileName");
  while (<IN>) {
    $line = $_;
    if ($line =~ /^>/) {
      ++$nSeq;
    }
  }
  return $nSeq;
}
