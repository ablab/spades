#!/usr/bin/env  perl

if ($#ARGV != 1) {
  print "usage: $0 macro_blocks_file blocks_file\n";
  exit(0);
}

$orderFileName = shift @ARGV;

open(orderFile, $orderFileName) or die "cannot open $orderFileName\n";


# assumptions about the input file:
# 1. there are two genomes
# 2. there is only one chromosome per genome

$found = 0;
while (<orderFile>){
  if ($_ =~ /genome2/) {
    $found = 1;
    last;
  }
}

if ($found == 0) {
  die "did not find genome 2\n";
}

# discard the name of the chromosome, but check
$chr = <orderFile>;
if ($chr !~ /# Chromosome 1/) {
  die "did not find Chromosome 1\n";
}

$permLine = <orderFile>;

@perm = split(/\s+/, $permLine);

@inverted = ();
# now find micro-inversions
for ($i = 1; $i < $#perm ; $i++) {
  $s = sign(@perm[$i]);
  if (@perm[$i-1] + 1 == -1*@perm[$i] and
      (-1 * @perm[$i]) + 1 == @perm[$i+1] and
     sign(@perm[$i]) == -1) {
    push @inverted, -1*@perm[$i];
  }
}



# now read what blocks were present

$blocksFileName = shift @ARGV;
open(blocksFile, $blocksFileName) or die "cannot open $blocksFileName\n";
@blocks = ();
while(<blocksFile>) {
  @vals = split(/\s+/, $_);
  push @blocks, [@vals];
}

$firstBlockIndex = $blocks[0][0];

# print the inverted blocks
for ($i =0; $i <= $#inverted; $i++) {
  print "@{@blocks[$inverted[$i] - $firstBlockIndex]}\n";
}


sub sign {
  my ($val) = @_;

  if ($val < 0) {
    return -1;
  }
  if ($val > 0) {
    return 1;
  }
  else {
    return 0;
  }
}


