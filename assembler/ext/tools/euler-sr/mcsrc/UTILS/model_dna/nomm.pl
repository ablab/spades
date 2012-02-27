#!/usr/bin/env perl

use POSIX;

if ($#ARGV < 0) {
  print "N-th Order Markov Model.\n";
  print "find the distribution of nucleotides given the previous n-nucleotides\n";
  print "usage: nomm.pl infile [order=1]\n";
  exit 0;
}

$inputFile = shift @ARGV;
$order = 1;
if ($#ARGV >= 0) {
  $order = shift @ARGV;
}

%nucs = {};

$nucs{'G'} = 0;
$nucs{'g'} = 0;
$nucs{'A'} = 1;
$nucs{'a'} = 1;
$nucs{'C'} = 2;
$nucs{'c'} = 2;
$nucs{'T'} = 3;
$nucs{'t'} = 3;

@revIndex =  ('g', 'a', 'c', 't');
$countSize = POSIX::pow(4, $order);

for $c (0 .. $countSize-1) {
  @count[$c] = 0;
}

open(IN, "$inputFile") or die "cannot open $inputFile\n";
$l = <IN>;
@seq = <IN>;
close IN;
chomp @seq;
$dna = join "", @seq;
$dna =~ tr/ACTG/actg/;
# compute the n-nucleotide frequencey
for $i (0 .. (length($dna) - $order - 1)) {
  $ostr = substr($dna, $i, $order);
  $index = GetIndex($ostr, \%nucs);
  @count[$index]++;
}

%margins = {};
$ldna=  length($dna);
$total = 0;
%dist = {};
for $idx (0 .. $countSize-1) {
  $str = IndexToStr($idx, $order, \@revIndex);
  $distr{$str} = @count[$idx] / $ldna;
}


CalcMarginalP(\%distr, \%margins);

$countSize = POSIX::pow(4,$order-1);
for $idx (0 .. $countSize-1) {
  $str = IndexToStr($idx, $order-1, \@revIndex);
  $seqg = $str . 'g';
  $seqa = $str . 'a';
  $seqc = $str . 'c';
  $seqt = $str . 't';
  printf("$str %.3f %.3f %.3f %.3f\n",$margins{$seqg}, $margins{$seqa}, $margins{$seqc}, $margins{$seqt});
}

sub IndexToStr {
  my($index, $order, $trans) = @_;
  my $str = "";
  my $i;
  while ($order > 0) {
    $i = $index % 4;
    $str = @{$trans}[$i] . $str;
    $index = $index / 4;
    $order--;
  }
  return $str;
}

sub GetIndex {
  my($dna, $trans) = @_;
  $len = length($dna);
  $index = 0;
  for $i  (0..$len-1) {
    $nuc = substr($dna, $i, 1);
    $index = $index * 4;
    $index = $index + ${$trans}{$nuc};
  }
  return $index;
}


sub CalcMarginalP {
  my ($prob, $margprob) = @_;
  
  @seqns = keys %{$prob};
  my $order = length(@seqns[0]);
  for $seq (@seqns) {
    $base = substr $seq, 0, $order-1;
    $marginal = 0;
    for $nuc ( ('a', 'c', 't', 'g')) {
      $margSeq = $base . $nuc;
      $marginal += ${$prob}{$margSeq};
    }
    if ($marginal != 0) {
      ${$margprob}{$seq} = ${$prob}{$seq} / $marginal;
    }
    else {
      ${$margprob}{$seq} = 0;
    }
  }
}
  
  
  

  

