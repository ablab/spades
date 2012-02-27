#!/usr/bin/env perl
use POSIX;

if ($#ARGV != 1) { 
  print "Make a biased markov model, with a qualitative bias level from 0 (no bias)\n";
  print "to 9 (lots of bias).  Bias is added iteratively by randomly selecting a \n";
  print "sequence, and increasng one of the transition probabilities (and decreasing \n";
  print " all others) with \n";
  print "the more likely transitions being more likely to be increased further\n";
  print "usage: $0 model-order bias-level (0-9)\n"; 
  exit(0);
}
$order = shift @ARGV;
$bias  = shift @ARGV;

# initialize the model
$modelSize = POSIX::pow(4,$order);
@model = ();
for $i (0 .. $modelSize-1) {
  @model[$i] = [];
  push @{@model[$i]}, 0.25, 0.25, 0.25, 0.25;
}

# go about making the model more biased

$step = 0.005;
for $i (0 .. $bias) {
  for $j (0 .. 100*$modelSize) {
    $transitionIndex = POSIX::floor(rand($modelSize));
    MakeDistMoreBiased(\@{@model[$transitionIndex]}, $step);
  }
}
for $j (0 .. $#model) {
  $str = IndexToStr($j, $order);
  printf("$str %.3f %.3f %.3f %.3f\n", @{@model[$j]}[0],@{@model[$j]}[1],@{@model[$j]}[2],@{@model[$j]}[3]);
}


sub IndexToStr {
  my($index, $order) = @_;
  my $str = "";
  my $i;
  my @trans = ('g', 'a', 'c', 't');
  while ($order > 0) {
    $i = $index % 4;
    $str = @trans[$i] . $str;
    $index = $index / 4;
    $order--;
  }
  return $str;
}

sub MakeCumDist {
  my ($dist, $cumDist) = @_;
  $cum = 0;
  @{$cumDist} = ();
  push @{$cumDist}, 0;
  for $d (@{$dist}) {
    $cum = $cum + $d;
    push @{$cumDist}, $cum;
  }
}

sub GenIndex {
  my ($dist) = @_;
  my $val = rand;
  my $i;
  for $i (0 .. $#{$dist}-1) {
#    print "$val @{$dist}[$i] @{$dist}[$i+1]\n";
    if ($val >=  @{$dist}[$i] and $val  <= @{$dist}[$i+1]) {
      return $i;
    }
  }
}

sub MakeDistMoreBiased {
  my ($dist, $step) = @_;
  if ($#{$dist} == 0) {
    return;
  }
  my @cumDist = ();
  MakeCumDist($dist, \@cumDist);
  my $increaseIndex = GenIndex(\@cumDist);
  my $i;
  my $increaseAmount = 0;
  $increaseAmount = (1 - @{$dist}[$increaseIndex]) * $step;
  @{$dist}[$increaseIndex] += $increaseAmount;

  $remainder = 0;
  for $i (0 .. $#{$dist}) {
    if ($i != $increaseIndex) {
      $remainder += @{$dist}[$i];
    }
  }
  for $i (0 .. $#{$dist}) {
    if ($i != $increaseIndex) {
      @{$dist}[$i] -= $increaseAmount * @{$dist}[$i] / $remainder;
    }
  }  
 # printf("new dist: %.3f %.3f %.3f %.3f %.3f\n", $increaseAmount, @{$dist}[0], @{$dist}[1], @{$dist}[2], @{$dist}[3]);
}
