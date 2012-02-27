#!/usr/bin/env perl

$0 =~ /.*\/([^\/].*)/;
$progname = $1;

if ($#ARGV < 2) {
  print "markov-model generate sequence\n";
  print "generate a sequence according to an n'th order markov-model\n";
  print "usage: $progname markov_model length outfile\n";
  exit(0);
}

$mmfile = shift @ARGV;
$length = shift @ARGV;
$outfile = shift @ARGV;

%model = ();
%cumDistModel = ();
ReadMarkovModel($mmfile, \%model);

GenCumDist(\%model, \%cumDistModel);
for $c (keys %cumDistModel) {
  print "cdm: $c  @{$cumDistModel{$c}} \n";
}

$seq = "";

GenerateSequence(\%cumDistModel, $length, \$seq);

PrintSeq($outfile, $seq, 60);

sub ReadMarkovModel {
  ($mmfile, $model) = @_;
  
  open(MMFILE, "$mmfile") or die "cannot open $mmfile\n";
  while (<MMFILE>) {
    $line = $_;
    $line =~ /([actg]+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
    $g = $2; $a = $3; $c = $4; $t = $5;
    @{${$model}{$1}} = ();
    push @{${$model}{$1}}, $g, $a, $c, $t;
  }
  close MMFILE;
  # quick sanity check
  @seqns = keys %{$model};
  if ($#seqns <= 0) {
    print "error processing $mmfile, no data\n";
    exit(0);
  }
  $l0 = length @seqns[0];
  for $seq (@seqns) {
    if (length($seq) != $l0) {
      print "error in $mmfile, length of $seq is different from @seqns[0]\n";
      exit(0);
    }
  }
}

sub GetOrder {
  my ($model) = @_;
  @seqns = keys %{$model};
  return length(@seqns[0]);
}

sub GenCumDist {
  my ($model, $cumDistModel) = @_;
  my $cum;
  for $seq (keys %{$model}) {
    @{${$cumDistModel}{$seq}} = ();
  $cum = 0;
  $l = $#{${$model}{$seq}}+1;
  for $i (0 .. $#{${$model}{$seq}}) {
      $cum = $cum + @{${$model}{$seq}}[$i];
      push @{${$cumDistModel}{$seq}}, $cum;
    }
  }
}

sub GenNuc {
  my ($dist) = @_;
  my $n = rand;
  if ($n >= 0 and $n <= @{$dist}[0]) {
    return 'g';
  }
  elsif ($n >= @{$dist}[0] and $n <= @{$dist}[1]) {
    return 'a';
  }
  elsif ($n >= @{$dist}[1] and $n <= @{$dist}[2]) {
    return 'c'
  }
  elsif ($n >= @{$dist}[2] and $n <= @{$dist}[3]) {
    return 't';
  }
}

sub GenerateSequence {
  ($model, $length ,$seq) = @_;
  @equal = (0.25, 0.50, 0.75, 1.0);
  my $order = GetOrder($model);
  # generate the part of the sequence that has no
  # corresponding markov model
  ${$seq} = "";
  my $n;
  for $i (0 .. $order-1) {
    $n = GenNuc(\@equal);
    ${$seq} = ${$seq} . GenNuc(\@equal);
  }
  $length = $length - $order;
  my $tail;
  for $i (0 .. $length) {
    $tail = substr ${$seq}, -1*$order;
    $n = GenNuc(\@{${$model}{$tail}});
    ${$seq} = ${$seq} . $n;
  }
}

  
  
sub PrintSeq {
  my ($outfile, $seq, $len) = @_;
  open (OUT, ">$outfile") or die "cannot open $outfile\n";
  my $index = 0;
  while ($index < length (${$seq})) {
    $subseq = substr ${$seq}, $index, $len;
    print OUT "$subseq\n";
    $index = $index + $len;
  }
  close OUT;
}
    
