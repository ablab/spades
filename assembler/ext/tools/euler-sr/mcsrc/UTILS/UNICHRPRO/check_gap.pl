#!/usr/bin/env perl
use IO::Handle;

if ($#ARGV != 1 && $#ARGV != 4) {
  print "usage: $0 tailsfile extractsfile [ lengthThresh minLength locThresh ]\n";
  exit(0);
}
$tailsFile = shift @ARGV;
$extractsFile = shift @ARGV;

$lengthThresh = 500;
$locThresh =  1200;
$minLength = 400;
if ($#ARGV == 2) {
  $lengthThresh = shift @ARGV;
  $minLength    = shift @ARGV;
  $locThresh    = shift @ARGV;
}

open(TAILS, "$tailsFile") or die "cannot open $tailsFile\n";
open(EXTRACTS, "$extractsFile") or die "Cannot open $extractsFile\n";

@tails = <TAILS>;
@extracts = <EXTRACTS>;

$tailFile = "tail.tmp.$$";
$extractFile = "extract.tmp.$$";
while ($#tails >= 0 && $#extracts >= 0) {
  ($tailTitle, $tailSeq) = GetSeq(\@tails); 
  ($extractTitle, $extractSeq) = GetSeq(\@extracts);

  open(TAIL, ">$tailFile") or die "cannot open tail.tmp\n";
  open(EXTRACT, ">$extractFile") or die "cannot open extract.tmp\n";
  
#  print "comparing $tailTitle $extractTitle\n";
  print TAIL $tailTitle;
  print TAIL $tailSeq;
  close TAIL;
  
  print EXTRACT $extractTitle;
  print EXTRACT $extractSeq;
  close EXTRACT;
  
  $bl2seqRes = `/home/mchaisso/blast/bl2seq -i $tailFile -j $extractFile -p blastn -D 1`;
  @bl2seqLines = split /\n/, $bl2seqRes;
  $possibleMatch = 1;
  $maxMatch = 0;
  $matchLocation = 0;
  foreach $bl2res (@bl2seqLines) {
    @bl2resVals = split /\s+/, $bl2res;
    $matchlen  = @bl2resVals[3];
    $tailEnd = @bl2resVals[7];
    $tailStart = @bl2resVals[6];
#    print "$matchlen $tailEnd $tailStart\n";
    if ($maxMatch < $matchlen) {
      $maxMatch = $matchlen;
      $matchLocation = $tailStart;
    }
    if ($matchlen > $lengthThresh) {
      if ($tailEnd > $locThresh) {
	print "$bl2res\n";
	$possibleMatch = 0;
      }
    }
  }
  if ($maxMatch < $minLength) {
    $possibleMatch = 0;
  }
  chomp $tailTitle;
  chomp $extractTitle;
  print "$possibleMatch $tailTitle|$extractTitle ($maxMatch, $matchLocation)\n";
}
  

sub GetSeq {
  my ($seqs) = @_;
  $title = shift @{$seqs};
  my($seq) = "";
  $continue = 1;
  while ($#{$seqs} >= 0 && $continue) {
    $line = shift @{$seqs};
    if ($line =~ /^>/) {
      unshift @{$seqs}, $line;
      $continue = 0;
    }
    else {
      $seq = $seq . $line;
    }
  }
  return ($title, $seq);
}

  

