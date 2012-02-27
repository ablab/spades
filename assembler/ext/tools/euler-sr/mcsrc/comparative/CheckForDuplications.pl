#!/usr/bin/env perl

use common;
use Bio::SeqIO;
use Bio::Seq;

$locusFile = shift @ARGV;
$dbDir = ".";
$seqDir = ".";
$verbose = 0;

while ($#ARGV >= 0) {
  $opt = shift @ARGV;
  if ($opt eq "-seqdir") {
    $seqDir = shift @ARGV;
  }
  if ($opt eq "-db") {
   $dbDir     = shift @ARGV;
  }
  elsif ($opt eq "-verbose") {
   $verbose = 1;
  }
}
my @loci;
common::ReadLocusFile($locusFile, \@loci);

$numDuplications = 0;
my @exc;
for($l = 0; $l <= $#loci; $l++) {
  $species = $loci[$l][0];
  $start   = $loci[$l][1];
  $end     = $loci[$l][2];
  $sbjctName = "$species.fasta";
  $seqIn = Bio::SeqIO->new('-format'=>"fasta",
			   '-file'=>"$seqDir/$sbjctName");
  $sequence = $seqIn->next_seq;

  if ($start < 0 or 
      $end > $sequence->length) {
    print "bad locus: $species $start $end\n";
  }
  else {
    $qrySeqStr = $sequence->subseq($start, $end);
    $qrySeqStr =~ tr/actg/N/;
    $qrySeq = Bio::Seq->new(-seq=>$qrySeqStr);
    $ql = $qrySeq->length;
    $qs = $qrySeq->seq;
    # now search for this
    $database = "$dbDir/$sbjctName";
    my @hits;
    @hits = ();
    $numHits = common::GetHitCoordinates($qrySeq, $database, "1e-10", \@hits);
    @coverage = ();
    $qsl = length $qrySeqStr;
    for ($c = 0; $c < $qsl; $c++ ) {
     push @coverage, 0;
    }
    for ($h = 0; $h < $numHits; $h++ ) {
      if ($hits[$h][0] < $hits[$h][1]) { 
        $hitStart = $hits[$h][0]; 
        $hitEnd   = $hits[$h][1];
      }
      else { 
        $hitStart = $hits[$h][1];
        $hitEnd   = $hits[$h][0];
      } 
      for ($p = $hitStart; $p <= $hitEnd; $p++) {
	$coverage[$p]++;
      }
    }
    # count the number  of alignable nucleotides
    $numAlignable = $qrySeqStr =~ tr/ACTG//;
    # compare the number alignable with the number masked
    $numAligned = 0;
    for ($p = 0; $p < $qrySeq->length; $p++ ) {
	$numAligned += $coverage[$p];
    }
    # generate a report
    $cRat = $numAligned *1.0 / $numAlignable;
    if ($verbose) {
	print "$locusFile $species $cRat $numAligned $numAlignable\n";
    }
    if ($cRat  > 1.25) {
      $numDuplications++;
      push @exc, $cRat;
     push @dup, $species[$l];
    }
  }
}
$numLoci = scalar @loci;
if ($numDuplications > 0) {
  print "WARNING, locus: $locusFile contains $numDuplications duplicated loci from $numLoci, it may not be reliable ";
  for ($i = 0; $i <= $#exc; $i++) {
    print "$dup[$i] $exc[$i] ";
  }
  print "\n";
}
