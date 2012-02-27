#!/usr/bin/env perl

use Bio::SeqIO;
use Bio::Seq;

if ($#ARGV < 1 ) {
  print "usage: JoinNISCContigs.pl alignfile outfile [displayid]\n";
  print "alignfile should reference files in the curernt directory\n";
  exit(0);
}
$alignfile = shift @ARGV;
$outfile   = shift @ARGV;


$displayId  = "";
if ($#ARGV >= 0) {
  $displayId = shift @ARGV;
}


$seqOut = Bio::SeqIO->new(-file=>">$outfile", -format=>"fasta");

if (! -e $alignfile or -z $alignfile) {
  print "warning: no alignments for $alignfile\n";
  if ($displayId != "") {
    @faNames = glob("*.$displayId.*.fasta");
    $first = 1;
    print "got sequences @faNames\n";
    foreach $seqName (@faNames) {
      $seqIn = Bio::SeqIO->new(-file=>$seqName, -format=>"fasta");
      if ($first == 1) {
	$fullSeq =  $seqIn->next_seq;
	$first = 0;
      }
      else {
	$next = $seqIn->next_seq;
	$fullSeq = Bio::Seq->new(-seq=>$fullSeq->seq . $gap . 
				 $next->seq);
      }
    }  
    $fullSeq->display_id($displayId);
    $seqOut->write_seq($fullSeq);
    exit(0);
  }
}


open(AF, $alignfile);

$firstSeq = 1;
$ref = 0;
$qry = 1;

$gap = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" .
  "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" .
  "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" .
  "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" .
  "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";


while(<AF>) {
  $line = $_;
  @alignSpec = split(/\s+/, $line);
  if ($firstSeq) {
    $seqIn = Bio::SeqIO->new(-file=>"$alignSpec[0]", -format=>"fasta");
    $fullSeq = $seqIn->next_seq;
    $firstSeq = 0;
  }
  $seqIn  = Bio::SeqIO->new(-file=>"$alignSpec[3]", -format=>"fasta");
  $qrySeq = $seqIn->next_seq;

  if ($alignSpec[1] == -1) {
    # the adjacent bacs have an overlap
    $fullSeq = Bio::Seq->new(-seq=>$fullSeq->seq . $gap . $qrySeq->seq);
  }
  else {
    # append the non-overlapping part of the query sequence
    if ($alignSpec[5] > $qrySeq->length) {
	$start = $qrySeq->length;
    }
    else {	
        $start = $alignSpec[5];
    }
    $fullSeq = Bio::Seq->new(-seq=>$fullSeq->seq .
			     $qrySeq->subseq($start,$qrySeq->length));
    if ($fullSeq) {
      $len = $fullSeq->length();
    }
  }
  # swap which is the ref and which is qry, cheap pointers
}

$fullSeq->display_id($displayId);
$seqOut->write_seq($fullSeq)
