#!/usr/bin/env perl

use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;


$gap = 2;

$match = -1;
$mismatch = 1;

%score = (
	  'G' => {
		  'G' => $match,
		  'g' => $match,
		  'a' => $mismatch,
		  'A' => $mismatch,
		  'C' => $mismatch,
		  'c' => $mismatch,
		  'T' => $mismatch,
		  't' => $mismatch,
		  'n' => $mismatch,
		  'N' => $mismatch,
		  '.' => $gap},
	  'g' =>  {
		   'G' => $match,
		   'g' => $match,
		   'a' => $mismatch,
		   'A' => $mismatch,
		   'C' => $mismatch,
		   'c' => $mismatch,
		   'T' => $mismatch,
		   't' => $mismatch,
		  'n' => $mismatch,
		  'N' => $mismatch,
		   '.' => $gap},
	  'A' => {
		  'G' => $mismatch,
		  'g' => $mismatch,
		  'a' => $match,
		  'A' => $match,
		  'C' => $mismatch,
		  'c' => $mismatch,
		  'T' => $mismatch,
		  'n' => $mismatch,
		  'N' => $mismatch,
		  't' => $mismatch,
		  '.' => $gap},
	  'a' => {
		  'G' => $mismatch,
		  'g' => $mismatch,
		  'a' => $match,
		  'A' => $match,
		  'C' => $mismatch,
		  'c' => $mismatch,
		  'T' => $mismatch,
		  't' => $mismatch,
		  'n' => $mismatch,
		  'N' => $mismatch,
		  '.' => $gap},
	  'C' => {
		  'G' => $mismatch,
		  'g' => $mismatch,
		  'a' => $mismatch,
		  'A' => $mismatch,
		  'C' => $match,
		  'c' => $match,
		  'T' => $mismatch,
		  't' => $mismatch,
		  'n' => $mismatch,
		  'N' => $mismatch,
		  '.' => $gap},
	  'c' => {
		  'G' => $mismatch,
		  'g' => $mismatch,
		  'a' => $mismatch,
		  'A' => $mismatch,
		  'C' => $match,
		  'c' => $match,
		  'T' => $mismatch,
		  't' => $mismatch,
		  'n' => $mismatch,
		  'N' => $mismatch,
		  '.' => $gap},
	  't' => {
		  'G' => $mismatch,
		  'g' => $mismatch,
		  'a' => $mismatch,
		  'A' => $mismatch,
		  'C' => $mismatch,
		  'c' => $mismatch,
		  'T' => $match,
		  't' => $match,
		  'n' => $mismatch,
		  'N' => $mismatch,
		  '.' => $gap},
	  'T' => {
		  'G' => $mismatch,
		  'g' => $mismatch,
		  'a' => $mismatch,
		  'A' => $mismatch,
		  'C' => $mismatch,
		  'c' => $mismatch,
		  'T' => $match,
		  't' => $match,
		  'n' => $mismatch,
		  'N' => $mismatch,
		  '.' => $gap},
	  '.' => {
		  'G' => $gap,
		  'g' => $gap,
		  'a' => $gap,
		  'A' => $gap,
		  'C' => $gap,
		  'c' => $gap,
		  'T' => $gap,
		  't' => $gap,
		  'n' => $mismatch,
		  'N' => $mismatch,
		  '.' => $match},
	  'N' => {
		  'G' => $mismatch,
		  'g' => $mismatch,
		  'a' => $mismatch,
		  'A' => $mismatch,
		  'C' => $mismatch,
		  'c' => $mismatch,
		  'T' => $mismatch,
		  't' => $mismatch,
		  'n' => $mismatch,
		  'N' => $mismatch,
		  '.' => $match},
	  'n' => {
		  'G' => $mismatch,
		  'g' => $mismatch,
		  'a' => $mismatch,
		  'A' => $mismatch,
		  'C' => $mismatch,
		  'c' => $mismatch,
		  'T' => $mismatch,
		  't' => $mismatch,
		  'n' => $mismatch,
		  'N' => $mismatch,
		  '.' => $match}
);

if ($#ARGV < 2) {
  print "usage: $0 refSeqFile qryBlastTable qrySeq [-verbose]\n";
  exit 0;
}

$refSeqFile = shift @ARGV;
$qryBlastTableFile = shift @ARGV;
$qrySeqFile = shift @ARGV;
$verbose = 0;
if ($#ARGV >= 0) {
  $opt = shift @ARGV;
  if ($opt eq "-verbose") {
    $verbose = 1;
  }
}

$refSeqIn = Bio::SeqIO->new('-format'=> 'fasta',
			    -file => $refSeqFile);

print "getting sequence \n";
$refSeq = $refSeqIn->next_seq();
print "done\n";
print "reading hits\n";

open (BLAST, $qryBlastTableFile) or die "cannot open $qryBlastTableFile\n";
while(<BLAST>) {
  $line = $_;
  if ($line =~ /^#/) {
    next;
  }
  @vals = split(/\s+/, $line);
  $title = @vals[0];
  $assign = 0;
  if (! exists $topQryHits{$title} ) {
    $assign = 1;
  }
  elsif ($topQryHits{$title}{bitScore} < $vals[11]) {
    $assign = 1;
  }
  if ($assign) {
    %{$topQryHits{$title}} = ( qryId => $title,
			       refId => $vals[1],
			       pctIdentity => $vals[2],
			       alignLength => $vals[3],
			       misMatches  => $vals[4],
			       gapOpen     => $vals[5],
			       qStart      => $vals[6],
			       qEnd        => $vals[7],
			       sStart      => $vals[8],
			       sEnd        => $vals[9],
			       eVal        => $vals[10],
			       bitScore    => $vals[11] );
  }
}


$qryBlastTable =  new Bio::SearchIO(-file   => $qryBlastTableFile,
				    -format => 'blast');


# now find the similarities of the hits
$qryIn = Bio::SeqIO->new('-format'=>'fasta',
			 '-file'=> $qrySeqFile);

while ($seq = $qryIn->next_seq) {
  $title = $seq->display_id();
  if (exists $topQryHits{$title}) {
    # gather statistics about the query

    # get the locations of the hit
    $qryStart = $topQryHits{$title}{qStart};
    $qryEnd   = $topQryHits{$title}{qEnd};
    $sbjStart = $topQryHits{$title}{sStart};
    $sbjEnd   = $topQryHits{$title}{sEnd};

    $dir = 0;
    if ($sbjEnd < $sbjStart) {
      $temp = $sbjEnd;
      $sbjEnd = $sbjStart;
      $sbjStart = $temp;
      $dir = 1;
    }
    # find the locations of the endpoints
    if ($dir == 0) {
      $sbjStart = $sbjStart - $qryStart + 1;
      $sbjEnd   = $sbjEnd   + $seq->length - $qryEnd;
    }
    else {
      $sbjStart = $sbjStart - ($seq->length - $qryEnd);
      $sbjEnd   = $sbjEnd   + $qryStart - 1;
    }
    # create an alignable subsequence, and align
    if ($dir == 0) {
      $sbjSubstr = $refSeq->subseq($sbjStart, $sbjEnd);
    }
    else {
      $sbjSubstr = $refSeq->subseq($sbjStart, $sbjEnd);
      $subSeq = new Bio::Seq('-seq'=>$sbjSubstr);
      $seqRc  = $subSeq->revcom();
      $sbjSubstr = $seqRc->seq();
    }

    $sbjSubSeq = new Bio::Seq('-seq'=>$sbjSubstr);
    $seq1 = $seq->seq();
    $seq2 = $sbjSubstr;
    if ($verbose) {
      print "seq1: $seq1\n";
      print "seq2: $seq2\n";
    }
    @map = ();
    SWAlign($seq1, $seq2, \@map);

    # now get the statistics of the whole sequence (ugh, finally)
    my $numMisMatch = 0;
    my $numIns = 0;
    my $numDel = 0;
    if ($verbose) {
      print "map: @map\n";
    }
    GetAlignSummary($seq1, $seq2, \@map, \$numMisMatch, \$numIns, \$numDel);
    print "$title: $numMisMatch, $numIns, $numDel\n";
    if ($numMisMatch > 0 or
	$numIns > 0 or
	$numDel > 0 ){ 
      print "seq1: $seq1\n";
      print "seq2: $seq2\n";
      print "map:  @map\n";
    }
  }
  else {
    print "not top query found for: $title\n";
    exit(0);
  } 
}

sub SWAlign {
  my ($seq1, $seq2, $map) = @_;

  @s1 = split(//, $seq1);
  @s2 = split(//, $seq2);

  @score = ();
  @path  = ();
#  @s1 = ('.', @s1);
#  @s2 = ('.', @s2);
  $ls1 = scalar @s1;
  $ls2 = scalar @s2;
  my $matchBoth = 0;
  my $gap1 = 1;
  my $gap2 = 2;
  my $init = -1;

  $len1 = scalar @s1;
  $len2 = scalar @s2;

#  print "$seq1\n";
#  print "$seq2\n";

  my ($i, $j);

  # initialize the matrix, first row and column
  my ($matchScore, $gap1Score, $gap2Score);

  @scoreRow = ();
  $scoreRow[0] = 0;
  for ($i = 1; $i < $len1; $i++ ) {
    @scoreRow[$i] = $scoreRow[$i-1]; # + $gap; for now use fitting seq
    @pathRow[$i] = $gap1;
  }
  my @scoreMat;
  my @pathMat;
  push @scoreMat, [@scoreRow];
  push @pathMat, [@pathRow];

  # create the dp tables
  for ($i = 1; $i <= $len2; $i++ ) {
    @scoreRow = ();
    @pathRow  = ();
    push @scoreMat, ();
    push @{$scoreMat[$i]}, $scoreMat[$i-1][0] ; #+ $gap; for now use fitting seq
    push @{$pathMat[$i]}, $gap2;
    for ($j = 1; $j <= $len1; $j++) {
      $sc = $score{$s1[$j-1]}{$s2[$i-1]};
      $match1Score = $score{$s1[$j-1]}{$s2[$i-1]} + $scoreMat[$i-1][$j-1];
      $gap1Score   = $gap + $scoreMat[$i-1][$j];
      $gap2Score   = $gap + $scoreMat[$i][$j-1];
#      print "($i,$j,$match1Score,$gap1Score,$gap2Score,$s1[$j-1],$s2[$i-1],$sc) ";
      if ($match1Score < $gap1Score and
	  $match1Score < $gap2Score) {
	push @{$scoreMat[$i]}, $match1Score;
	push @{$pathMat[$i]}, $matchBoth;
      }
      elsif ($gap1Score <= $match1Score and
	     $gap1Score <= $gap2Score) {
	push @{$scoreMat[$i]}, $gap1Score;
	push @{$pathMat[$i]}, $gap1;
      }
      elsif($gap2Score <= $match1Score and
	    $gap2Score <= $gap1Score) {
	push @{$scoreMat[$i]}, $gap2Score;
	push @{$pathMat[$i]}, $gap2;
      }
      else {
	print "huh?\n";
	exit(0);
      }
    }
#    print "\n";
  }

  # now create the alignment map
#  for ($i = 0; $i <= $len2; $i++ ) {
#    for ($j = 0; $j <= $#{$pathMat[$i]}; $j++) {
#      print "$pathMat[$i][$j]";
#    }
#    print "\n";
#  }
#
#  print "align matrix: $#pathMat\n";
#  for ($i = 0; $i <= 40; $i++) {
#    printf STDOUT "%3d", $i;
#  }
#  printf "\n";
#  print "align matrix: $#pathMat\n";
#  for ($i = 0; $i <= $len2; $i++ ) {
#    for ($j = 0; $j <= 40; $j++) {
#      printf STDOUT "%3d", $scoreMat[$i][$j];
#    }
#    print "\n";
#  }

  $i = $len2;
  $j = $len1;
  @$map = ();
  for ($i = 0; $i < $len1-1; $i++ ){
    push @$map, $init;
  }
  # find the minimum start position
  $startRow = $len2-1;
  $startCol = $len1-1;
  $minRow   = $startRow;
  $minCol   = $startCol;
  $minScore = $scoreMat[$len2-1][$len1-1]; 
  for ($i = 1; $i <= $len2; $i++ ) {
    if ($minScore > $scoreMat[$i][$startCol]) {
      $minScore = $scoreMat[$i][$startCol];
      $minRow   = $scoreMat[$i][$startCol];
    }
  }
  for ($j = 1; $j <= $len1; $j++ ) {
    if ($minScore > $scoreMat[$startRow][$j]) {
	$minScore = $scoreMat[$startRow][$j];
 	$minRow = $startRow;
	$minCol = $j;
    }
  }

  $i = $minRow; 
  $j = $minCol;
  while ($i > 0 and $j > 0) {
    if ($pathMat[$i][$j] == $matchBoth) {
      $map[$j-1] = $i-1;
      $i--;
      $j--;
    }
    elsif ($pathMat[$i][$j] == $gap1) {
      $i--;
    }
    else {
      if ($pathMat[$i][$j] != $gap2) {
	print "error, path matrix is $pathMat[$i][$j] but it should be $gap2\n";
	exit 0;
      }
      $j--;
    }
  }
}

sub GetAlignSummary {
  my ($seq1, $seq2, $map,
      $numMismatch, $numIns, $numDel )  = @_;
#  print "map: @$map\n";
#  print "s1: $seq1\n";
#  print "s2: $seq2\n";
  @s1 = split(//, $seq1);
  @s2 = split(//, $seq2);

  my $i;

  $$numMismatch = 0;
  $$numIns      = 0;
  $$numDel      = 0;
  while($$map[$i] == -1) {
    $i++;
  }
  my $end = scalar @s1;
  while ($end > 0 and $map[$end-1] == -1) {
    $end--;
  }
  for ($i = 0; $i < $end; $i++ ) {
    if ($s1[$i] < 0) {
      $$numDel++;
#      print "del $i\n";
    }
    elsif ($i > 0 and $$map[$i] != ($$map[$i-1] + 1)) {
      $$numIns++;
    }
    elsif ($s1[$i] ne $s2[$$map[$i]]) {
#      print "match $i $s1[$i] $s2[$$map[$i]]\n";
      $$numMismatch++;
    }
  }
}
