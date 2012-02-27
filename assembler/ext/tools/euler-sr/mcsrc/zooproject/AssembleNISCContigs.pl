#!/usr/bin/env perl
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use zooproject;

if (@ARGV < 1) {
  print "usage: assemble_contigs.pl targettable target nametable overlapname [speciesid] [-noalign]\n";
  exit(0);
}

$targetTableName = shift @ARGV;
$target = shift @ARGV;
$nameTableName = shift @ARGV;
$ovpName = shift @ARGV;
#set up default options
$speciesId = -1;
$align = 1;
#read in options
while ($#ARGV >= 0) {
  $option = shift @ARGV;
  if ($option eq "-noalign") {
    $align = 0;
  }
  else {
    $speciesId = $option;
  }
}

print "$align $speciesId $ovpName\n";
my @targetTable;
zooproject::ParseTargetTable($targetTableName, \@targetTable);


# Read in name ids
my %niscToName;
my %niscToId;
my %nameToId;
my @idToName;
zooproject::ParseNameTable($nameTableName,
			   \%niscToName, \%niscToId, \%nameToId, \@idToName);

$bl2seqOut = "blastout".$$;
$blastSA = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn',
						 'outfile'=>"blah",
						 'm'=>'T');


open (OVP, ">$ovpName") or die "cannot open $ovpName\n";

# assemble each species
while ($#targetTable >= 0) {
  while ($#targetTable >=0 and $targetTable[0]{"Status"} ne "Sequenced") {
    shift @targetTable;
  }
  $clone1Name = $targetTable[0]{"Organism"};
  $clone1Acc = $targetTable[0]{"GenbankId"};
  $cloneId = $niscToId{$clone1Name};
  if ($speciesId < 0 or $cloneId == $speciesId) {
    system("mkdir -p $target.$cloneId");
    shift @targetTable;
    while ($#targetTable >= 0 and
	   $targetTable[0]{"Organism"} eq $clone1Name) {
      if ($targetTable[0]{"Status"} ne "Sequenced") {
        shift  @targetTable;
	next;
      }
      $clone2Acc = $targetTable[0]{"GenbankId"};

      # prepare the names of the sequences to align
      $seq1Name = "$target.$cloneId.$clone1Acc.fasta";
      $seq2Name = "$target.$cloneId.$clone2Acc.fasta";
      $alignFile = "$target.$cloneId/$clone1Acc.$clone2Acc.lav";

      # only align if necessary (specified on command line)
      if ($align == 1) {
	system("blastz $seq1Name $seq2Name B=0 M=2  T=0 W=12 > $alignFile ");
      }


      # Find the highest scoring alignment.
#      print "opening $alignFile\n";
      open(AF, $alignFile);
      @alignment = <AF>;
      @alignScores = ();
      ($refLen, $qryLen) = GetLengths(\@alignment);
      StoreAlignmentScores(\@alignment, \@alignScores, 100000, 10000);
      ($maxQryScore, $maxQryIndex, $maxQryPos) =
	 FindBestFrontOverlap(\@alignScores, $refLen, $qryLen);
      # initialize to default of no-overlap
      $nalign = scalar @alignScores;
      print "got $nalign align scores\n";
      $refStart = -1; $refEnd = -1; $qryStart = -1; $qryEnd = -1;
      if ($maxQryIndex >= 0) {
	$refStart = $alignScores[$maxQryIndex][1];
	$refEnd   = $alignScores[$maxQryIndex][2];
	$qryStart = $alignScores[$maxQryIndex][3];
	$qryEnd   = $alignScores[$maxQryIndex][4];
	print "$maxQryScore lengths: $refLen $qryLen maxalignment: $refStart $refEnd $qryStart $qryEnd\n";

	# determine where the alignment specifies the sequences ovelap
	# qryStart should be 1
	$refAlignStart = $refStart;
	$refAlignEnd   = $refLen;
	$refOvpLen = $refLen - $refAlignStart;
	$qryOvpLen = $refOvpLen;

	$refStart = $refAlignStart;
	$qryStart = 1;
	$refEnd   = $refAlignEnd;
	$qryEnd   = $qryOvpLen;
      }

      print OVP "$target.$cloneId.$clone1Acc.fasta $refStart $refEnd $target.$cloneId.$clone2Acc.fasta $qryStart $qryEnd\n";
      print "$target.$cloneId.$clone1Acc.fasta $refStart $refEnd $target.$cloneId.$clone2Acc.fasta $qryStart $qryEnd\n";

      $clone1Acc = $clone2Acc;
      print "\n\n\n";
      shift @targetTable;
    }
  }
  else {
    shift @targetTable;
  }
}
close OVP;



sub GetLengths {
  my ($alignment) = @_;
  my $rl = -1;
  my $ql = -1;
  my $i;
  while ($i <=  $#$alignment and
	$rl == -1) {
    if ($$alignment[$i] =~ /^s {/) {
      $$alignment[$i+1] =~ /\".*\" (\d+) (\d+) (\d+) (\d+)/;
      $rl = $2 - $1;
      $$alignment[$i+2] =~ /\".*\" (\d+) (\d+) (\d+) (\d+)/;
      $ql = $2 - $1;
    }
    $i++;
  }
  return ($rl, $ql);
}


sub StoreAlignmentScores {
  my ($alignment, $alignscores, $minalignscore, $maxstartpos) = @_;
  my $i;
  $i = 0;
  while ($i <= $#$alignment) {
    if ($$alignment[$i] =~ /^  s (\d+)/) {
      $score = $1;
      $i++;
      $$alignment[$i] =~ /b (\d+) (\d+)/;
      $br = $1; $bq = $2;
      $i++;
      $$alignment[$i] =~ /e (\d+) (\d+)/;
      $er = $1; $eq = $2;
#      print "$score $br $bq $er $eq ";
      if ($score > $minalignscore and $bq < $maxstartpos) {
	push @$alignscores, [$score, $br, $er, $bq, $eq];
#	print "ok\n";
      }
      else {
#	print "no\n";
      }

    }
    $i++;
  }
}

sub FindBestFrontOverlap {
  # find the highest scoring alignment in the first half
  my ($alignment, $refLength, $qryLength) = @_;
  my $i;
  my $refMiddle = $refLength / 2;
  my $qryMiddle = $qryLength / 2;
  my $maxRefScore = -9999999;
  my $maxRefIndex = -1;
  my $maxRefPos = -1;
  my $maxQryScore = -9999999;
  my $maxQryIndex = -1;
  my $maxQryPos = -1;
  $i = 0;
  $nalign = $#$alignment;
  while ($i <= $#$alignment) {
    if ($$alignment[$i][1] < $refMiddle and
	$$alignment[$i][0] > $maxRefScore) {
      $maxRefScore = $$alignment[$i][0];
      $maxRefIndex = $i;
      $maxRefPos = $$alignment[$i][1];
    }
    if ($$alignment[$i][3] < $qryMiddle and
	$$alignment[$i][0] > $maxQryScore) {
      $maxQryScore = $$alignment[$i][0];
      $maxQryIndex = $i;
      $maxQryPos = $$alignment[$i][3];
    }
    $i++;
  }
  return ($maxQryScore, $maxQryIndex, $maxQryPos);
}


sub FindSequenceOverlap {
  my ($seq1Name, $seq2Name) = @_;

  $seqio = Bio::SeqIO->new(-file => $seq1Name);
  $seq1  = $seqio->next_seq;

  $seqio = Bio::SeqIO->new(-file => $seq2Name);
  $seq2  = $seqio->next_seq;
  $qryLen = 10000;
  if ($seq2->length < 10000) {
    $qryLen = $seq2->length;
  }
  $seq2qrystr = $seq2->subseq(1,$qryLen-1);
  $seq2qry = Bio::Seq->new(-seq=>$seq2qrystr);
  $bl2seqReport = $blastSA->bl2seq($seq2, $seq1);

  $refOvpStart = -1;
  $refOvpEnd   = -1;
  $qryOvpStart = -1;
  $qryOvpStart = -1;
  # get the highest scoring hit
  if ( $result = $bl2seqReport->next_result ) {
    if ($hit = $result->next_hit ) {
#      print "got a hit\n";
      if ($hsp = $hit->next_hsp ) {
	($rBegin, $rEnd) = $hsp->range('sbjct');
	($qBegin, $qEnd) = $hsp->range('query');
	$rlen = $seq1->length;
#	print "got blast result: $rBegin, $rEnd ($rlen) $qBegin, $qEnd\n";
	# determine roughly the overlapping side
      }
    }
  }
  return ($refOvpStart, $refOvpEnd, $qryOvpStart, $qryOvpEnd);
}


