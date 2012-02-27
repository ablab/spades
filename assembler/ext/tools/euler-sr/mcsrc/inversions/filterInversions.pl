#!/usr/bin/env perl

use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;

if ($#ARGV < 2) {
  print "usage: $0 invFileName seqdir dbdir \n";
  exit(1);
}

$invFileName = shift @ARGV;
$seqDir = shift @ARGV;
$dbDir  = shift @ARGV;
$scoreFile = "";

if ($#ARGV >= 0) {
  $scoreFile = shift @ARGV;
}

$binDir = $ENV{"HOME"}. "/projects/mcsrc/". $ENV{"MACHTYPE"} . "/bin";

open(INVFILE, $invFileName) or die "cannot open $invFileName\n";

# parse the output of dbinvcheck

%species = {};
$invNumber = 0;
while (<INVFILE>) {
  $line = $_;
  # look to see if this line is a header
  if ($line =~ /\// or $line =~ /\.lav/) {
    $refSpec = ""; $qrySpec = "";
    $refSeq  = ""; $qrySeq =  "";
    if ($line =~ /\//) {
      $line =~ /\.*\/((\w+)\.(\w+)\.(\w+))\.((\w+)\.(\w+)\.(\w+))\..*/;
      $refSpec = $2; $qrySpec = $6;
      $refSeq = $3; $qrySeq = $7;
      $refFile = $1; $qryFile = $5;
    }
    else {
      $line =~ /((\w+)\.(\w+)\.(\w+))\.((\w+)\.(\w+)\.(\w+))\..*/;
      $refSpec = $2; $qrySpec = $6;
      $refSeq = $3; $qrySeq = $7;
      $refFile = $1; $qryFile = $5;
    }
  }
  else {
    $line =~/(\d+) (\d+) (\d+) (\d+).*/;
    $refStart = $1; $refEnd = $2;
    $qryStart = $3; $qryEnd = $4;
    if (!exists $species{$refSpec}) {
      $species{$refSpec} = {};
    }
    if (!exists ${$species{$refSpec}}{$qrySpec}) {
      @{$species{$refSpec}{$qrySpec}} = ();
    }
    push @{$species{$refSpec}{$qrySpec}}, [$refStart, $refEnd, $qryStart, $qryEnd, $invNumber, 0, 0, $refFile, $qryFile];
  }
  ++$invNumber;
}


close INVFILE;

# Configure blast2seq
$blastOut = "blast.out";
$home = $ENV{'HOME'};
$BLASTDIR = $home . "/projects/software/blast-2.2.12/bin";
$blastSA = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn',
						 'outfile'=>$blastOut,
						 'q'=>'-2',
						 'W'=>'10'
						);

# Make a directory where all of the results will go

if (! -e "tmp") {
  `mkdir tmp`;
}

# Extract all of the inverted sequences
foreach $refKey (keys %species) {
  foreach $qryKey (keys %{$species{$refKey}}) {
    $firstInv = @{$species{$refKey}{$qryKey}}[0];
    $refInvFile = "tmp/$refKey.$qryKey.ref";
    $cmd = "$binDir/extinv $invFileName $refKey $qryKey $seqDir/@{$firstInv}[7] $refInvFile ";
    if (! -e $refInvFile) {
      `$cmd`;
      if ($? != 0) {
	print "FAILED: $cmd \n";
      }
    }
    $qryInvFile = "tmp/$refKey.$qryKey.qry";
    $cmd = "$binDir/extinv $invFileName $refKey $qryKey $seqDir/@{$firstInv}[8]  $qryInvFile -q";
    if (! -e $qryInvFile) {
      `$cmd`;
      if ($? != 0) {
	print "FAILED: $cmd \n";
      }
    }
  }
}


# Extract all of the inverted sequences
foreach $refKey (keys %species) {
  foreach $qryKey (keys %{$species{$refKey}}) {
    $firstInv = @{$species{$refKey}{$qryKey}}[0];
    $refInvFile = "tmp/$refKey.$qryKey.ref";
    $cmd = "$binDir/mskrep $refInvFile $refInvFile.msk";
    if (! -e "$refInvFile.msk") {
      `$cmd`;
      if ($? != 0) {
	print "FAILED: $cmd \n";
      }
    }
    if (! -e "$qryInvFile.msk") {
      $qryInvFile = "tmp/$refKey.$qryKey.qry";
      $cmd = "$binDir/mskrep $qryInvFile $qryInvFile.msk ";
      `$cmd`;
      if ($? != 0) {
	print "FAILED: $cmd \n";
      }
    }
  }
}

%invHSPScores = {};

# Find the scores of the HSPs for each 'reference' inversion.
@alignedInv = ();
$invNum = 0;
if ($scoreFile eq "") {
  foreach $refKey (keys %species) {
    %{$invHSPScores{$refKey}} = {};
    foreach $qryKey (keys %{$species{$refKey}}) {
      %{$invHSPcores{$refKey}{$qryKey}} = {};
      $firstInv = @{$species{$refKey}{$qryKey}}[0];
      $refSeqFile = "tmp/$refKey.$qryKey.ref.msk";
      $qrySeqFile = "tmp/$refKey.$qryKey.qry.msk";

      $refSeqio = Bio::SeqIO->new(-file => $refSeqFile,
				  -format => "fasta" );
      
      $qrySeqio = Bio::SeqIO->new(-file => $qrySeqFile,
				  -format => "fasta" );
      $pairInvNum = 0;
      while ($refSeq = $refSeqio->next_seq and 
	     $qrySeq = $qrySeqio->next_seq) {
	#      print "blasting $refSeq $qrySeq\n";
	$seqout = Bio::SeqIO->new(-file=>">tmp/ref.fasta",
				  -format=>"fasta");
	$seqout->write_seq($refSeq);
	$seqout = Bio::SeqIO->new(-file=>">tmp/qry.fasta",
				-format=>"fasta");
	$seqout->write_seq($qrySeq);
	if ($scoreFile eq "") {
	  $bl2seqReport = $blastSA->bl2seq($refSeq, $qrySeq);
	  $score = GetTotalScore(\$bl2seqReport, 0.1, 0);
	  #     print "got sum bitscore: $score\n";
	  $alignedInv[$invNum] = [$refKey, $qryKey, $pairInvNum, $score, $refSeq, $qrySeq];
	if ($scoreFile eq "") {
	  print "@{$alignedInv[$invNum]}[0] @{$alignedInv[$invNum]}[1] @{$alignedInv[$invNum]}[2] @{$alignedInv[$invNum]}[3] \n";
	}
	}
	else {
	  push @{$alignedInv[$invNum]}, ($refSeq, $qrySeq, $refSeqFile, $qrySeqFile);
	  #	print "@{$alignedInv[$invNum]} \n";
	}
	$invNum++;
	++$pairInvNum;
      }
    }
  }
  # can continue processing on a different run.
  exit(0);
}
else {
  open (SF, "$scoreFile");
  while (<SF>) {
    @vals = split /\s+/;
    $alignedInv[$invNum] = [@vals];
    $invNum++;
  }
  for $invNum (0 .. $#alignedInv) {
    $refKey = @{$alignedInv[$invNum]}[0];
    $qryKey = @{$alignedInv[$invNum]}[1];
    $refSeqFile = "tmp/$refKey.$qryKey.ref.msk";
    $qrySeqFile = "tmp/$refKey.$qryKey.qry.msk";
#    print "reading sequences $refSeqFile $qrySeqFile\n";
    $refSeqio = Bio::SeqIO->new(-file => $refSeqFile,
				-format => "fasta" );
    $qrySeqio = Bio::SeqIO->new(-file => $qrySeqFile,
				-format => "fasta" );

    while ($refSeq = $refSeqio->next_seq and 
	   $qrySeq = $qrySeqio->next_seq) {
      push @{$alignedInv[$invNum]}, ($refSeq, $qrySeq, $refSeqFile, $qrySeqFile);
    }
  }
}


@bins = ();
@binSpecies = ();
@invBins = ();

for $invNum (0 .. $#alignedInv) {
  $binFound = 0;
  print "searching for match for @{$alignedInv[$invNum]}[0]  @{$alignedInv[$invNum]}[1]  @{$alignedInv[$invNum]}[2]\n";
  for $binNum (0 .. $#bins) {
    # try to find a match of this sequence to the bin
    print "  comparing sequences in bin: $binNum\n";
    for $binSeq ( 0 .. $#{$bins[$binNum]}) {
      $seqIndex = @{$bins[$binNum]}[$binSeq];
      # try aligning ref seq to binned qry seq
      $refSeq = \@{$alignedInv[$invNum]}[4];
      $qrySeq = \@{$alignedInv[$seqIndex]}[4];
      $seqout = Bio::SeqIO->new(-file=>">tmp/ref.fasta",
				-format=>"fasta");
      $seqout->write_seq($$refSeq);
      $seqout = Bio::SeqIO->new(-file=>">tmp/qry.fasta",
				-format=>"fasta");
      $seqout->write_seq($$qrySeq);
      $bl2seqReport = $blastSA->bl2seq(@{$alignedInv[$invNum]}[4],
				       @{$alignedInv[$seqIndex]}[4]);
      $refLen = ${@{$alignedInv[$invNum]}}[4]->length();
      $qryLen = ${@{$alignedInv[$seqIndex]}}[4]->length();
      $refSeqName = ${@{$alignedInv[$invNum]}}[4]->display_id();
      $qrySeqName = ${@{$alignedInv[$seqNum]}}[4]->display_id();
      $score = GetTotalScore(\$bl2seqReport, 0.1, 0);
      print "  ref bininv: $seqIndex, ref vs bin: $score, binscore: @{$alignedInv[$seqIndex]}[3] rl: $refLen ql: $qryLen $refSeqName $qrySeqName.\n";
#      exit(0);
       if (SameOrder($score, @{$alignedInv[$seqIndex]}[3], 5)) {
	print "    found a bin\n";
	$binFound = 1;
	last;
      }
      # try aligning ref seq to binned qry seq
      $bl2seqReport = $blastSA->bl2seq(@{$alignedInv[$invNum]}[4],
				       @{$alignedInv[$seqIndex]}[5]);
      $score = GetTotalScore(\$bl2seqReport, 0.1, 0);
      $refLen = ${@{$alignedInv[$invNum]}}[4]->length();
      $qryLen = ${@{$alignedInv[$seqIndex]}}[5]->length();
      $refFile = ${@{$alignedInv[$invNum]}}[4]->display_id();
      $qryFile = ${@{$alignedInv[$seqNum]}}[5]->display_id();
      print "  qry bininv: $seqIndex, ref vs bin: $score, binscore: @{$alignedInv[$seqIndex]}[3] rl: $refLen ql: $qryLen $refSeqName $qrySeqName.\n";
      if (SameOrder($score, @{$alignedInv[$seqIndex]}[3], 5)) {
	print "  found a bin\n";
	$binFound = 1;
	last;
      }
    }
    if ($binFound == 1) {
      print "  adding $invNum to bin: $binNum\n";
      push @{$bins[$binNum]}, $invNum;
      $binIndexNum = $binNum;
      last;
    }
  }
  if ($binFound == 0) {
    # no bin was found forthis, create a new bin
    print "creating a bin @{$alignedInv[$invNum]}[0] @{$alignedInv[$invNum]}[1] @{$alignedInv[$invNum]}[3]  \n";
    push @bins, [$invNum];
    $binIndexNum = $#bins;
  }
  print "assigning inversion $invNum @{$alignedInv[$invNum]}[0] @{$alignedInv[$invNum]}[1] to bin: $binIndexNum\n";
  push @invBins, $binIndexNum;
}

sub SameOrder {
  my ($val1, $val2, $fact) = @_;
  if ($val1 > $val2) { 
    $p = $fact * $val2;
#  print "so $val1, $val2, $p, $fact\n";
    if ( $p > $val1) {
      return 1;
    }
    else {
      return 0;
    }
  }
  else {
    $p = $fact * $val1;
#  print "so $val1, $val2, $p, $fact\n";
    if ( $p > $val2) {
      return 1;
    }
    else {
      return 0;
    }
  }
}

sub GetTotalScore {
  my ($blastReport, $eThresh, $orientation) = @_;
  $totalScore = 0;
  $numHsps = 0;
  while ( $result = $$blastReport->next_result ) {
    while ($hit = $result->next_hit ) {
      while ($hsp = $hit->next_hsp ) {
	($rBegin, $rEnd) = $hsp->range('query');
	($qBegin, $qEnd) = $hsp->range('sbjct');
	$score = $hsp->score();
	$evalue = $hsp->evalue();
	if ($evalue < $eThresh) {
	  print "    got range: $rBegin, $rEnd, $qBegin, $qEnd $score $evalue\n";
	  $totalScore += $score;
	}
	$numHsps++;
      }
    }
  }
#  print "got $numHsps hsps\n";
  return $totalScore;
}
