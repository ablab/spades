#!/usr/bin/env perl

use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use common;
use DBI;

if ($#ARGV < 6) {
  print "usage: checkBin.pl dbname binSeqFile speciesListFile regionName  seqDir dbDir orthFileName confirmedseqName eThresh\n";
  exit(0);
}

$sqlDb           = shift @ARGV;
$binSeqFileName  = shift @ARGV;
$seqSetFileName  = shift @ARGV;
$regionName      = shift @ARGV;
$seqDir          = shift @ARGV;
$dbDir           = shift @ARGV;
$orthFileName    = shift @ARGV;
$confirmedSeqName= shift @ARGV;
$eThresh         = shift @ARGV;

# non required options
$dataFile = "";
while ($#ARGV >= 0) {
  $opt = shift @ARGV;
  if ($opt eq "-d") {
    $dataFile = shift @ARGV;
  }
}

$socket = "";

$dbConnect = "dbi:mysql:$sqlDb";
if (exists $ENV{"GSCK"}) {
  $socket = $ENV{"GSCK"};
  $dbConnect .= ";mysql_socket=$socket";
  print "connect now: $dbConnect\n";
}

print "starting with dbName                 = $sqlDb\n";
print "starting with option binSeqFileName  = $binSeqFileName\n";
print "starting with option seqSetFileName  = $seqSetFileName\n";
print "starting with option regionName      = $regionName\n";
print "starting with option seqDir          = $seqDir\n";
print "starting with option dbDir           = $dbDir\n";
print "starting with option orthFileName    = $orthFileName\n";
print "starting with option confirmedSeqName = $confirmedSeqName\n";
print "starting with option eThresh         = $eThresh\n";

if (-e $confirmedSeqName) {
  print "skipping $confirmedSeqName\n";
  exit(0);
}

%data = ();
print "df: $dataFile\n";
if ($dataFile ne "") {
  print "reading data fiel \n";
  common::ReadDataFile($dataFile, \%data);
}

$dbh = DBI->connect($dbConnect, "mchaisso", "");

$invSeqReader = Bio::SeqIO->new(-file => $binSeqFileName,
				-format => "fasta" );


# Read the binned fasta files.
@sequences = ();
while ($invSeq = $invSeqReader->next_seq) {
  push @sequences, $invSeq;
  $did = $invSeq->display_id();
  $len = $invSeq->length();
  print "got seq: $did $len\n";
}


# read what specie are present.
open (SPECIES_LIST, $seqSetFileName) or die "cannot open $seqSetFileName \n";
$line = <SPECIES_LIST>;
@species = split(/\s+/, $line);
close SPECIES_LIST;
%species = ();
foreach $s (@species) {
  $species{$s} = 0;
}


# read the coordinates of every inversion in human
%orthoPos = ();
common::ReadOrthData($orthFileName, \%orthoPos);


# find out what species are represented in the bin
for ($s = 0; $s < scalar @sequences; $s++ ) {
  $key = $sequences[$s]->display_id();
  $key =~ /(\w+)\..*/;
  $spec = $1;
  $species{$1} = 1;
}


# find out what species are NOT represeneted in the bin.
@without = ();
@with = ();
foreach $s (keys %species) {
  if ($species{$s} == 0) {
    push @without, $s;
  }
  else {
    push @with, $s;
  }
}

@without = sort @without;

print "species not represented are: @without \n";
print "species represented: @with \n";
@newSeq = ();
$sequencesAdded = 1;
$newIndex = 0;
$oldLast = 0;
while ($sequencesAdded == 1) {
  $sequencesAdded = 0;
  $w = 0;
  $oldLast = $#sequences + 1;

  while ($w <= $#without) {
    # try to search each inversion against the database where without is
    print "trying to fit @without[$w] of $#without starting search at $newIndex\n";
    # no chicken
    if (@without[$w] eq "chicken") {
      splice @without, $w, 1;
      next;
    }
    $qryDBName = @without[$w] . ".fasta";
    my $qrySpecies = @without[$w];

    foreach $withIndex ($newIndex .. $#sequences) {
      my $refSeqName = @sequences[$withIndex]->display_id();
      $refSeqName =~ /(\w+)\..*/;
      $refSpecies = $1;
      print "with $refSeqName $refSpecies\n";
      my ($refHumanStartPos, $refHumanEndPos);
      $refHumanStartPos = ${$orthoPos{$refSpecies}}[2];
      $refHumanEndPos = ${$orthoPos{$refSpecies}}[3];
      my $foundForward = 0;
      my $foundReverse = 0;
      my $humanStartPos;
      my $humanEndPos;
      my @hsps = ();
      my $qrySubSeqStart;
      my $qrySubSeqEnd;
      ($evalue) = FindInDB(@sequences[$withIndex], $qryDBName,
			   $regionName, $dbDir, $eThresh,
			   \@hsps);
      if ($#hsps >= 0) {
	# found several blast hits.  We don't want to take this one if 
	# it maps to both orientations at the same locus.
	my $matchBoth = 0;
	# Find the coordinates of the query sequence.
	my $refLength = @sequences[$withIndex]->length;
	print "found $#hsps hsps for $refSeqName\n";
	my $maxhsps = common::min($#hsps, 10);
	for ($h = 0; $h <= $maxhsps; $h++ ) {
	  print "  testing hsp: $h\n";
	  # Store a bunch of information from the blast hit.
	  my $refStart  = $hsps[$h]->start('query');
	  my $refEnd    = $hsps[$h]->end('query');
	  my $qryStart  = $hsps[$h]->start('sbjct');
	  my $qryEnd    = $hsps[$h]->end('sbjct');
	  my $qrySeq    = $hsps[$h]->seq('sbjct');
	  my $qryName   = $hsps[$h]->seq('sbjct')->display_id();
	  my $qryStrand = $hsps[$h]->strand('sbjct');
	  my $refStrand = $hsps[$h]->strand('query');
	
	  # Given the position of the HSP find the actual endpoints
	  # of the inversion in the query strand.  Since the HSP
	  # represents a local alignment, use length offsets to calculate
	  # the endpoints.
	  if ($qryStrand == 1) {
	    $qrySubSeqStart = $qryStart - $refStart;
	    $qrySubSeqEnd   = $qryEnd + ($refLength - $refEnd);
	  }
	  else {
	    print "   assuming reverse match\n";
	    $qrySubSeqStart = $qryStart - ($refLength - $refEnd);
	    $qrySubSeqEnd   = $qryEnd + $refStart;
	  }
	
	  print "   got refStrand: $refStrand qryStrand: $qryStrand pos $qryStart -> $qrySubSeqStart $qryEnd -> $qrySubSeqEnd\n";
	  # find the corresponding location in human
	  $humanStartPos = common::GetOrthoPosOnly(@without[$w],
						   "human",
						   $regionName,
						   $qrySubSeqStart, 0, $sqlDb, 1);

	  $humanEndPos = common::GetOrthoPosOnly(@without[$w],
						 "human",
						 $regionName,
						 $qrySubSeqEnd, 0, $sqlDb, 1);
	
	  print "   human orthologous position for $qryStart: $humanStartPos \n";
	  print "   human orthologous position for $qryEnd: $humanEndPos \n";
	  # look to see if the positions map to the same location
	  my ($startDiff, $endDiff);
	  # Coordinates of inversion in reference sequence.
	  my $refStartPos;
	  my $refEndPos;
	  # human orthologous position of ref
	
	  $startDiff = abs($refHumanStartPos - $humanStartPos);
	  $endDiff = abs($refHumanEndPos - $humanEndPos);
	  print "   comparing $seqName \n\t$refHumanStartPos - $humanStartPos = $startDiff " .
	    " and\n\t$refHumanEndPos - $humanEndPos = $endDiff\n";
	  $overlap = abs(common::ComputeOverlap($refHumanStartPos, $refHumanEndPos,
						$humanStartPos, $humanEndPos));
	  # use abs here since sometimes the reference positions are out of order.
	  my $length = abs($refEndPos - $refStartPos);
	  print "   overlap: $overlap, length: $length\n";
	  #	if ($overlap > 0.5 * $length) {
	  if (($startDiff < 2000 and $endDiff < 2000) or $overlap > 0.5*$length) {
	    # found a hit that is close enough to the inversion to be considered 
	    # a match.
	    if ($qryStrand == -1) {
	      print "  found reverse match\n";
	      $foundReverse = 1;
	    }
	    if ($qryStrand == 1) {
	      print "  found forward match\n";
	      $foundForward = 1;
	    }
	  }
	} # end looping over blast hits.
	
	# Found 0 or more blast hits.  Hopefully only one direction has been 
	# found.  If not, the inversion is actually palindromic, or an inverted
	# duplication, and it should be discarded.  Also, the species should not
	# be considered again for searching.
	if ($foundReverse == 1 and $foundForward == 1) {
	  # found a match to both strands.
	  print "ambiguous match, removing element at $w of $#without\n";
	  splice @without, $w, 1;
	  $w--;
	  last;
	}
	elsif ($foundReverse == 1 or $foundForward == 1) {
	  # found a match to only one strand
	  # consider the first (highest scoring) match only.
	  print "matched only one strand: reverse($foundReverse) forward($foundForward)\n";
	  $hspRef = $hsps[0];
	  my $refStart  = $hspRef->start('query');
	  my $refEnd    = $hspRef->end('query');
	  my $qryStart  = $hspRef->start('sbjct');
	  my $qryEnd    = $hspRef->end('sbjct');
	  my $qrySeq    = $hspRef->seq('sbjct');
	  my $qryName   = $hspRef->seq('sbjct')->display_id();
	  my $qryStrand = $hspRef->strand('sbjct');
	  my $refStrand = $hspRef->strand('query');
	  if ($qryStrand == 1) {
	    $qrySubSeqStart = $qryStart - $refStart;
	    $qrySubSeqEnd   = $qryEnd + ($refLength - $refEnd);
	  }
	  else {
	    print "   assuming reverse match\n";
	    $qrySubSeqStart = $qryStart - ($refLength - $refEnd);
	    $qrySubSeqEnd   = $qryEnd + $refStart;
	  }
	  print "got stats for inversion: $qryStart $qryEnd $qryStrand\n";
	  $orthoPos{@without[$w]} = [$qrySubSeqStart, $qrySubSeqEnd, $humanStartPos, $humanEndPos];
	  print "storing $humanStartPos, $humanEndPos for $qryName\n";
	  $sequencesAdded  = 1;
	  # Found the inverted sequence, store it.
	  $hspRef->seq('sbjct')->display_id(@without[$w]);
	  my $matchStrand = $hspRef->strand();
	  $evalue = $hspRef->evalue();
	  print "found inverted sequence $qrySpecies $qryStart $qryEnd $evalue\n";
	  # Now extract the query sequence.
	  #
	
	  my $qrySeqName = "$seqDir/$qrySpecies.fasta";
	  my $fullSeqIO = Bio::SeqIO->new(-file => $qrySeqName,
					  -format => "fasta" );
	  my $qryFullSeq = $fullSeqIO->next_seq;
	  if ($qrySubSeqStart < 0) {
	    $qrySubSeqStart = 0;
	  }
	  if ($qrySubSeqEnd >= $qryFullSeq->length()) {
	    $qrySubSeqEnd = $qryFullSeq->length()-1;
	  }
	  my $qrySubSeqLen = $qrySubSeqEnd - $qrySubSeqStart;
	  print "extracting $qrySubSeqStart $qrySubSeqEnd $qrySubSeqLen\n";

	  # Mask the repetitive sequence in the match so that no new
	  # spurrious matches are found when trekking.
	  my $qrySubStringSeq = $qryFullSeq->subseq($qrySubSeqStart, $qrySubSeqEnd);
	  $qrySubStringSeq =~ tr/actg/N/;
	  my $qrySubSeq = Bio::Seq->new(-seq=>$qrySubStringSeq,
				      -display_id=>$qryName);
	  $qrySubSeq->display_id("@without[$w].$qrySubSeqStart.$qrySubSeqEnd");
	  if ($qryStrand == -1) {
	    $qrySubSeqRC = $qrySubSeq->revcom;
	    $qrySubSeq = $qrySubSeqRC;
	  }
	  push @sequences, $qrySubSeq;
	  # don't bother checking this species again.
	  print "found match, removing element at $w of $#without\n";
	  splice @without, $w, 1;
	  $w--;
	  print "\n\n\n";
	  last;
	}
      }
      print "\n";
    }
    $w++;
  }
  $newIndex = $oldLast;
}

# Now print a list of orientations of inversions.
# for each sequence known to be with
$binWriter = Bio::SeqIO->new(-file=>">$confirmedSeqName",
			     -format=>"fasta");
foreach $seqIndex (0 .. $#sequences ) {
  $binWriter->write_seq(@sequences[$seqIndex]);
}


sub FindInDB {
  my ($sequencePtr, $qrySpecies, $seqName, $dbDir, $eThresh, $hsps) = @_;
  my $dbName = $dbDir . "/". $qrySpecies;
  my $seqName = $sequencePtr->display_id();
  my $blastSeq = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn',
						       'database'=>$dbName,
						       'q'=>'-1',
						       'W'=>'7',
						       'e'=>"$eThresh",
						      );
  my $seqPrinter = Bio::SeqIO->new(-format=>"fasta");

  my $tmpOut = Bio::SeqIO->new(-file=>">tmpcf.fasta");

  $blastRes  = $blastSeq->blastall($sequencePtr);
  $minEValue = GetHSPS(\$blastRes, $hsps);
  if ($minEValue < 0 or
      $minEValue > $eThresh or
      $#$hsps < 0) {
    print "no match \n";
    return (-1);
  }
  else {
    print "match with evalue: $minEValue \n";
    return ($minEValue);
  }
}

sub FindWithout {
  my ($bin, $speciesNames, $invList, $with, $without) = @_;

  # get which names are present from the bins
  my @invNames = ();
  my @binKeys  = keys %{$bin};
  my $nbk = $#binKeys + 1;
  foreach my $invIndex (keys %{$bin}) {
    my @invDes = @{${$invList[$invIndex]}[0]};
    print "adding @invDes[0]\n";
    push @invNames, @invDes[0];
    print "adding @invDes[3]\n";
    push @invNames, @invDes[3];
  }

  my $invName;
  my %withSet = ();
  my %withoutSet = ();
  my ($refName, $qryName, $start, $end);
  #  print "got invNames: @invNames\n";
  # generate a list of unique names with the inversion.
  foreach my $invName (@invNames) {
    $withSet{$invName} = 1;
  }
  @{$with} = keys %withSet;
  my $specName;
  foreach $specName (@{$speciesNames}) {
    if (! exists $withSet{$specName} ) {
      push @{$without}, $specName;
    }
  }
}

sub GetHSPS {
  my ($blastResPtr, $hsps) = @_;
  $minEValue = -1;
  while ( $result = $$blastResPtr->next_result ) {
    while ($hit = $result->next_hit ) {
      while ($hsp = $hit->next_hsp ) {
	$v = $hsp->start('sbjct');
	if ($minEValue == -1 or $minEValue > $hsp->evalue()) {
	  $minEValue = $hsp->evalue();
	}
	push @{$hsps}, $hsp;
      }
    }
  }
  return $minEValue;
}

