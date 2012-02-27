#!/usr/bin/env perl

use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use common;
use DBI;

if ($#ARGV < 5) {
  print "usage: $0 invFileName seqdir dbFile seqname orthoFile binOut [dbname] [socket]\n";
  print "example: $0 inv/ENm001.inversions inv ENm001.all.fasta ENm001 ortho.txt/\n";
  exit(1);
}

$invFileName = shift @ARGV;
$seqDir      = shift @ARGV;
$dbFile      = shift @ARGV;
$regionName  = shift @ARGV;
$orthoFile   = shift @ARGV;
$binFileOut  = shift @ARGV;
$socket="";
$dbName = "zoo";
if ($#ARGV >= 0) {
  $dbName = shift @ARGV;
}
if ($#ARGV >= 0) {
  $socket      = shift @ARGV;
}
if ($socket ne "") {
  $dbConnect = "dbi:mysql:$dbName;mysql_socket=$socket";
}
else {
  $dbConnect = "dbi:mysql:$dbName";
}

print "starting with invFileName: $invFileName\n";
print "starting with seqDir      $seqDir  \n";
print "starting with dbFile      $dbFile      \n";
print "starting with regionName  $regionName  \n";
print "starting with orthoFile   $orthoFile   \n";
print "starting with binFileOut  $binFileOut  \n";
print "starting with socket      $socket\n";

$dbh = DBI->connect($dbConnect, "mchaisso", "");

$binDir = $ENV{"HOME"}. "/projects/mcsrc/". $ENV{"MACHTYPE"} . "/bin";

%pairs = ();
@refSpecies = ();

common::ReadInvFile($invFileName, \@refSpecies, \%pairs);
print "got ref species @refSpecies \n";
# Configure blast2seq
$blastOut = "blast.out";
$home = $ENV{'HOME'};

$blastall = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn',
						  'database'=>$dbFile,
						  'outfile'=>$blastOut,
						  'q'=>'-1',
						  'W'=>'8',
						  'e'=>'0.001',
						 );

# Make a directory where all of the results will go

if (! -e "tmp") {
  `mkdir tmp`;
}

# Read in all of the sequences corresponding to the inversions.
@invList = ();
$invIndex = 0;

# invSet is a map from title to sequence.
%nameToSequence = ();
%nameToIndex  = ();

common::ReadFiles(\@refSpecies, \%pairs, \@invList, \%nameToSequence, \%nameToIndex, ".");

# create a list of orthologous positions to speed things up

%orthoPos = ();
common::ReadOrthoCoordinates($orthoFile, \%orthoPos);

# each bin contians a list of sequenes that align and 
# map to a similar location on the human sequence.

%bins = ();

# Create a list of indices of inversions
%nameToBin = ();
foreach $invName (keys %nameToSequence) {
  $nameToBin{$invName} = -1;
}

# find out which inversions correspond to the same sequence
$binIndex = 0;
$invSeqIndex = 0;
print "binning sequences \n";
%binnedSequences = ();

# find the entry I'm interested in

#$test = "Chinese_muntjac_758863_760540_Baboon_olive"; #"Armadillo_nine_banded_Lemur_gray_mouse_486442_487080";
#$invSeqIndex = $nameToIndex{$test};

#print "got invseqindex: $invSeqIndex\n";

#for ($invSeqIndex = 0; $invSeqIndex <= $#invList; $invSeqIndex++) {
#  if ($invList[invSeqIndex][0][0][0] eq ) {
#    last;
#  }
#}
print "got $#invList elements in invlist \n";
while ($invSeqIndex <= $#invList) {
 # blast every inversion to the db.
  $refSeq_ptr = \@{$invList[$invSeqIndex]}[-1];
  $qrySeq_ptr = \@{$invList[$invSeqIndex+1]}[-1];
  # fish with the reference sequence
  $refName = @{$invList[$invSeqIndex]}[-3];
  $qryName = @{$invList[$invSeqIndex]}[-2];
  $refSeqStr = $$refSeq_ptr->seq();
  $refSeqLen = $$refSeq_ptr->length();
  $qrySeqStr = $$qrySeq_ptr->seq();
  $qrySeqLen = $$qrySeq_ptr->length();

  my ($toStartPos, $toEndPos);
  $toStartPos = @{$orthoPos{$invList[$invSeqIndex][-3]}}[0];
  $toEndPos   = @{$orthoPos{$invList[$invSeqIndex][-3]}}[1];
  print "looking up name: $invList[$invSeqIndex][-3]  start: $toStartPos end: $toEndPos $refSeqLen\n";
  print "translating refseqstr: $refSeqStr\n";
  my $refNumMasked  = ($refSeqStr =~ tr/Nactg//);
  my $qryNumMasked  = ($qrySeqStr =~ tr/Nactg//);
  
  my ($refMaskedRatio, $qryMaskedRatio);
  if ($refSeqLen > 0) {
    $refMaskedRatio = $refNumMasked / $refSeqLen;
  }
  else {
    $refMaskedRatio = -1;
  }
  if ($qrySeqLen > 0) {
    $qryMaskedRatio = $qryNumUnmasked / $qrySeqLen;
  }
  else {
    $qryMaskedRatio = -1;
  }
  print "unmasked $refMaskedRatio $qryMaskedRatio\n";
  if (! exists $binnedSequences{$refName} and
      $refSeqLen > 0 and 
      $refMaskedRatio < 0.90) {
    my $out = Bio::SeqIO->new('-file'=>">tmp.fasta");
    $out->write_seq($$refSeq_ptr);
    $blastReport = $blastall->blastall($$refSeq_ptr);
    %matches = ();
    my ($startHumanPosition, $endHumanPosition);
    my ($fromSeq, $toSeq, $toStrand);

    GetMatchList(\$blastReport, 0.001,
		 $toStartPos, $toEndPos,
		 \%orthoPos, $regionName, \%matches,
		 \@invList, \%nameToIndex);

    $binUpdatedRef = BinMatchList(\%matches, \%nameToIndex,
				  \@bins, \%nameToBin,
				  $binIndex, \%binnedSequences,
				  \%orthoPos, $toStartPos, $toEndPos,
				  \@invList);
    print "after binning match list: $#bins bins exist \n";
  }

  # fish with the query sequence

  if (! exists $binnedSequences{$qryName} and
       $qrySeqLen > 0 and 
      (1.0*$qryNumMasked ) / $qrySeqLen > 0.25) {
    $blastReport = $blastall->blastall($$qrySeq_ptr);
    $qryInvStart = @{$orthoPos{$invList[$invSeqIndex+1][-2]}}[0];
    $qryInvEnd   = @{$orthoPos{$invList[$invSeqIndex+1][-2]}}[1];

    # print "mapped from $invList[$invSeqIndex][0][3] ($qryInvStartFor)
    # to $toStartPos and $qryInvEndFor $toEndPos\n";

    GetMatchList(\$blastReport, 0.005,
		 $toStartPos, $toEndPos,
		 \%orthoPos, $regionName, \%matches, \@invList, \%nameToIndex);

    $binUpdatedQuery = BinMatchList(\%matches,
				    \%nameToIndex,
				    \@bins,
				    \%nameToBin,
				    $binIndex,
				    \%binnedSequences,
				    \%orthoPos, $toStartPos, $toEndPos,
				    \@invList);

    print "after binning match list: $#bins bins exist \n";
    ++$binIndex;
  }
  # debugging information
  # $refId = $$refSeq_ptr->display_id();

  # matches are returned as the names of the sequences hit
  # try to link this cluster with another cluster.  Each bin is a set
  # of svequence names.  There is a match any time one sequence is found
  # in a bin.

  print "--------------------------------------------------------------\n";
  $invSeqIndex+=2;
  # determine what species are missing from each bin, and try
  # to locate them using blast.
}

# print the bins to a file
foreach $binNumber (0 .. $#bins ) {
  $binSequenceName = "$binFileOut.$binNumber.fasta";
  $binDataName     = "$binFileOut.$binNumber.data";

  print "writing $binNumber: $binSequenceName, $binDataName\n";

  open (BINDATA, ">$binDataName") or die "cannot write to $binDataName";
  $binWriter = Bio::SeqIO->new(-file=>">$binSequenceName",
			       -format=>"fasta");
  foreach $seqIndex (keys %{$bins[$binNumber]}) {
    print BINDATA "@{$invList[$seqIndex][0]} $invList[$seqIndex][1] $invList[$seqIndex][2]\n";
    $binWriter->write_seq($invList[$seqIndex][-1]);
  }
  close BINDATA;
}


# done
exit(0);


################################################################################
# functions
################################################################################

sub GetMatchList {
  # Given an inverted sequence, find all other inversions that have a similar
  # sequence, and map to a similar location on the human genome.
  # Each inversion has a sequence in a reference and query species.  This considers
  # matches to either the reference sequence, or query sequence, of an inversion. 
  my ($blastReport, $eThresh, $origStart, $origEnd, $orthoPos,
      $region, $match, $invList, $nameToIndex) = @_;
  $hitNumber = 0;
  if ($origEnd < $origStart) {
    my $temp = $origEnd;
    $origEnd = $origStart;
    $origStart = $temp;
  }
  my $refSeqLength = $origEnd - $origStart;
  while ( $result = $$blastReport->next_result ) {
    while ($hit = $result->next_hit ) {
      while ($hsp = $hit->next_hsp ) {
	print "checking hit: $hitNumber\n";
	++$hitNumber;
	($rBegin, $rEnd) = $hsp->range('query');
	($qBegin, $qEnd) = $hsp->range('sbjct');
	$score = $hsp->score();
	$evalue = $hsp->evalue();
	if ($evalue < $eThresh) {
	  $seqName = $hsp->seq('sbjct')->display_id();
	  if (! exists ${$match}{$seqName}) {
	    my $strand = $hsp->strand('sbjct');
	    #	    print "got strand match $seqName: $strand\n";
	    my ($humanStartPos, $humanEndPos);
	    $humanStartPos = @{$$orthoPos{$seqName}}[0];
	    $humanEndPos =  @{$$orthoPos{$seqName}}[1];
	    # add the position if they map to the same area on the human chrom
	    my ($startDiff, $endDiff);
	    $startDiff = abs($origStart - $humanStartPos);
	    $endDiff = abs($origEnd - $humanEndPos);
	    print "seqname: $seqName\n";
	    print "startDiff $startDiff ($origStart - $humanStartPos) endDiff: $endDiff ($origEnd - $humanEndPos)\n";
	    if ((abs($startDiff) < 1000) and
		(abs($endDiff) < 1000)) {

	      # compare the lengths of the two to check that we're still in the same ballpark
	      
	      print "$seqName\t\t$origStart\t$humanStartPos\t$origEnd\t$humanEndPos\n";
	      $$match{$seqName} = $hsp;
	      # add the aligned sequence to the match list.
	      my $hitIndex;
	      $hitIndex = $$nameToIndex{$seqName};

	      my $qrySeqLength = $$invList[$hitIndex][-1]->length();
	      $okLength = 0;
	      print "got refseq length: $refSeqLength qry: $qrySeqLength\n";
	      
	      if ($refSeqLength > 0 and 
		  $qrySeqLength / $refSeqLength < 3 and
		  $qrySeqLength > 0 and
		  $refSeqLength / $qrySeqLength < 3) {
		$okLength = 1;
	      }
	      print "comparing '$$invList[$hitIndex][-3]' and '$seqName' \n";
	      if ($$invList[$hitIndex][-3] eq $seqName) {
		# this is the reference sequence, add the query
		my $qryName = $$invList[$hitIndex][-2];
		if (! exists $$match{$qryName} ) {
		  print "adding inversion aligned sequence $qryName \n";
		  $$match{$qryName} = 1;
		}
	      }
	      else {
		print "comparing '$$invList[$hitIndex][-2]' and '$seqName' \n";
		if ($$invList[$hitIndex][-2] eq $seqName ) {
		  my $refName = $$invList[$hitIndex][-2];
		  if (! exists $$match{$refName} ) {
		    print "adding inversion aligned sequence $refName \n";
		    $$match{$refName} = 1;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

sub BinMatchList {
  my ($matches, $nameToIndex, $bins, $nameToBin,
      $binIndex, $binnedInversions, $orthoPos, $toStart, $toEnd, $invList) = @_;
  my $matchKey;

  # matches is a list of names where the inversion
  # was found.  assign all of the names the bin index.
  # It is possible that a name was already assigned a bin
  # in this case the bin should be reused.
  # However, if this matchlist corresponds to bins that have
  # more than one bin, that's a problem.  Try to find out that here.
  @matchNames = keys %{$matches};


  # matches contains a hash of sequences that align to a query sequence, and map
  # to the same location on the human sequence.  Each bin contains a list of 
  # sequences that have the same properties, but may not have been found using the 
  # original query.  The process of binning inversions involves linking together 
  # one match list if any of its elments are found in another match list (bin).
  my $binFound = 0;
  my $binUpdated = 0;
  foreach $matchName (keys %{$matches}) {
    my $seqIndex = $$nameToIndex{$matchName};
    foreach $binIndex (0 .. $#{$bins}) {
      @binkeys = keys %{$$bins[$binIndex]};
      if (exists ${$$bins[$binIndex]}{$seqIndex}) {
	$existingName = ${$$bins[$binIndex]{$seqIndex}}[1];
	print "found an existing bin $binIndex $matchName ($seqIndex), $existingName \n";
	# we're about to link two bins together. The two sequences should all
	# map to the same orthologous position.  Double check that here.
	
	my $binCloseEnough = 1;
	foreach $binnedSeq (@binkeys) {
	  my ($binHumanStart, $binHumanEnd);
	  $binSeqName = ${$$bins[$binIndex]{$binnedSeq}}[1];
#	  print "binSeqName is $binSeqName\n";
	  $binHumanStart = @{$$orthoPos{$binSeqName}}[0];
	  $binHumanEnd   = @{$$orthoPos{$binSeqName}}[1];
	  my ($startDiff, $endDiff);
	  $startDiff = abs($binHumanStart - $toStart);
	  $endDiff   = abs($binHumanEnd   - $toEnd);
#	  print "$binnedSeq checking mapped coordinates $toStart $binHumanStart $startDiff " .
	#    " $toEnd $binHumanEnd $endDiff\n";
	  if ($startDiff > 2000 or $endDiff > 2000) {
	    $binCloseEnough = 0;
	  }
	}
	if ($binCloseEnough == 1 ) {
	  # found a bin for one of the inversions.  Add
	  # all other inversions to it.
	  $binFound = 1;
	  my $evalue = $$matches{$matchName};
	  print "found existing bin for $seqIndex in bin $binIndex with value $evalue\n";
	  foreach $matchName2 (keys %${matches}) {
	    my $seqIndex2 = $$nameToIndex{$matchName2};
	    # @{bins[$binIndex]} is a bin
	    # @{bins[$binIndex]}{matchName2} is a reference to a matched
	    # sequence inside the bin.
	    @{$$bins[$binIndex]{$seqIndex2}} = ($matches{$matchName}, $matchName2);
	    $$binnedInversions{$matchName2} = 1;
	    $binUpdated = 1;
	  }
	  last;
	}
      }
    }
    if ($binFound == 1) {
      last;
    }
  }
  if ($binFound == 0) {
    # No bin found for this cluster. Create one.
    $lastIndex = $#{$bins} + 1;
    print "no bins found, creating one at index $lastIndex\n";
    my @mk = keys %matches;
    my $nMatches = scalar @mk;
    print "nmatches: $nMatches\n";
    if ($nMatches > 0) {
      push @{$bins}, {};
      foreach $matchName (keys %matches) {
	$seqIndex = $$nameToIndex{$matchName};
	@{$$bins[-1]{$seqIndex}} = ($matches{$matchName}, $matchName);
	$$binnedInversions{$matchName} = 1;
      }
    }
    $binUpdated = 1;
  }
  return $binUpdated;
}

