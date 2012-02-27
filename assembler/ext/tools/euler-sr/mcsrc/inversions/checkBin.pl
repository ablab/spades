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
#$binDataFileName = shift @ARGV;
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
#print "starting with option binDataFileName = $binDataFileName\n";
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

@sequences = ();
while ($invSeq = $invSeqReader->next_seq) {
  push @sequences, $invSeq;
  $did = $invSeq->display_id();
  $len = $invSeq->length();
  print "got seq: $did $len\n";
}
#@data = ();
#open(BINDATA, $binDataFileName) or die "cannot open $binDataFileName\n";

#while (<BINDATA>) {
#  $line = $_;
#  @values = split(/\s+/, $line);
#  # the extra value is a place holder for the orientation of the match
#  push @data, [@values, 0];
#}
#close BINDATA;


open (SPECIES_LIST, $seqSetFileName) or die "cannot open $seqSetFileName \n";
$line = <SPECIES_LIST>;
@species = split(/\s+/, $line);
close SPECIES_LIST;
%species = ();
foreach $s (@species) {
  $species{$s} = 0;
}

%orthoPos = ();
common::ReadOrthoCoordinates($orthFileName, \%orthoPos);

#if ($#data != $#sequences) {
#  print "$#data $#sequences\n";
#  print "error: there should be just as many data items as sequences\n";
#  exit(1);
#}

# hack to get rid of duplicates
%seqNames = ();
$seqIndex = 0;
while  ($seqIndex < $#sequences) {
  $seqName = $sequences[$seqIndex]->display_id();
  if (exists $seqNames{$seqName}) {
    # this sequence is a duplicate, remove it.
    splice @sequences, $seqIndex, 1;
    splice @data, $seqIndex, 1;
  }
  else {
    $seqNames{$seqName} = $seqIndex;
    $seqIndex++;
  }
}

# find out what species are present in the inversion
for ($s = 0; $s < scalar @sequences; $s++ ) {
  $key = $sequences[$s]->display_id();
  print "logging $key: $data{$key}[0]\n";
  $species{@{$data{$key}}[0]} = 1;
}

#foreach $dataIndex (0 .. $#data) {
#  $species{$data[$dataIndex][0]} = 1;
#  $species{$data[$dataIndex][3]} = 1;
#}
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

    foreach $withIndex ($newIndex .. $#sequences) {
      $seqName = @sequences[$withIndex]->display_id();
      print "with $seqName \n";
      $qryDBName = @without[$w] . ".fasta";
      ($evalue, $hspRef) = FindInDB(@sequences[$withIndex], $qryDBName,
				    $regionName, $dbDir, $eThresh);
      if ($hspRef != 0) {
	my $refStart = $$hspRef->start('query');
	my $refEnd   = $$hspRef->end('query');
	my $qryStart = $$hspRef->start('sbjct');
	my $qryEnd   = $$hspRef->end('sbjct');
	my $qrySeq   = $$hspRef->seq('sbjct');
	my $qryName  = $$hspRef->seq('sbjct')->display_id();
	my $qryStrand = $$hspRef->strand('sbjct');
	my $refSeqName  = $data[$seqIndex][6];
	my ($ref, $qry, $rs, $re) = common::ParseName($refSeqName, \@sequences);
	my $refSpecies = $ref;
	my $qrySpecies = @without[$w];

	my ($humanStartPos, $humanEndPos);
	my $refIndex = $seqNames{$seqName};
	my $refLength = @sequences[$refIndex]->length;

	my $qrySubSeqStart = $qryStart - $refStart;
	my $qrySubSeqEnd   = $qryEnd + ($refLength - $refEnd);

	print "got query strand: $qryStrand pos $qryStart -> $qrySubSeqStart $qryEnd -> $qrySubSeqEnd\n";
	$humanStartPos = common::GetOrthoPosOnly(@without[$w],
						 "human",
						 $regionName,
						 $qrySubSeqStart, 0, $sqlDb);

	$humanEndPos = common::GetOrthoPosOnly(@without[$w],
					       "human",
					       $regionName,
					       $qrySubSeqEnd, 0, $sqlDb);

	print "human orthologous position for $qryStart: $humanStartPos \n";
	print "human orthologous position for $qryEnd: $humanEndPos \n";
	# look to see if the positions map to the same location
	my ($startDiff, $endDiff);
	my ($refStartPos, $refEndPos);
	$refStartPos = ${$orthoPos{$seqName}}[0];
	$refEndPos = ${$orthoPos{$seqName}}[1];
	
	$startDiff = abs($refStartPos - $humanStartPos);
	$endDiff = abs($refEndPos - $humanEndPos);	print "comparing $seqName \n\t$refStartPos - $humanStartPos = $startDiff ".
	  " and\n\t$refEndPos - $humanEndPos = $endDiff\n";
	if ($startDiff < 2500 and
	    $endDiff < 2500) {
	  @{$orthoPos{$qryName}} = ($humanStartPos, $humanEndPos);
	  print "storing $humanStartPos, $humanEndPos for $qryName\n";
	  $sequencesAdded  = 1;
	  # Found the inverted sequence, store it.
	  $$hspRef->seq('sbjct')->display_id(@without[$w]);
	  $evalue = $$hspRef->evalue();
	  print "found inverted sequence $qrySpecies $qryStart $qryEnd $evalue\n";
	  # Now extract the query sequence.
	  # 
	  my $qrySubSeqStart = $qryStart - $refStart;
	  my $qrySubSeqEnd   = $qryEnd + ($refLength - $refEnd);
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
	  my $qrySubStringSeq = $qryFullSeq->subseq($qrySubSeqStart, $qrySubSeqEnd);
	  $qrySubStringSeq =~ tr/actg/N/;
	  my $qrySubSeq = Bio::Seq->new(-seq=>$qrySubStringSeq,
					-display_id=>$qryName);
	  push @sequences, $qrySubSeq;
	  push @data, [$refSpecies, $refStart, $refEnd, $qrySpecies, $qryStart, $qrySeq,
		       $seqName, $qryName, $qrySubSeq];
	  # don't bother checking this species again.
	  print "removing element at $w of $#without\n";
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
  my ($sequencePtr, $qrySpecies, $seqName, $dbDir, $eThresh) = @_;
  my $dbName = $dbDir . "/". $qrySpecies;
  my $seqName = $sequencePtr->display_id();
  print "looking in $dbName\n";
  my $blastSeq = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn',
						       'database'=>$dbName,
						       'q'=>'-1',
						       'W'=>'8',
						       'e'=>"$eThresh",
						      );
  my $seqPrinter = Bio::SeqIO->new(-format=>"fasta");

  my $tmpOut = Bio::SeqIO->new(-file=>">tmpcf.fasta");

  $blastRes  = $blastSeq->blastall($sequencePtr);
  my $hsp = 0;
  $minEValue = GetMinEValue(\$blastRes, \$hsp);
  if ($minEValue < 0 or
      $minEValue > $eThresh or
      $hsp eq 0) {
    print "no match \n";
    return (0,0);
  }
  else {
    print "match with evalue: $minEValue \n";
    return ($minEValue, \$hsp);
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

sub GetMinEValue {
  my ($blastResPtr, $hspRef) = @_;
  $minEValue = -1;
  while ( $result = $$blastResPtr->next_result ) {
    while ($hit = $result->next_hit ) {
      while ($hsp = $hit->next_hsp ) {
	if ($minEValue == -1 or $minEValue > $hsp->evalue()) {
	  $minEValue = $hsp->evalue();
	  $$hspRef = $hsp;
	}
      }
    }
  }
  return $minEValue;
}

