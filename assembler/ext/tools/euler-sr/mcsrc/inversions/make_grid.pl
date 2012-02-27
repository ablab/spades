#!/usr/bin/env perl

use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use common;
use DBI;

if ($#ARGV < 4) {
  print "usage: make_grid.pl fastafile sequences orthoPosFile seqDir outFile [-s seqNumber] [-n binNumber] [-v (verbose)]\n";
  exit(0);
}
print "running make_grid\n";
$seqFileName = shift @ARGV;
$specNameFile = shift @ARGV;
$orthoPosFile = shift @ARGV;
$seqDir  = shift @ARGV;
$outFile = shift @ARGV;

print STDERR "option seqFileName: $seqFileName \n";
print STDERR "option specNameFile: $specNameFile\n";
print STDERR "option orthoPosFile: $orthoPosFile\n";
print STDERR "option seqDir: $seqDir\n";
print STDERR "option outFile: $outFile\n";
$binNumber = 0;
$printSequences = 0;
$verbose=0;
$seqNumber = 0;
$dataFile   = "$orthoPosFile.data";
$refSeqFile = "$orthoPosFile.fasta";


while ($#ARGV >= 0) {
  $val = shift @ARGV;
  if ($val eq "-s") {
    $seqNumber = shift @ARGV;
  }
  if ($val eq "-n") {
    $binNumber = shift @ARGV;
  }
  if ($val eq "-v") {
    $verbose++;
  }
}

open (CHARFILE, ">>$outFile") or die "cannot open $outFile\n";
@species = ();
if ($specNameFile ne "") {
  common::ReadSpeciesFile($specNameFile, \@species);
}
$nspec = $#species;
#print "read $nspec species\n";

$refSeqFile   = "$orthoPosFile.fasta";

$refSeqWriter = Bio::SeqIO->new(-file=>">>$refSeqFile",
				-format=>"fasta");

$invSeqReader = Bio::SeqIO->new(-file => $seqFileName,
				-format => "fasta" );


#open (DATAOUT, ">>$dataFile") or die "Cannot open $dataFile \n";



@sequences = ();
while ($invSeq = $invSeqReader->next_seq) {
  push @sequences, $invSeq;
}

%orthoCoords = ();
common::ReadOrthoCoordinates($orthoPosFile, \%orthoCoords);

# build a list of representative sequences for each inversion.
%presentSeq = ();
%presentOrientation = ();
%presentCoords = ();

# list of full length sequences
%fullSequences = ();

my $orientation;
for $seqIndex (0 .. $#sequences) {
  $seqName = @sequences[$seqIndex]->display_id();
  if (common::IsNamePair($seqName)) {
    if ($verbose) {
      print "parsing $seqName\n";
    }
    if (common::IsNameQuery($seqName) == 0) {
      ($ref) = common::ParseName($seqName, \@species);
      if ($verbose) {
	print "$seqName setting orientation $ref 1 \n";
      }
      $orientation = 1;
    }
    else {
      ($query, $ref) = common::ParseName($seqName,\@species);
      $orientation = -1;
    }
    $presentCoords{$ref} =  [ @{$orthoCoords{$seqName}}] ;
    if ($verbose) {
      print "got coords for $ref: @{$presentCoords{$ref}} \n";
    }
  }
  else {
#    print "$seqName is not a pair\n";
    ($ref) = common::ParseNCBIName($seqName);
    $orientation = 1;
    @{$presentCoords{$ref}} = (0, 0);
  }
  if (! exists $presentSeq{$ref}) {
    $presentSeq{$ref} = $seqIndex;
    $presentOrientation{$ref} = $orientation;
  }
}

open (NAMEFILE, "$specNameFile") or die "cannot open $specNameFile\n";
$line = <NAMEFILE>;
@allSequences = split(/\s+/, $line);
close NAMEFILE;
if ($verbose) {
  foreach $ref (@allSequences) {
    if (exists $presentCoords{$ref} ) {
      print "got present coords: $ref @{$presentCoords{$ref}} $presentOrientation{$ref}\n";
    }
  }
}

$blastOut = "tmp.txt";
if ($verbose > 1) {
  $blastOut = 0;
}

$eThresh = 0.001;
$blastOut = "tmpOut";
$blastSA = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn',
						 'outfile'=>$blastOut,
						 'q'=>'-1',
						 'W'=>'7',
						 'e'=>"$eThresh",
						);

#print "@allSequences\n";
@k = keys  %presentSeq;
foreach $ref (keys %presentSeq) {
#  print "searching ref seq: $ref\n";
  @orientations = ();
  $refSeq = @sequences[$presentSeq{$ref}];
  $refName = $refSeq->display_id();
  $refSeq->display_id("$seqNumber");
  $refSeqWriter->write_seq($refSeq);
#  print DATAOUT "$seqNumber $refName @{$presentCoords{$ref}}\n";
  my ($rSpec, $qSpec, $refStart, $refEnd) = common::ParseName($refName);

  foreach $qry (@allSequences) {
    if ($seq ne $ref) {
      $orientation = 2;
      $qrySeq = @sequences[$presentSeq{$qry}];
      if (exists $presentSeq{$qry}) {
	my ($rl, $ql);
	$rl = $refSeq->length();
	$ql = $qrySeq->length();
	print STDERR "rl: $rl, ql: $ql\n";
	$rw = Bio::SeqIO->new(-file=>">ref.fasta", -display_id=>"ref");
	$rw->write_seq($refSeq);
	$qw = Bio::SeqIO->new(-file=>">qry.fasta", -display_id=>"qry");
	$qw->write_seq($qrySeq);
	

	$qrySeq = $sequences[$presentSeq{$qry}];

	$bl2seqReport = $blastSA->bl2seq($refSeq, $qrySeq);
	while ( $result = $bl2seqReport->next_result and $orientation == 2) {
	  while ($hit = $result->next_hit and $orientation == 2) {
	    while ($hsp = $hit->next_hsp and $orientation == 2) {
	      ($rBegin, $rEnd) = $hsp->range('query');
	      ($qBegin, $qEnd) = $hsp->range('sbjct');
	      $score  = $hsp->score();
	      $evalue = $hsp->evalue();
	      if ($verbose) {
		print "$ref $qry align score : $evalue \n";
	      }
	      if ($evalue < $eThresh) {
		$strand = $hsp->strand('sbjct');
		if ($verbose) {
		  print "ro: $ref $presentOrientation{$ref}  $qry $presentOrientation{$qry} $strand \n";
		}
		$orientation = RelativeOrientation($presentOrientation{$ref},
						   $presentOrientation{$qry},
						   $strand);

		$refStr = $hsp->seq_str('query');
		$qryStr = $hsp->seq_str('sbjct');
	      }
	    }
	  }
	}
      }
      # the orientation here is still 2 try to decide if the 
      # sequence was deleted.
      # load up the full reference sequence to try to find in 
      # other sequence
      if (! exists $fullSequences{$ref} ) {
	my $seqFileName = "$seqDir/$ref.fasta";
#	print "loading $seqFileName\n";
	my $refSeqReader = Bio::SeqIO->new(-file => $seqFileName,
					   -format => "fasta" );
	my $fullRefSeq = $refSeqReader->next_seq;
	$fullSequences{$ref} = $fullRefSeq;
      }
      my $fullRefSeq = $fullSequences{$ref};
      my $qrySeqName = "$seqDir/$qry.fasta";
      my ($refStart, $refEnd);
      $refStart = $presentCoords{$ref}[0];
      $refEnd   = $presentCoords{$ref}[1];
#      print "checking for deletion got coords: $refStart, $refEnd.\n";
      if ($orientation == 2 and $refStart != 0 and $refEnd != 0) {
	$refNamehere = $fullRefSeq->display_id();
#	print "$ref, $refNamehere checking deletion in $qry\n";
	($before, $after) = (0,0);
#	($before, $after) = common::CheckForDeletion($fullRefSeq, "$seqDir/blastdb/$qry.fasta",
#						     $refStart, $refEnd, 200, 400, 0.001);

	if ($before != 0 and $after != 0) {
	  if (($after - $before) < ($refEnd - $refStart)) {
	    $orientation = 3;
	  }
	  elsif (($after - $before ) > 2 * ($refEnd - $refStart)) {
	    $orientation = 4;
	  }
	}
      }
    }
    else {
      $orientation = 1;
    }
    push @orientations, $orientation;
  }
  for $o (0 .. $#orientations) {
    if ($orientations[$o] == -1) {
      $orientations[$o] = 0;
    }
  }
  @{$presentCharList{$ref}} = @orientations;
}


if (exists $presentCoords{"human"} ) {
  my ($start, $end);
  $start = @{$presentCoords{"human"}}[0];
  $end = @{$presentCoords{"human"}}[1];
  print CHARFILE "bin: $binNumber $start to $end\n";
}
else {
  print CHARFILE "bin: $binNumber 0 to 0\n";
}


foreach $seq (@allSequences) {
  if (exists $presentCoords{$seq}) {
    $len = @{$presentCharList{$seq}};
    printf CHARFILE "%18s\t@{$presentCoords{$seq}}[0]\t@{$presentCoords{$seq}}[1]\t $len", $seq ;
    foreach $c (0 .. $#{$presentCharList{$seq}}) {
      printf CHARFILE " %2d", ${$presentCharList{$seq}}[$c];
    }
    printf CHARFILE "\n";
  }
}

close CHARFILE;
sub RelativeOrientation {
  my ($refStrand, $queryStrand, $relative) = @_;

  if ($refStrand == $queryStrand) {
    return $relative;
  }
  else {
    if ($relative == 1) { return -1; }
    else { return 1; }
  }
}
