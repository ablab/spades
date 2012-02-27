#!/usr/bin/env perl

use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use common;
use DBI;

if ($#ARGV < 2) {
  print "usage: make_grid2.pl fastafile sequences outFile [-s seqNumber] [-b binNumber] [-v (verbose)]\n";
  exit(0);
}
print "running make_grid\n";
$seqFileName = shift @ARGV;
$specNameFile = shift @ARGV;
$outFile = shift @ARGV;

print STDERR "option seqFileName: $seqFileName \n";
print STDERR "option specNameFile: $specNameFile\n";
print STDERR "option outFile: $outFile\n";
$binNumber = 0;
$printSequences = 0;
$verbose=0;
$seqNumber = 0;

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
if ($verbose) {
	print "got species @species\n";
}

$invSeqReader = Bio::SeqIO->new(-file => $seqFileName,
				-format => "fasta" );

@sequences = ();
while ($invSeq = $invSeqReader->next_seq) {
  push @sequences, $invSeq;
}

# build a list of representative sequences for each inversion.
%presentSeq = ();
%presentOrientation = ();

# list of full length sequences
%fullSequences = ();

my $orientation;
for $seqIndex (0 .. $#sequences) {
  $seqName = @sequences[$seqIndex]->display_id();
  $seqName =~/(.*)\.(\d+)\.(\d+)/;
  $ref = $1;
  if (! exists $presentSeq{$ref}) {
    $presentSeq{$ref} = $seqIndex;
  }
}

open (NAMEFILE, "$specNameFile") or die "cannot open $specNameFile\n";
$line = <NAMEFILE>;
@allSequences = split(/\s+/, $line);
close NAMEFILE;

$blastOut = "tmp.txt";
if ($verbose > 1) {
  $blastOut = 0;
}

$eThresh = 0.001;
$blastOut = "tmpOut";
$blastSA = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn',
						 'outfile'=>$blastOut,
						 'q'=>'-2',
						 'W'=>'7',
						 'e'=>"$eThresh",
						);

#print "@allSequences\n";
@k = keys  %presentSeq;
foreach $ref (keys %presentSeq) {
#  print "searching ref seq: $ref\n";
  @orientations = ();
  $refSeq = @sequences[$presentSeq{$ref}];
  $refSeq->display_id("$seqNumber");
  foreach $qry (@allSequences) {
    if ($seq ne $ref) {
      $orientation = 2;
      if (exists $presentSeq{$qry}) {
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
		$length = $rEnd - $rBegin;
		print "$ref $qry align score : $evalue, length: $length\n";
	      }
	      if ($evalue < $eThresh) {
		$strand = $hsp->strand('sbjct');
		if ($verbose) {
		  print "ro: $ref $qry $strand \n";
		}
		if ($strand == 1) {
		  $orientation = 1;
		}
		else {
		  $orientation = 0;
		}
	      }
	    }
	  }
	}
      }
    } 
    else {
      # $seq == $ref, they are in the same orientation for sure.
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

print CHARFILE "bin: $binNumber 0 to 0\n";



foreach $seq (@allSequences) {
  if (exists $presentSeq{$seq}) {
    $len = scalar @{$presentCharList{$seq}};
    printf CHARFILE "%18s\t0\t0\t $len", $seq ;
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
