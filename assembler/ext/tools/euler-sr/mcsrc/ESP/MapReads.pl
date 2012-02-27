#!/usr/bin/env perl

use ReadLibrary;
use IO::Handle;
if ($#ARGV < 2) {
  print "\nusage: $0 genomedb readsfile mappedOut [rawBlastOutput]\n\n";
  print "Maps the location of every read in readsfile to genomedb\n";
  print "mappedOut contains the map locations of each read.\n\n";
  exit(0);
}
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;

$genomeDB = shift @ARGV;
$readsFile = shift @ARGV;
$mappedFile = shift @ARGV;
open(MAP, ">$mappedFile") or die "cannot open $mappedFile\n";
$blastIO = 0;
$rawBlastOutput = "out.blast";
if ($#ARGV >= 0) {
  $rawBlastOutput = shift @ARGV;
}
$seqReader = Bio::SeqIO->new(-file => $readsFile,
			     -format => "fasta" );

$pid = $$;
$blastOut = "blast.tmp.$pid";
# configure blast parameters and interface
$eThresh = 0.00001;
$blastall = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn',
						  'outfile'=>$blastOut,
						  'database'=>$genomeDB,
						  'W'=>'11',
						  'E'=>'1',
						  'e'=>"$eThresh"
						 );

# Determine the trimming of every read
print "determining read trimming ";
%readTrimming = ();
%readLength   = ();
$nreads = 0;
STDOUT->autoflush(1);
while ($read = $seqReader->next_seq) {
  $read->alphabet("dna");
  ($startSeq, $endSeq) = ReadLibrary::FindTrimStart($read->seq());
  $endSeq = $read->length - $endSeq;
  ($a,$b,$c,$name) = ReadLibrary::ParseABITitle($read->display_id);
  @{$readTrimming{$name}} = ($startSeq, $endSeq);
  $readLength{$name} = $read->length;
  if ($nreads % 500 == 499) {
    print ".";
  }
  $nreads++;
}
print "\n";

print "starting blast ...\n";
$blastReport = $blastall->blastall($readsFile);
print "done!\n";


while ( $result = $blastReport->next_result ) {
  $nResults = $result->num_hits;
  print "parsing $nResults results\n";
  $seqTitle = $result->query_name;
  while ($hit = $result->next_hit ) {
    # Find what read this hit corresponds to

    print "got title : $seqTitle\n";
    ($baseName,$dir,$type,$seqName) = ReadLibrary::ParseABITitle($seqTitle);
    ($trimStart, $trimEnd) = ($readTrimming{$seqName}[0],$readTrimming{$seqName}[1]);

    @qryStartLocations = ();
    @qryEndLocations = ();
    @refStartLocations = ();
    @refEndLocations = ();
    @scores = ();
    @strands = ();
    @identity = ();
    GetMatchList(\$hit, $eThresh, \@refStartLocations, \@refEndLocations,
		 \@qryStartLocations, \@qryEndLocations, \@scores, \@strands,
		 \@identity);
    if ($#qryStartLocations >= 0) {
      print "found a match\n";
      print MAP ">$seqTitle\n";
      for ($i = 0; $i <= $#qryStartLocations; $i++) {
	$refLen = $refEndLocations[$i] - $refStartLocations[$i];
	print MAP "$refLen\t$refStartLocations[$i]\t$refEndLocations[$i]";
	print MAP "\t$qryStartLocations[$i]\t$qryEndLocations[$i]";
	print MAP "\t$scores[$i]\t$strands[$i]\t$trimStart\t$trimEnd\t$identity[$i]\t";
	print MAP "$type\t$readLength{$seqName}\n";
      }
    }
    if ($blastIO != 0) {
      $blastIO->write_result($blastReport);
    }
  }
}

close MAP;

sub GetMatchList {
  # Given an inverted sequence, find all other inversions that have a similar
  # sequence, and map to a similar location on the human genome.
  # Each inversion has a sequence in a reference and query species.  This considers
  # matches to either the reference sequence, or query sequence, of an inversion.
  my ($hit, $eThresh,
      $refStartLocations, $refEndLocations,
      $qryStartLocations, $qryEndLocations, $scores, $strands,
     $identity) = @_;
  $hitNumber = 0;
  if ($origEnd < $origStart) {
    my $temp = $origEnd;
    $origEnd = $origStart;
    $origStart = $temp;
  }
  my $refSeqLength = $origEnd - $origStart;
  while ($hsp = $$hit->next_hsp ) {
    ++$hitNumber;
    ($rBegin, $rEnd) = $hsp->range('query');
    ($qBegin, $qEnd) = $hsp->range('sbjct');
    $score = $hsp->score();
    $evalue = $hsp->evalue();
    if ($evalue < $eThresh) {
      push @$refStartLocations, $rBegin;
      push @$refEndLocations, $rEnd;
      push @$qryStartLocations, $qBegin;
      push @$qryEndLocations, $qEnd;
      push @$scores, $evalue;
      push @$strands, $hsp->strand('sbjct');
      push @$identity, $hsp->frac_identical('query');
    }
  }
}


system("rm $blastOut");
