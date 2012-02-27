#!/usr/bin/env perl
use common;
use Bio::SeqIO;
BEGIN{
  if (! -e "/state/partition1/bl2seq") {
    `cp /home/fwg.mchaisso/software/blast-2.2.13/bin//bl2seq /state/partition1/bl2seq`;
    $ENV{"BLASTDIR"} = "/state/partition1";
  }
}

use Bio::Tools::Run::StandAloneBlast;
use common;

if ($#ARGV < 2) {
  print "usage: $0 fragmentsFile blastdb outname\n";
  exit(1);
}

$queryFile = shift @ARGV;
$sbjctFile = shift @ARGV;
$outName   = shift @ARGV;
# these parameters took a while to tune.
$blastOut = common::CreateTempFileName(".tophit");
$blastObj = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn',
						  'outfile'=>$blastOut,
						  'q'=>'-1',
						  'W'=>'8',
						  'e'=>'0.00001',
						  'U'=>'T'
						 );

$exepath = $blastObj->executable('bl2seq');

eval {
  $report = $blastObj->bl2seq($queryFile, $sbjctFile);
};
if ($@) {
 exit 1;
}
($rb, $re, $qb, $qe, $ev) = GetBestHitLocation(\$report);
open(OUT, ">$outName") or die "cannot open $outName\n";
print OUT "$rb $re $qb $qe $ev\n";
close OUT;
`rm -f $blastOut`;

###############################################################################
# functions
###############################################################################

sub GetBestHitLocation {
  my ($blastReport) = @_;
  my ($rb, $re, $qb, $qe, $ev);
  $rb = -1; $re = -1; $qb = -1; $qe = -1;
  $numHits = 0;
  while ( $result = $$blastReport->next_result ) {
    while ($hit = $result->next_hit ) {
      while ($hsp = $hit->next_hsp ) {
	($rb, $re) = $hsp->range('query');
	($qb, $qe) = $hsp->range('sbjct');
	$ev = $hsp->evalue();
	$numHits++;
      }
    }
  }
  if ($numHits == 1) {
    return ($rb, $re, $qb, $qe, $ev);
  }
  else {
    return (-1,-1,-1,-1);
  }
}



