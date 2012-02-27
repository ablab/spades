#!/usr/bin/env perl

use Bio::Seq;
use Bio::SeqIO;
use common;

if ($#ARGV < 2) {
  print "usage: $0 binCoordsFile seqdir outbsae\n";
  exit(0);
}

$binCoordFile = shift @ARGV;
$seqDir       = shift @ARGV;
$outBase      = shift @ARGV;
print "got outbase: $outBase\n";
$outFile = "";
if ($outBase eq "-f") {
  $outFile = shift @ARGV;
  print "using outfile: $outFile, not outbase\n";
}
else {
  print "using outbase: $outBase\n";
}

%sequences = ();

open (BINFILE, "$binCoordFile") or die "cannot open $binFile\n";


$binNumber = -1;
if ($outFile ne "") {
  $binOut = new Bio::SeqIO('-file'=>">$outFile");
  print "opening binout to $outFile\n";
}
while (<BINFILE>) {
  $line = $_;
  if ($line =~ /bin/) {
    if ( $outFile eq "") {
      # start a new bin
      $binNumber++;
      $binOut = new Bio::SeqIO('-file'=>">$outBase.$binNumber.fasta");
    }
  }
  else {
    $line =~ /\s*(\S+)\s+(\d+)\s+(\d+)/;
    $spec = $1; $start = $2; $end = $3;
    print "subseq for $spec, $start, $end: $line";
    $seq = GetSequence($spec, $seqDir, \%sequences);
    $invSeqStr = $seq->subseq($start, $end);
    $invSeq = new Bio::Seq(-seq=>$invSeqStr);
    $id = "$spec.$start.$end";
    $invSeq->display_id($id);
    $binOut->write_seq($invSeq);
  }
}



sub GetSequence {
  my ($speciesName, $seqDir, $sequences) = @_;
  if (exists $$sequences{$speciesName}) {      
    return $$sequences{$speciesName};
  }
  else {
    $seqName = "$seqDir/$speciesName" . ".fasta";
    print "getting sequence $seqName\n";
    $seqReader = new Bio::SeqIO('-file'=>"$seqName");
    $$sequences{$speciesName} = $seqReader->next_seq;
    return $$sequences{$speciesName};
  }
}
