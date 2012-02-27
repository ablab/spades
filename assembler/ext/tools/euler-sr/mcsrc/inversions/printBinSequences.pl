#!/usr/bin/env perl

use common;
use Bio::Seq;
use Bio::SeqIO;

if ($#ARGV != 2) {
  print "usage: $0 invfilename binfile binseqbase\n";
  exit(0);
}

$invFileName = shift @ARGV;
$binFile = shift @ARGV;
$binSeqBase = shift @ARGV;



@bins = ();

common::ReadBinFile($binFile, \@bins);


#for ($b = 0; $b <= $#bins; $b++ ) {
#  print "bin: $b\n";
#  for ($bi = 0; $bi <= $#{$bins[$b]}; $bi++) {
#    print "@{$bins[$b][$bi]}\n";
#  }
#}


@refSpecies = ();
@invList = ();
%pairs = ();



common::ReadInvFile($invFileName, \@refSpecies, \%pairs);
common::ReadFiles(\@refSpecies, \%pairs, \@invList, \%nameToSequence, \%nameToIndex, ".");
for ($b = 0; $b <= $#bins; $b++ ) {
  $out = "$binSeqBase.$b.fasta";
  $seqOut_obj = Bio::SeqIO->new(-file => ">$out");

  print "bin: $b\n";
  for ($bi = 0; $bi <= $#{$bins[$b]}; $bi++) {
    $index = $nameToIndex{$bins[$b][$bi][0]};
    $seqOut_obj->write_seq($invList[$index][-1]);
  }
}
