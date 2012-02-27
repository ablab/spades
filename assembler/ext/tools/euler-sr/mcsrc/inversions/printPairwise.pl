#!/usr/bin/env perl

use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use common;
use DBI;


if ($#ARGV != 1) {
  print "usage: $0 infile.fasta speciesfile\n";
  exit(1);
}
$infile = shift @ARGV;
$speciesFile = shift @ARGV;
$invSeqReader = Bio::SeqIO->new(-file => $infile,
				-format => "fasta" );

@sequences = ();
while ($invSeq = $invSeqReader->next_seq) {
  push @sequences, $invSeq;
}
@data = ();


open (SPECIES_LIST, $speciesFile) or die "cannot open $speciesFile \n";
$line = <SPECIES_LIST>;
@species = split(/\s+/, $line);
close SPECIES_LIST;
%species = ();
foreach $s (@species) {
  $species{$s} = 0;
}

%refSpecies = ();
%qrySpecies = ();
%refSequences = ();
%foundSpecies = ();
foreach $seq (@sequences) {
  $seqName = $seq->dipslay_id();
  my ($ref, $qry, $start, $end);
  if (IsNamePairwise($seqName)) {
    ($ref, $qry, $start, $end) = common::ParseName($seqName);
    if (common::IsNameQuery($seqName) == 0) {
      $refSpecies{$ref} = $seq;
      $foundSpecies{$ref} = 1;
    }
    else {
      $qrySpecies{$qry} = $seq;
      $foundSpecies{$qry} = 1;
    }
  }
  else {
    $ref = GetSpeciesFromName($seqName);
    $refSpecies{$ref} = $seq;
    $foundSpecies{$ref} = 1;
  }
}

$blastSA = Bio::Tools::Run::StandAloneBlast->new('program'=>'blastn',
						 'outfile'=>$blastOut,
						 'q'=>'-2',
						 'W'=>'10'
						);

@fs = keys %foundSpecies;

foreach $spec (@fs) {
  print "$spec ";
  if (exists $refSpecies{$spec}) {
    foreach $compSpec (@species) {
      if (exists $refSpecies{$compSpec}) {
	# compare the two species as both reference 
      }
      else if (exists $qrySpecies{$compSpec}) {
	# compare the species as one ref one query
      }
    }
  }
  elsif (exists $qrySpecies{$spec}) {
    foreach $compSpec (@species) {
      if (exists $refSpecies{$compSpec}) {
	# compare the two species one query, the other ref
      }
      else if (exists $qrySpecies{$compSpec}) {
	# compare the species as both query
      }
    }
  }
}




sub IsNamePairwise {
  my ($name) = @_;
  return (($name =~ /(\w+)_(\d+)_(\d+)_(\w+)/) or
	  ($name =~ /(\w+)_(\w+)_(\d+)_(\d+)/));
}

sub GetSpeciesFromName {
  my ($name) = @_;
  @vals = split(/\|/, $name);
  return @vals[0];
}
