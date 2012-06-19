#!/usr/bin/perl -w

use strict;
use Bio::Perl;
use Bio::Tools::GFF;
use Bio::DB::Fasta;
use Bio::PrimarySeq;

my $fasta_file = $ARGV[0];

my $gffio = Bio::Tools::GFF->new(-fh => \*STDIN, -gff_version => 2);
my $db = Bio::DB::Fasta->new($fasta_file);
my $out = Bio::SeqIO->new('-format' => 'Fasta','-fh' => \*STDOUT);
my $feature;

while($feature = $gffio->next_feature()) {
	my $start = $feature->start;
	my $end = $feature->end;
	my $strand = $feature->strand;
	my $seq_id = $feature->seq_id;
	my $ref = $db->get_Seq_by_id($seq_id);
	
	my $seq;

	if ($strand >= 0) {
		$seq = Bio::PrimarySeq->new( 
			-seq => $ref->subseq($start,$end),
                        -id  => "$seq_id:$start-$end");
	} else {
		$seq = Bio::PrimarySeq->new( 
			-seq => $ref->subseq($start,$end),
                        -id => "$seq_id:$start-$end");
	}

	$out->write_seq($seq);
}

$gffio->close();
