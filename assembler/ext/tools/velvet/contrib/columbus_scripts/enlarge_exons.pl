#!/usr/bin/perl -W

use strict;
use Bio::Perl;
use Bio::Tools::GFF;

my $gffio = Bio::Tools::GFF->new(-fh => \*STDIN, -gff_version => 2);
my $out = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 2);
my $feature;

#Record entries
while($feature = $gffio->next_feature()) {
	$feature->start($feature->start - 100);
	$feature->end($feature->end + 100);
	$out->write_feature($feature);	
}
