#!/usr/bin/perl -W

use strict;
use Bio::Perl;
use Bio::Tools::GFF;

my $gffio = Bio::Tools::GFF->new(-fh => \*STDIN, -gff_version => 2);
my $out = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 2);
my $feature;
my @features = ();

#Record entries
while($feature = $gffio->next_feature()) {
	push @features, $feature; 
}
$gffio->close();

#Sort them
my @sorted_features = sort {
	if ($a->seq_id gt $b->seq_id) {return 1;}
	elsif ($b->seq_id gt $a->seq_id) {return -1;}
	else {return $a->start <=> $b->start;}
} @features;

#Merge redundant elements
my $current_index = 0;
my $next_index = 1; 
while ($next_index <= $#sorted_features) {
	my $current_feature = $sorted_features[$current_index];
	my $next_feature = $sorted_features[$next_index];

	if ($next_feature->seq_id gt $current_feature->seq_id
	    || $next_feature->start > $current_feature->end) {
		$current_index = $next_index;
		$next_index++;
		next;	
	}

	if ($next_feature->end > $current_feature->end) {
		$current_feature->end($next_feature->end);
	}

	undef $sorted_features[$next_index];
	$next_index++;
}

#Output
my $index;
for ($index = 0; $index < @sorted_features; $index++) {
	if (defined $sorted_features[$index]) {
		$out->write_feature($sorted_features[$index]);	
	}
}

