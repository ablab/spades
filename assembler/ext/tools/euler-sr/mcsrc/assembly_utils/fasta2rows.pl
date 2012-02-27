#!/usr/bin/perl -w

use strict;
use warnings;

my $fasta_file = shift @ARGV;

my $first_row = 0;
open(FASTA,$fasta_file) or die('cannot read fasta file $fasta_file\n'); 
while (<FASTA>) {
    chomp;
    if (m/^>/) { 
	chomp;
	s/^>//; 
	print "\n" if $first_row; 
	$first_row = 1;
#	my @l = split; 
#	print join(" ",@l),"\t"; 
	print "$_\t";
    }
    else { print; }
} close(FASTA);
