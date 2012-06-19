#!/usr/bin/perl -w
##
##       extractContigReads.pl
##
##       Copyright 2009 Daniel Zerbino <zerbino@ebi.ac.uk>
##
##       This program is free software; you can redistribute it and/or modify
##       it under the terms of the GNU General Public License as published by
##       the Free Software Foundation; either version 2 of the License, or
##       (at your option) any later version.
##
##       This program is distributed in the hope that it will be useful,
##       but WITHOUT ANY WARRANTY; without even the implied warranty of
##       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##       GNU General Public License for more details.
##
##       You should have received a copy of the GNU General Public License
##       along with this program; if not, write to the Free Software
##       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
##       MA 02110-1301, USA.
##
##       A script to split out a single contig's reads from an LastGraph file
##       produced by the Velvet assembler.  The output is in fasta format.
##
##       Usage:  ./extractContigReads.pl <contig number> <directory> 
##
##
##       Where:  <contig number> is the number of the contig of interest.
##		 <directory> is the Velvet directory 
##
######################################################################

use strict;

my $usage = "$0 <contig number> <directory>";

my $contig_no=$ARGV[0];
my $directory=$ARGV[1];

unless($directory){die "Usage: $usage\n"};

my $graphfile = "$directory/LastGraph";

unless(-e $directory){die "$directory does not exist\n"};
unless(-e $graphfile){die "$graphfile does not exist, please re-run Velvet with the reads tracking on"};

open GRAPH, $graphfile;

$_ = <GRAPH>;
chomp;
my @data = split /\t/;
my %reads = ();
my $read_no;
my $readingSeqPath = 0;
my $recordReadIDs = 0;

unless ($data[0] >= $contig_no) {die "Contig $contig_no does not exist in the LastGraph file.\n"};

while(<GRAPH>) {
	if (/SEQ/) {
		chomp;
		@data = split /\t/;
		$read_no = $data[1];
		$readingSeqPath = 1;
		next;
	}

	if (/NR/) {
		$readingSeqPath = 0;
	}

	if ($readingSeqPath == 1) {
		chomp;
		@data = split /\t/;
		if (abs($data[0]) == abs($contig_no)) {
			$reads{$read_no} = 1;
		}
	}

	if (/^NR\t-?$contig_no\t/) {
		$recordReadIDs = 1;
	} elsif (/^NR/) {
		$recordReadIDs = 0;
	}

	if ($recordReadIDs == 1) {
		chomp;
		@data = split /\t/;
		$read_no = $data[0];
		$reads{$read_no} = 1;
	}
}

close GRAPH;

my $seqfile = "$directory/Sequences"; 
unless (-e $seqfile) {die "$seqfile does not exist, exiting.\n"};
open SEQ, $seqfile;
my $printRead = 0;

while(<SEQ>) {
	if (/>/) {
		chomp;
		@data = split /\t/;
		if ($reads{$data[1]}) {
			$printRead = 1;
			print "${data[0]}\n";
		} else {
			$printRead = 0;
		}
	} elsif ($printRead) {
		print $_;
	}
}

close SEQ;
