#!/usr/bin/env perl

############################################################################
# Title:          Illumina2Fastq.pl
# Author:         Glenn Tesler
# Created:        2009
# Last modified:  01/23/2010
#
# Copyright (c) 2009-2010 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
############################################################################

use Getopt::Long;

sub print_usage {
		print <<"END_usage";
$0 inputfile outputfile [options]
$0 inputfile1 inputfile2 outputfile [options]

    Specify two input files if they are paired reads in the same order;
    the paired reads will be interleaved on output.

    inputfile(s) are in an Illumina format (Eland, qseq, GA): tab-sep columns
        machine, run \#,               [omitted in some variations]
        lane, tile, x, y, index, read \#, read, quality string,
        ...,   [additional columns for alignment; omitted in some variations]
        filter

    outputfile is in fastq or fasta format, depending on options.

    --suffix S   (typically "1" or "2"):
    Normally it uses the "Read number" column (column 6 or 8) to produce the
    suffix.  If "-suffix S" is given, the read number column is ignored.
    For two files, use "--suffix S1,S2"

   --keepbad:  Keep reads that failed the quality filter.
               Default is to discard them.

   --append:   Append to the output file if it already exists.
               Default is to clobber the file.

   --outfasta: Output in fasta format.
   --outfastq: Output in fastq format [default].
END_usage
}

# Future options:
#   --ineland   input formats
#   --inqseq
# specify which columns should be concatenated together for the name;
# which has the paired read number;
# the read;
# the QV string;
# and the quality filter
####   --insff   this is a binary format for 454, keep it in its own converter


sub main {
		my $suffix;
		my ($filter, $keepbad);
		my $append;
		my $help;

		my $outfasta;
		my $outfastq;

		my $error;

		GetOptions("help" => \$help,
							 "suffix=s" => \$suffix,
#							 "filter" => \$filter,
							 "keepbad" => \$keepbad,
							 "append" => \$append,
							 "outfasta" => \$outfasta,
							 "outfastq" => \$outfastq,
							 );
		$filter = !$keepbad;

		if ($outfasta && $outfastq) {
				print "Can only specify one of --outfasta or --outfastq\n";
				$error++;
		} elsif (!$outfasta && !$outfastq) {
				$outfastq++;
		}

		if ($help || $error || scalar(@ARGV) < 2 || scalar(@ARGV) > 3) {
				print_usage();
				exit(0);
		}

		my $outputfile = pop @ARGV;
		my $inputfile1 = shift @ARGV;
		my $inputfile2 = shift @ARGV;  # undef if only one input file given

		my $suffix1;
		my $suffix2;
		if ($suffix =~ /^(.*),(.*)$/) {
				$suffix1 = $1;
				$suffix2 = $2;
		} else {
				$suffix1 = $suffix;
		}

		my $outmode = $append ? ">>" : ">";
		open (OFILE, $outmode, $outputfile)
				or die "Cannot open output file $outputfile";

		open (IFILE1, "<", $inputfile1)
				or die "Cannot open input file $inputfile1";

    # get column numbers
		my ($readnum_col, $filter_col) = deduce_format(\*IFILE1);
		my @cols = ($readnum_col, $filter_col);

		# statistics:
		# $stats->[0] = num. reads in input
		# $stats->[1] = num. bases in input
		# $stats->[2] = num. reads in output
		# $stats->[3] = num. bases in output
		my $stats = [0,0,0,0];

		# For paired reads:
		# $pairs->[0] = num. pairs both reads discarded
		# $pairs->[1] = num. pairs with read 1 only in output
		# $pairs->[2] = num. pairs with read 2 only in output
		# $pairs->[3] = num. pairs in output
		my $pairs = [0,0,0,0];


		if ($inputfile2) {
				open (IFILE2, "<", $inputfile2)
						or die "Cannot open input file $inputfile2";

				while (!(eof(IFILE1) && eof(IFILE2))) {
						my ($got1,$got2);
						if (!eof(IFILE1)) {
								$got1 = getread(\*IFILE1, \*OFILE, $filter, $suffix1, $outfastq, @cols, $stats);
						}
						if (!eof(IFILE2)) {
								$got2 = getread(\*IFILE2, \*OFILE, $filter, $suffix2, $outfastq, @cols, $stats);
						}
						$pairs->[$got1 + 2*$got2] ++;
				}

				close(IFILE2)
						or die "Cannot close input file $inputfile2";

		} else {
				while (!eof(IFILE1)) {
						getread(\*IFILE1, \*OFILE, $filter, $suffix1, $outfastq, @cols, $stats);
				}
		}

		close(IFILE1)
				or die "Cannot close input file $inputfile1";

		close(OFILE)
				or die "Cannot close output file $outputfile";

		print <<"ENDstats";
Input:  $stats->[0] reads with $stats->[1] bases
Output: $stats->[2] reads with $stats->[3] bases
ENDstats

    if ($inputfile2) {
				print <<"ENDpairstats";

$pairs->[0] pairs with both reads discarded
$pairs->[1] pairs with read 1 retained, read 2 discarded
$pairs->[2] pairs with read 1 discarded, read 2 retained
$pairs->[3] pairs with both reads retained
ENDpairstats
		}

}

sub getread {
		my ($fh_in,$fh_out,$filter,$suffix,$outfastq,
				$readnum_col, $filter_col,
				$stats
				) = @_;

		my $line = <$fh_in>;
		return if (!defined $line);

		chomp $line;
		my @f = split /\t/, $line;

		$stats->[0]++; # number of reads in input
		$stats->[1] += length($f[$readnum_col + 1]); # number of bases in input

		if ($filter) {
				my $pass = $f[$filter_col];
				return 0 if ($pass eq '0' || $pass eq 'N');
		}

		my ($readnum, $seq, $qual) = @f[$readnum_col .. $readnum_col + 2];
		$stats->[2]++; # number of reads in output
		$stats->[3] += length($seq); # number of bases in output

		if (defined $suffix) { $readnum = $suffix; }

		my $name = join ":", @f[0 .. $readnum_col-1];
		$name = "${name}/${readnum}";

		$seq =~ s/\./N/g;
		if ($outfastq) {
				print $fh_out "\@${name}\n${seq}\n+${name}\n${qual}\n";
		} else {
				print $fh_out ">${name}\n${seq}\n";
		}

		return 1;
}

# Variations in Illumina files include:
#  20 or 22 column versions of Eland format
#  11 column qseq format
# Read the first line from file (and then rewind).
# Deduce the column numbers to use:
#  0 .. $readnum_col-1:   columns with coordinates
#  $readnum_col:          start of 3 columns with read #, read, QV string
#  $filter_col:           did read pass quality filtering?

sub deduce_format {
		my ($fh_in) = @_;

		my $line = <$fh_in>;
		if (!defined $line) {
				die "Cannot read first file";
		}
		seek($fh_in, 0, 0); # rewind

		chomp $line;
		my @f = split /\t/, $line;

		my ($readnum_col, $filter_col);
		foreach my $col (5, 7) {
				# see if columns $col, $col+1, $col+2
				# have correct format for read number, read, quality string
				my ($readnum, $read, $qv) = @f[$col .. $col+2 ];
				next if (!($readnum eq ''         # not mate paired
									 || $readnum eq '1'     # mate paired, read 1 or 2
									 || $readnum eq '2'));

				# verify lengths seem to be for reads,
				# that the read length and qv length match,
				# and that the characters are valid

				next if (length($read) < 20);
				next if (length($read) != length($qv));

				next if ($read !~ /^[acgtnACGTN\.]*$/);
				next if ($qv !~ /^[\x{40}-\x{7f}]*$/);

				# Consistent!
				$readnum_col = $col;
				last;
		}

		if (!defined $readnum_col) {
				die "Cannot locate read number, read, quality value columns";
		}

		$filter_col = scalar(@f)-1;
		my $pass = @f[$filter_col];
		if ($pass !~ /^[01NYny]$/) {
				die "Cannot locate quality filter column";
		}

		return ($readnum_col, $filter_col); 

}

main();
1;

###############################################################################
# File formats
###############################################################################

###############################################################################
# Known input formats:
###############################################################################

# .qseq input file format:
#	Column 1: Machine name
#	Column 2: Run number
#	Column 3: Lane number
#	Column 4: Tile number
#	Column 5: X coordinate of the cluster
#	Column 6: Y coordinate of the cluster
#	Column 7: Index
#	Column 8: Read number (1 or 2)
#	Column 9: Sequence (using "." instead of "N")
#	Column 10: Quality string
#	Column 11: Did the read pass quality filtering? 0 = No; 1 = Yes

# Eland format: same as GA format below, but sometimes
# includes machine name and run number, sometimes doesn't.
# When it doesn't, all other column #'s are shifted down 2.

# GA format:
#
#    1. Machine (Parsed from Run Folder name)
# 
#    2. Run Number (Parsed from Run Folder name)
# 
#    3. Lane
# 
#    4. Tile
# 
#    5. X Coordinate of cluster
# 
#    6. Y Coordinate of cluster
# 
#    7. Index string (Blank for a non-indexed run)
# 
#    8. Read number (1 or 2 for paired-read analysis, blank for a
# single-read analysis)
# 
#    9. Read
# 
#   10. Quality string--In symbolic ASCII format (ASCII character code =
# quality value + 64)
# 
#   11. Match chromosome--Name of chromosome match OR code indicating why
# no match resulted
# 
#   12. Match Contig--Gives the contig name if there is a match and the
# match chromosome is split into contigs (Blank if no match found)
# 
#   13. Match Position--Always with respect to forward strand, numbering
# starts at 1 (Blank if no match found)
# 
#   14. Match Strand--"F" for forward, "R" for reverse (Blank if no match
# found)
# 
#   15. Match Descriptor--Concise description of alignment (Blank if no
# match found)
# 
#           * A numeral denotes a run of matching bases
# 
#           * A letter denotes substitution of a nucleotide: For a 35 base
# read, "35" denotes an exact match and "32C2" denotes substitution of a
# "C" at the 33rd position
# 
#   16. Single-Read Alignment Score--Alignment score of a single-read
# match, or for a paired read, alignment score of a read if it were
# treated as a single read. Blank if no match found; any scores less than
# 4 should be considered as aligned to a repeat
# 
#   17. Paired-Read Alignment Score--Alignment score of a paired read and
# its partner, taken as a pair. Blank if no match found; any scores less
# than 4 should be considered as aligned to a repeat
# 
#   18. Partner Chromosome--Name of the chromosome if the read is paired
# and its partner aligns to another chromosome (Blank for single-read
# analysis)
# 
#   19. Partner Contig--Not blank if read is paired and its partner aligns
# to another chromosome and that partner is split into contigs (Blank for
# single-read analysis)
# 
#   20. Partner Offset--If a partner of a paired read aligns to the same
# chromosome and contig, this number, added to the Match Position, gives
# the alignment position of the partner (Blank for single-read analysis)
# 
#   21. Partner Strand--To which strand did the partner of the paired read
# align? "F" for forward, "R" for reverse (Blank if no match found, blank
# for single-read analysis)
# 
#   22. Filtering--Did the read pass quality filtering? "Y" for yes, "N"
# for no




###############################################################################
# Output formats:
###############################################################################

# .fastq output file format:
# Line 1: @seqname
# Line 2: sequence
# Line 3: +seqname
# Line 4: quality values as ASCII string of characters 64 + log(1/p)

# .fasta output file format:
#	Line 1: >seqname
#	Line 2: sequence
