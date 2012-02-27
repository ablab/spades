#!/usr/bin/perl -w

###########################################################################
# Title:          seqcount.pl
# Author:         Glenn Tesler
# Description:    Read a Fasta sequence file, and count the number
#                 of sequences, lengths, etc.
# Created:        02/10/2002
# Last modified:  02/10/2002
#
# Copyright (c) 2001-2004 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
###########################################################################

sub slurpfile {
    my $counts = init_counts();

    my $len = 0;

    while (<>) {
	if (/\A>/) {
	    register_len($counts,$len);
	    $len = 0;
	} else {
	    chomp;
	    $len += length;
	}
    }

    register_len($counts,$len);

    return $counts;
}

sub init_counts {
    my $counts = {
	'min' => undef,
	'max' => undef,
	'numchar' => 0,
	'numseq' => 0,
    };

    return $counts;
}

sub register_len {
    my ($counts,$len) = @_;

    return        if (!$len);

    $counts->{'max'} = max2($counts->{'max'},$len);
    $counts->{'min'} = min2($counts->{'min'},$len);
    $counts->{'numchar'} += $len;
    $counts->{'numseq'}++;
}


sub printstats {
    my ($counts) = @_;

    if (!$counts->{'numseq'}) {
	print "Empty file!\n";
    } else {
	printf
	    "%d seqs, %d bp, %d-%d bp/seq, %.2f avg\n",
	    $counts->{'numseq'},
	    $counts->{'numchar'},
	    $counts->{'min'},
	    $counts->{'max'},
	    $counts->{'numchar'} / $counts->{'numseq'};
    }
}


# max2($a,$b)
# min2($a,$b)
# value undef is permitted for $a; forces returning $b

sub max2 {
    my ($a,$b) = @_;
    (defined $a)
	? ( ($a>$b) ? $a : $b )
	    : $b;
}

sub min2 {
    my ($a,$b) = @_;
    (defined $a)
	? ( ($a<$b) ? $a : $b )
	    : $b;
}


sub main {
    my $counts = slurpfile();
    printstats($counts);
}

main();
