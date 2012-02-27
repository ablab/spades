#!/usr/bin/perl -w

############################################################################
# Title:          ColorEdgesByTypeOfReadCoverage.pl
# Author:         Boyko Kakaradov
# Created:        2009
# Last modified:  12/06/2009
# 
# Copyright (c) 2009 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
############################################################################


use List::Util qw[min max sum];
use strict;
use warnings;

if ($#ARGV < 2) {
    print "usage ColorEdgesByTypeOfRead.pl [-noZeroMult]  readFile intervalFile gvzFile [indexCutoff]\n";
    exit(0);
}

my $reads_file = shift @ARGV;
my $showZeroMultiplicityEdges = 1;
if ($reads_file eq '-noZeroMult'){
    $showZeroMultiplicityEdges = 0;
    $reads_file = shift @ARGV;
}
my $interval_file = shift @ARGV;
my $gvz_file = shift @ARGV;

my ($index_cutoff, $common_prefix, $first_lines) = (0,"",3);
open(READS,$reads_file) or die("cannot open reads file $reads_file\n");
while(<READS>) {
    next if not m/^>/;
    $index_cutoff ++;
    if ($first_lines > 0) {
	$common_prefix = $_ if $common_prefix eq "";
	chop $common_prefix while (not m/^\Q$common_prefix\E/);
	$first_lines --;
#	print "common prefix = ",$common_prefix,"\n";
    } else {
	next if m/^\Q$common_prefix\E/;
#	print "index cutoff = ",$index_cutoff,"\n\n";
	last;
    }
}
close(READS);

$index_cutoff = shift @ARGV if $#ARGV > -1;
print STDERR "index cutoff = ",$index_cutoff,"\n\n";

my ($edgeIDX, $len, $multiplicity) = (-1,-1,-1) ;
my ($readIDX, $color) = (-1,"none") ;

my @edgeColor, my @edgeMultiplicity;
my $reference = "1 1 1"; # red, covered only by first set of reads
my $contigs   = "0.33 1 1"; # green, covered by only second set of reads
my $both      = "0.66 1 1"; # blue, covered by both sets of reads
my $neither   = "0 0 0"; # black, covered by no reads (multiplicity = 0) ???
die("$interval_file is not a .intv file.\n") if not ($interval_file =~ m/\.intv/);
open(INTV,"$interval_file") or die("cannot open intervals file $interval_file\n");
while(<INTV>) {
    my @l = split;
    if (m/^EDGE/) {
	$edgeColor[$edgeIDX] = $color  if $edgeIDX > -1;
	($edgeIDX, $len, $multiplicity) = ($l[1],$l[3],$l[5]);
	$color = $neither;
	$edgeMultiplicity[$edgeIDX] = $multiplicity;
    } else {
	$readIDX = $l[1];
	if (($readIDX / 2) < $index_cutoff) {
	    $color = $reference if $color eq $neither;
	    $color = $both if $color eq $contigs;
	} else {
	    $color = $contigs if $color eq $neither;
	    $color = $both if $color eq $reference;
	}
    }
}
close(INTV);
$edgeColor[$edgeIDX] = $color; # for the last edge

die("$gvz_file is not a GraphViz file.\n") if not ($gvz_file =~ m/\.(gvz|dot)/);
open(GVZ,"$gvz_file") or die("cannot open GraphViz file $gvz_file\n");
while(<GVZ>) {
    my @l = split;
    print and next if ($#l < 1) or (not $l[1] eq "->");
    my ($foo,$edgeID) = split('"',$l[3]);
    next if $showZeroMultiplicityEdges == 0 and $edgeMultiplicity[$edgeID] == 0;
    my ($multiplicity,$bar) = split('"',$l[5]);
    $l[5] = "m=" . $edgeMultiplicity[$edgeID] . "\" , color=\"" . $edgeColor[$edgeID] . "\"]" ;
    print "\t",join(" ",@l), "\n";
#    die( "gvz mult = $1 but contig mult = " . $edgeMultiplicity[$edgeID] ) if $1 != $edgeMultiplicity[$edgeID];
}
close(GVZ);
