#! /usr/bin/perl

#***************************************************************************
# Title:          errcorr-lapper-pair-mem.pl
# Author:         Haixu Tang
# Created:        Jun. 2002
# Last modified:  May. 2004
#
# Copyright (c) 2001-2004 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
#***************************************************************************/

$MYHOME = $ENV{'EULERBIN'};

my $allsteps = 4;
my $min_id = 0.95;
my $overlen = 20;
my $qualfile;
my $qual;
my $num;

$num = @ARGV;
if($num > 3)	{
	$overlen = $ARGV[2];
	$min_id = $ARGV[3];
	$qualfile = $ARGV[1];
	$qual = 1;
} else {
	if($num > 2)	{
		$overlen = $ARGV[1];
		$min_id = $ARGV[2];
	}
	$qual = 0;
}

if($allsteps < 1)	{
	print "No error correction: check iterations\n";
	die;
}

$alnfile=$ARGV[0].".aln";

if($qual == 1)	{
	system("$MYHOME/overlapper-all -i $ARGV[0] -q $qualfile -o $alnfile -l $overlen -H ");
	system("$MYHOME/errcorr_pair_mem -i $alnfile -s $ARGV[0] -q $qualfile -o reads.tmp -c 1 -k 5 -H");
} else	{
	system("$MYHOME/overlapper-all -i $ARGV[0] -o $alnfile -l $overlen -H ");
	system("$MYHOME/errcorr_pair_mem -i $alnfile -s $ARGV[0] -o reads.tmp -c 1 -k 5 -H ");
}
$i = 1;
while($i <= $allsteps)	{
	system("mv reads.tmp reads.txt");
	system("mv reads.tmp.score reads.txt.score");
	system("$MYHOME/overlapper-all -i reads.txt -q reads.txt.score -o reads.aln -l $overlen -p 200 -d $min_id -H");
	$i = $i + 1;
	system("$MYHOME/errcorr_pair_mem -i reads.aln -s reads.txt -q reads.txt.score -o reads.tmp -k 5 -H ");
}

system("mv reads.tmp reads.txt");
system("mv reads.tmp.score reads.txt.score");
system("$MYHOME/overlapper-all -i reads.txt -q reads.txt.score -o reads.aln -l $overlen -d $min_id -u 10 -H");

system("mv reads.txt $ARGV[0].new");
system("mv reads.txt.score $ARGV[0].new.score");
system("mv reads.aln $ARGV[0].new.aln");
