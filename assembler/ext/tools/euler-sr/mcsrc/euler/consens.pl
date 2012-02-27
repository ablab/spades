#! /usr/bin/perl

#***************************************************************************
# Title:          consens.pl
# Author:         Vagisha Sharma
# Created:        Jun. 2002
# Last modified:  May. 2004
#
# Copyright (c) 2001-2004 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
#**************************************************************************/

my $RAHOME=$ENV{'EULERBIN'};

`$RAHOME/create_alignment.pl $ARGV[0] $ARGV[1] 80 3 realigner.inp`;
`$RAHOME/realigner/converter -hV < realigner.inp| $RAHOME/realigner/ReAligner | $RAHOME/realigner/converter -vHC > realigner.out`;
`$RAHOME/realign2contig.pl realigner.out $ARGV[2]`;
