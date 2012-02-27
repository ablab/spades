#!/usr/bin/env perl

if ($#ARGV != 0) {
  print "usage: addheader.pl source\n";
  exit(0);
}
$source = $ARGV[0];

$pid = ??;
if (! -e $source) {
  print "$source does not exist\n";
  exit(0);
}

if ($source =~ /\/([^\/]+)$/) {
		$base = $1;
}
else {
		$base = $source;
}

`grep -q "Title:" $source`;
if ($? == 0) {
		exit(0);
}

open(TMP, ">$pid.header.tmp") or die "cannot open tmp file \n";
print TMP <<END_HEADER;
/***************************************************************************
 * Title:          $base 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
END_HEADER
`cat $pid.header.tmp $source > $pid.tmp`;
`mv $pid.tmp $source`;
`rm $pid.header.tmp`;
