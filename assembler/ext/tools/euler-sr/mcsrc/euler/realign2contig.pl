#!/usr/bin/perl -w

#***************************************************************************
# Title:          realign2contig.pl
# Author:         Vagisha Sharma
# Created:        Jun. 2002
# Last modified:  May. 2004
#
# Copyright (c) 2001-2004 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
#**************************************************************************/

use strict;

my $mafile = shift;
my $contigfile = shift;

my $output = "";

open (MA, "$mafile") or die "could not open $mafile\n";

my $header = "";
my $outline;
my $i;
my $seq = "";
my $gapseq = "";
my $len = 0;
while (<MA>) {
  my $line = $_;
  chomp ($line);

  if ($line =~ m/^% Contig/) {
    my @tokens = split (/\s+/, $line);
    my $connum = $tokens[2];
    #print "Contig num is : $connum\n";
    if ($header ne "") {
      my $temp = $seq;
      $temp =~ s/\n//g;
      $len = length ($temp);
# Modified by H. Tang
      $seq =~s/\n//g;
      $len = length($seq);
      $output .= "$header $len\n";
      $i = 0;
      while($i < $len)        {
        $outline = substr($seq, $i, 60);
        $output .= "$outline\n";
        $i += 60; 
      }
    }
    my $connum1 = $connum+1;
    $header = ">Contig$connum1";
    $seq = "";
    $gapseq = "";
  }
  elsif ($line =~ m/^-/ && $line =~ m/-$/ && $line !~ m/a|A|c|C|g|G|t|T/) {
    #print "$line\n";
    $line = <MA>;
    last  if (eof (MA));
    my $thisline = $line;
    #print "$thisline\n";
    $gapseq .= $thisline if ($line !~ m/^\s|%/);
    $thisline =~ s/-//g;
    #print "$thisline\n";

    #$thisline =~ tr/AGCT/agct/;
    $seq .= $thisline if ($line !~ m/^\s|%/);
  }

}

close MA;

# Modified by H. Tang
$seq =~s/\n//g;
$len = length($seq);
$output .= "$header $len\n";
$i = 0;
while($i < $len)        {
        $outline = substr($seq, $i, 60);
        $output .= "$outline\n";
        $i += 60; 
}

#output sequence

open (OUT, ">$contigfile") or die "could not open $contigfile\n";
print OUT "$output";
close OUT;

exit;
