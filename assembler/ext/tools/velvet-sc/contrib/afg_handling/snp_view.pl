#!/usr/bin/perl -w
#
#   snp_view.pl
#       
#   Copyright 2008 Simon Gladman <gla048@localhost.localdomain>
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#   
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#      
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#   MA 02110-1301, USA.
#
#	A script to take a single contig AFG file and a position in that
#	contig and show the reads aligned with the consensus sequence 
#	centred on that position.  It outputs a text file called "snp_view.txt"
#
#	Usage: ./snp_view.pl <afg file> <position>
#
#	Where: 	<afg file> is a single contig AFG file extracted from the 
#			velvet_assy.afg file (produced by Velvet) using asmbly_splitter.pl
#			<position> is the position in the contig about which the view
#			will be centred.  The position must be >= a full read length from
#			start of the contig and must be <= a full read length from
#			the end of the contig.
#
#	Other parameters:
#			$readlength: by changing the size of this parameter, it changes
#			the width of the view.  Note that the minimum value is the read length.
#
#			$outfile: set this parameter to whatever output filename you want.
#
# Modified by Daniel Zerbino (Aug 28 2009) to display reverse strand reads in the 
# snp_view_reads file and to better handle the display of snps near the end of contigs

use strict;
use List::Util qw[min max];

my $outfile = "snp_view.txt";
my $readlength = 50;

my $usage = "snp_view.pl\nCopyright 2008 Simon Gladman - CSIRO\nUsage: $0 <single_contig_afg_file> <integer_position_in_contig>\n";

my $afg_file = $ARGV[0];
my $posn = $ARGV[1];

unless($afg_file) { die $usage; }
unless($posn) {die $usage; }

my $low = $posn - $readlength;
my $high = $posn + $readlength;
my $con;
my $seq = "";
my $incontig = 0;
my %offs;
my %seqs;
my %revs;
my $off;
my $src;

unless(open IN, $afg_file){die "No such $afg_file\n$!\n";}
unless(open OUT, ">$outfile"){die "No such $afg_file\n";}
unless(open OUT2, ">snp_view_ref.fa"){die "Couldn't open fasta outfile\n";}
unless(open OUT3, ">snp_view_reads.fa"){die "Couldn't open reads fasta file\n";}

while(<IN>){
	if(/{CTG/){
		$con = <IN>;
		chomp($con);
		$con =~ s/^iid://;
		$incontig = 1;
	}
	if(/^seq/ && $incontig){
		while(<IN>){
			chomp;
			if(/\./){
				last;
			}
			else {
				$seq .= $_;
			}
		}
	}
	if(/{TLE/){
        $_ = <IN>;
        chomp;
        my @temp = split /:/, $_;
        $src = $temp[1];
        $_ = <IN>;
        chomp;
        @temp = split /:/, $_;
        $off = $temp[1];
        if($off >= $low && $off <= $posn){
            #we want this one..
            #print "$low, $off, $high\n";
            $offs{$src} = $off;
        }
        $_ = <IN>;
        chomp;
        if(/clr:[0-9]+,0/){
            $revs{$src} = "r";
        }
        if(/clr:0,[0-9]+/){
            $revs{$src} = "f";
        }
    }
}

print OUT "Contig $con, position of X below is: $posn\n";
my $wantSeqStart = max($low, 0);
my $wantSeqFinish = min($high, length $seq);
my $wantSeq = substr($seq, $wantSeqStart, $wantSeqFinish - $wantSeqStart);
print OUT "Position:\t\t$low";
my $num_spaces = 2 * $readlength - length($low) -1;
for(my $i = 0; $i < $num_spaces; $i ++){
	if($i == $readlength-length($low)){
		print OUT "X";
	} else {
		print OUT " ";
	}
}
print OUT "$high\n";
print OUT "Read #\t\tdir\t";
if ($low < 0) { 
	for(my $i = 0; $i > $low; $i--){
	print OUT ".";
	}
}
print OUT "$wantSeq\n";
print OUT2 ">Reference\n";
print OUT2 "$wantSeq\n";
close IN;

unless(open IN, $afg_file){die "No such $afg_file\n$!\n";}

while(<IN>){
    if(/{RED/){
        $_ = <IN>;
        chomp;
        my @temp = split /:/, $_;
        if(defined $offs{$temp[1]}){
            $_ = <IN>;
            $_ = <IN>;
            $_ = <IN>;
            chomp;
            $seqs{$temp[1]} = $_;
        }
    }
    if(/{CTG/){
    	last;
	}
}

foreach my $key (sort { $offs{$a} <=> $offs{$b} } keys %offs){
    my $dots = $offs{$key} - $low;
    print OUT "$key\t";
    if($key < 10000000){ print OUT "\t";}
    print OUT $revs{$key} ."\t";
    for(my $i = 0; $i < $dots; $i ++){
        print OUT ".";
    }
    if($revs{$key} eq "f"){
        print OUT $seqs{$key} . "\n";
        print OUT3 ">$key\n";
        print OUT3 $seqs{$key} . "\n";
    }
    else {
        my $rev = reverse($seqs{$key});
        $rev =~ tr /ATCG/TAGC/;
        print OUT "$rev\n";
        print OUT3 ">$key\n";
        print OUT3 $seqs{$key} . "\n";
    }
}
