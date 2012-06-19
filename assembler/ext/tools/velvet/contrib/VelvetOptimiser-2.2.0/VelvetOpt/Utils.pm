#
#       VelvetOpt::Utils.pm
#
#       Copyright 2008,2009,2010 Simon Gladman <simon.gladman@csiro.au>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

#		Version 2.1.3

#	Changes for Version 2.0.1
#	Added Mikael Brandstrom Durling's numCpus and freeMem for the Mac.
#
#	Changes for Version 2.1.0
#	Fixed bug in estExpCov so it now correctly uses all short read categories not just the first two.
#
#	Changes for Version 2.1.2
#	Fixed bug in estExpCov so it now won't take columns with "N/A" or "Inf" into account
#
#	Changes for Version 2.1.3
#	Changed the minimum contig size to use for estimating expected coverage to 3*kmer size -1 and set the minimum coverage to 2 instead of 0.
#	This should get rid of exp_covs of 1 when it should be very high for assembling reads that weren't ampped to a reference using one of the standard read mapping programs


package VelvetOpt::Utils;

use strict;
use warnings;
use POSIX qw(ceil floor);
use Carp;
use List::Util qw(max);
use Bio::SeqIO;

# 	num_cpu
#	It returns the number of cpus present in the system if linux.
#	If it is MAC then it returns the number of cores present.
#	If the OS is not linux or Mac then it returns 1.
#	Written by Torsten Seemann 2009 (linux) and Mikael Brandstrom Durling 2009 (Mac).

sub num_cpu {
    if ( $^O =~ m/linux/i ) {
        my ($num) = qx(grep -c ^processor /proc/cpuinfo);
        chomp $num;
        return $num if $num =~ m/^\d+/;
    }
	elsif( $^O =~ m/darwin/i){
		my ($num) = qx(system_profiler SPHardwareDataType | grep Cores);
		$num =~ /.*Cores: (\d+)/;
		$num =$1;
		return $num;
	}
    return 1;
}

#	free_mem
#	Returns the current amount of free memory
#	Mac Section written by Mikael Brandstrom Durling 2009 (Mac).

sub free_mem {
	if( $^O =~ m/linux/i ) {
		my $x       = `free | grep '^Mem:' | sed 's/  */~/g' | cut -d '~' -f 4,7`;
		my @tmp     = split "~", $x;
		my $total   = $tmp[0] + $tmp[1];
		my $totalGB = $total / 1024 / 1024;
		return $totalGB;
	}
	elsif( $^O =~ m/darwin/i){
		my ($tmp) = qx(vm_stat | grep size);
		$tmp =~ /.*size of (\d+) bytes.*/;
		my $page_size = $1;
		($tmp) = qx(vm_stat | grep free);
		$tmp =~ /[^0-9]+(\d+).*/;
		my $free_pages = $1;
		my $totalGB = ($free_pages * $page_size) / 1024 / 1024 / 1024;
		return $totalGB;
	}
}

#	estExpCov
#   it returns the expected coverage of short reads from an assembly by
#   performing a math mode on the stats.txt file supplied..  It looks at
#	all the short_cov? columns..  Uses minimum contig length and minimum coverage.
#	needs the stats.txt file path and name, and the k-value used in the assembly.
#	Original algorithm by Torsten Seemann 2009 under the GPL.
#	Adapted by Simon Gladman 2009.
#	It does a weighted mode...

sub estExpCov {
    use List::Util qw(max);
    my $file   = shift;
    my $kmer   = shift;
    my $minlen = 3 * $kmer - 1;
    my $mincov = 2;
    my $fh;
    unless ( open IN, $file ) {
        croak "Unable to open $file for exp_cov determination.\n";
    }
    my @cov;
    while (<IN>) {
        chomp;
        my @x = split m/\t/;
		my $len = scalar @x;
        next unless @x >= 7;
        next unless $x[1] =~ m/^\d+$/;
        next unless $x[1] >= $minlen;
		
		#add all the short_cov columns..
		my $cov = 0;
		for(my $i = 5; $i < $len; $i += 2){
			if($x[$i] =~ /\d/){
				$cov += $x[$i];
			}
		}
        next unless $cov > $mincov;
        push @cov, ( ( int($cov) ) x $x[1] );
    }

    my %freq_of;
    map { $freq_of{$_}++ } @cov;
    my $mode = 0;
    $freq_of{$mode} = 0;            # sentinel
    for my $x ( keys %freq_of ) {
        $mode = $x if $freq_of{$x} > $freq_of{$mode};
    }
    return $mode;
}

#	estVelvetMemUse
#	returns the estimated memory usage for velvet in GB

sub estVelvetMemUse {
	my ($readsize, $genomesize, $numreads, $k) = @_;
	my $velvetgmem = -109635 + 18977*$readsize + 86326*$genomesize + 233353*$numreads - 51092*$k;
	my $out = ($velvetgmem/1024) / 1024;
	return $out;
}

#	getReadSizeNum
#	returns the number of reads and average size in the short and shortPaired categories...

sub getReadSizeNum {
	my $f = shift;
	my %reads;
	my $num = 0;
	my $currentfiletype = "fasta";
	#first pull apart the velveth string and get the short and shortpaired filenames...
	my @l = split /\s+/, $f;
	my $i = 0;
	foreach (@l){
		if(/^-/){
			if(/^-fasta/){
				$currentfiletype = "fasta";
			}
			elsif(/^-fastq/){
				$currentfiletype = "fastq";
			}
			elsif(/(-eland)|(-gerald)|(-fasta.gz)|(-fastq.gz)/) {
				croak "Cannot estimate memory usage from file types other than fasta or fastq..\n";
			}
		}
		elsif(-r $_){
			my $file = $_;
			if($currentfiletype eq "fasta"){
				my $x = `grep -c "^>" $file`;
				chomp($x);
				$num += $x;
				my $l = &getReadLength($file, 'Fasta');
				$reads{$l} += $x;
				print STDERR "File: $file has $x reads of length $l\n";
			}
			else {
				my $x = `grep -c "^@" $file`;
				chomp($x);
				$num += $x;
				my $l = &getReadLength($file, 'Fastq');
				$reads{$l} += $x;
				print STDERR "File: $file has $x reads of length $l\n";
			}
		}
		$i ++;
	}
	my $totlength = 0;
	foreach my $k (keys %reads){
		$totlength += ($reads{$k} * $k);
	}
	
	
	my @results;
	push @results, floor($totlength/$num);
	push @results, ($num/1000000);
	printf STDERR "Total reads: %.1f million. Avg length: %.1f\n",($num/1000000), ($totlength/$num);
	return @results;
}

# getReadLength - returns the length of the first read in a file of type fasta or fastq..
#
sub getReadLength {
	my ($f, $t) = @_;
	my $sio = Bio::SeqIO->new(-file => $f, -format => $t);
	my $s = $sio->next_seq() or croak "Something went bad while reading file $f!\n";
	return $s->length;
}

return 1;

