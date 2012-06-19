#!/usr/bin/perl -w
#
#       asmbly_splitter2.pl
#
#       Copyright 2008 Simon Gladman <simon.gladman@csiro.au>
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
#
#       A script to split out a single contig's data from an AFG file
#       produced by Daniel Zerbino's Velvet assembler.  The output is
#       in AFG format and is a subset of the data from the original file.
#       Because of the size of the AFG file, this script can take quite a
#       while to run.
#
#       Usage:  ./asmbly_splitter2.pl <contig number> <afg file>
#
#
#       Where:  <contig number> is the number of the contig of interest.
#               <afg file> is the AFG file produced by Velvet.
#
#####################################################################
#
# Modified by Daniel Zerbino to allow for scaffolding and to preserve read
# information
use strict;

my $usage = "Usage: $0 <contig number> <afg file>\n";

sub toggle {
    my $x = shift;
    if($x == 1){
        $x = 0;
    }
    else {
        $x = 1;
    }
    return $x;
}

my $contig = $ARGV[0];
my $file = $ARGV[1];

unless($contig){die "Usage: $usage \n";}

unless($file){die "Usage: $usage \n";}

unless(-e $file){die "$0: File $file doesn't exist.\n"};

my $outfile = $file;
$outfile =~ s/\.afg$/_$contig\.afg/;

open IN, $file;
open OUT, ">$outfile";

#need to find the contig first and store the read data..
my @reads;
my $readcount = 0;
my $count = 0;
my $line;
my $iid_line;
my $contigString;

print STDERR "Searching for start of contig $contig\n";

#first pass of the file.  Can't store reads..  too much info. Therefore,
#need two passes of the file...
#first pass gets the contig information.

$_ = <IN>;
while($_){
    $count++;
    if($count % 1000000 == 0){
        print STDERR ".";
    }

    if($_ =~ /^\{CTG/){
        #found the contigs
        $iid_line = <IN>;
        $line = <IN>;
	if (!$line) {
		print STDERR $_;
		print STDERR $iid_line;
		exit;
	}	
        if($line =~ /^eid:$contig(-[0-9]*)?$/){
            #found the right contig..
            print STDERR "\nFound a part of contig $contig\n";
	    $contigString .= "{CTG\n";
	    $contigString .= $iid_line;
            $contigString .= $line;
            while(<IN>){
                if($_ =~ /^\{CTG/ || $_ =~ /^\{SCF/){
                    #found the end of the contig..
                    last;
                }
                else {
                    $contigString .= $_;
                    if(/^src:/){
                        #this is a read number, store it..
                        chomp;
                        my @tmp = split /:/, $_;
                        push @reads, $tmp[1];
                        $readcount ++;
		    }
                }
            }
        } else {
		$_ = <IN>;
	}
    }
    elsif($_ =~ /\{SCF/){
        #found a scaffold 
        $line = <IN>;
        if($line =~ /^eid:$contig$/){
            #found the right scaffold..
	    $contigString .= "{SCF\n";
            $contigString .= $line;
            while(<IN>){
                if($_ =~ /^\{CTG/ || $_ =~ /^\{SCF/){
                    #found the end of the scaffold..
                    last;
                }
                else {
                    $contigString .= $_;
                }
            }
	    last;
        } else {
            $_ = <IN>;
	}
    }
    else {
	$_ = <IN>;
    }
}

close IN;

print STDERR "Found the end of contig $contig\nNumber of reads = $readcount\nStarting to sort reads\n";

#remove double entries in read list
my @uniq = keys %{{ map { $_ => 1 } @reads }};
#sort the reads.
my @sorted = sort { $a <=> $b } @uniq;
@reads = @sorted;
my @reads_copy = @reads;
print STDERR "Finished sorting reads\nStarting to look for fragments\n";

# Second pass through the file to print out all the appropriate
# library, fragment, read, contig and scaffold info

open IN, $file;
my $currentread = shift @reads;
my $currentread_copy = shift @reads_copy;
my $foundreads = 0;
my @foundfrags = ();  

print STDERR "Reads to find = " . ($readcount - $foundreads) . "\tCurrent read: $currentread\n";

while(<IN>) {
	# All library info is printed out
	if(/\{LIB/){
		print OUT $_;
		while(<IN>){
			print OUT $_;
			if(/\}/){
				last;
			}
		}
		$_ = <IN>;
		print OUT $_;

	}

	# All fragment blocks that contain two selected reads are printed out
	elsif(/\{FRG/){
		my $frgString = $_;
		$frgString .= <IN>;

		$_ = <IN>;
		chomp;
		my @tmp = split /[:,]/, $_;
		my $read = $tmp[1];
		my $read2 = $tmp[2];

		while (defined($currentread_copy) && $currentread_copy < $read) {
			$currentread_copy = shift @reads_copy;
		}

		if (!defined($currentread_copy) || $currentread_copy > $read) {
			next;
		} else {
			$currentread_copy = shift @reads_copy;
		}

		if (!defined($currentread_copy) || $read2 < $currentread_copy) {
			next;
		}

		# If we got to this line then both read IDs were recognized!
		push @foundfrags, $read;
		push @foundfrags, $read2;
	
		print OUT $frgString;
		print OUT "$_\n";
		while(<IN>){
			print OUT $_;
			if(/\}/){
				last;
			}
		}
	}

	# All selected read blocks are printed out
	elsif(/\{RED/){
		my $line = <IN>;
		if(defined $currentread && $line =~ /^iid:$currentread/){
			$foundreads ++;
			print OUT "\{RED\n";
			print OUT $line;
			while(<IN>){
				if (/frg:/) {
					if (defined $foundfrags[0] && $currentread != $foundfrags[0]) {
						next;
					} else {
						shift @foundfrags;
					}
				}
				print OUT $_;
				if(/\}/){
					last;
				}
			}
			if (scalar( @reads) != 0) {
			  $currentread = shift @reads;

		
			  if($foundreads % 100 == 0){
			    print STDERR "Reads to find = " . ($readcount - $foundreads) . "\tCurrent read: $currentread\n";
			  }
			  if($foundreads == $readcount){
			    last;
			  }
		}
			else {
			  print STDERR 'No reads left in RED!';
			  next;
			}
			      

	}
}
      }
print "\n";


print OUT $contigString;
close OUT;

print STDERR "Finished!\n";
