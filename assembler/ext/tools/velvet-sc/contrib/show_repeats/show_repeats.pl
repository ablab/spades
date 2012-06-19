#!/usr/bin/perl -w

#    show_repeats.pl             
#
#    Show the distribution of repeated contigs from a velvet assembly. The repeats
#    are inferred from the coverage distribution as presented in stats.txt
#
#    Copyright (C) 2009 Ken Doig, Torsten Seemann
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

use strict;
use List::Util qw(max);

my(@Options, $minlen, $mincov, $width, $binsize, $repcut, $csvout);
setOptions();

# cmdline parameter sanity

$width = $width || $ENV{'COLUMNS'} || 80;
$width -= 20; # for labels

# MAIN

my $cov    = read_stats(\*ARGV);
my $expcov = max_cov($cov);

if ( ! $csvout )
{
  print_histogram($cov,$expcov);
  print "Predicted expected coverage: $expcov\n";
  print "velvetg parameters: -exp_cov $expcov -cov_cutoff $mincov\n";
}
else
{
  print_csv($cov,$expcov);
}

# END

#
# Read in a Velvet stats.txt file
#
# Velvet stats.txt format
# column 0   1     2    3  4         5           6            7           8
# value  ID  lgth  out  in long_cov  short1_cov  short1_Ocov  short2_cov  short2_Ocov
#

sub read_stats
{
  my($fh) = @_;
  my %cov;

  my $line = <$fh> || die $!;	# check stats.txt file is there

  while (<$fh>)
  {
    chomp;
    my @x = split m/\t/;
    next unless @x >= 7;
    next unless $x[1] =~ m/^\d+$/;
    next unless $x[1] >= $minlen;
    my $cov = $x[5] + $x[7];    # add both short read channels
    next unless $cov > $mincov;
    $cov{int($cov)} += $x[1];   # sum the contig lengths for each integer coverage
  }

  return \%cov;
}

#
# Find the coverage with the greatest length of contigs
#

sub max_cov
{
  my($hash_ref) = @_;
  my %covh = %{$hash_ref};
  my $maxcov = 0;
  my $maxlen = 0;

  for my $key (keys %covh)
  {
    if ( $covh{$key} > $maxlen )
    {
      $maxlen = $covh{$key};
      $maxcov = $key;
    }
  }
  return $maxcov; 
}

#
# Print a bar graph of the contig lengths for binned coverages
#

sub print_histogram
{
  my($hash_ref,$expcov) = @_;
  my %covh;

  #
  # Sum the contig lengths for a range of coverages
  #
  for my $key (keys %{$hash_ref})
  {
    next if ( $key < $expcov * $repcut );    # skip contigs with insufficient coverage

    $covh{int($key/($expcov*$binsize))} += ${$hash_ref}{$key};
  }

  #
  # Print an ascii graph line for each range of coverages
  #
  my $maxv = max values %covh;
  my $maxc = max keys   %covh;

  for (my $c=0; $c <= $maxc; ++$c)
  {
    my $f = $covh{$c};
       $f = 0 if ( ! defined $f );
    my $N = int( $f  * $width / $maxv );
    my $coverage = sprintf "X %d", $c*$binsize;
    printf "%8s | %6d | %-${width}s\n", $coverage, $f, ('*'x$N);
  }
}

#
# Output csv of the coverage and contig lengths for all contigs
#

sub print_csv
{
  my($hash_ref,$expcov) = @_;

  print "Coverage,Contig_lengths\n";

  for my $key (sort { $a <=> $b } keys %{$hash_ref})
  {
    next if ( $key < $repcut * $expcov);
    print $key,",",${$hash_ref}{$key},"\n";
  }  
}

#----------------------------------------------------------------------
# Option setting routines

sub setOptions
{
  use Getopt::Long;

  @Options = (
    {OPT=>"help",        VAR=>\&usage,                  DESC=>"This help"},
    {OPT=>"l|minlen=i",  VAR=>\$minlen,  DEFAULT=>50,   DESC=>"Minimum size contigs to include"},
    {OPT=>"c|mincov=i",  VAR=>\$mincov,  DEFAULT=>0,    DESC=>"Minimum kmer coverage to include"},
    {OPT=>"b|binsize=i", VAR=>\$binsize, DEFAULT=>1,    DESC=>"Bin size for coverage"},
    {OPT=>"r|repeat=f",  VAR=>\$repcut,  DEFAULT=>1.5,  DESC=>"Minimum X coverage for repeats"},
    {OPT=>"w|width=i",   VAR=>\$width,   DEFAULT=>0,    DESC=>"Width of output graph (0=auto)"},
    {OPT=>"x|xls!",      VAR=>\$csvout,  DEFAULT=>0,    DESC=>"CSV output"}
  );

  (!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options)
  {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}}))
    {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage
{
  print "Usage: $0 [options] stats.txt\n";
  
  foreach (@Options)
  {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
