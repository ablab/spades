#!/usr/bin/perl -w

#    velvet-estimate-exp_cov.pl             
#
#    Estimates the expected k-mer coverage parameter (-exp_cov) for velvetg
#    by finding the mode of coverage distribution as presented in stats.txt
#
#    Copyright (C) 2009 Torsten Seemann
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

my(@Options, $verbose, $minlen, $mincov, $width);
setOptions();

# cmdline parameter sanity

$width = $width || $ENV{'COLUMNS'} || 80;
$width -= 20; # for labels

# MAIN

my $cov = read_stats(\*ARGV);
print_histogram($cov);
my $expcov = math_mode($cov);
print "Predicted expected coverage: $expcov\n";
print "velvetg parameters: -exp_cov $expcov -cov_cutoff $mincov\n";

# END

# Velvet stats.txt format
#0       1       2       3       4               5               6               7               8
#ID      lgth    out     in      long_cov        short1_cov      short1_Ocov     short2_cov      short2_Ocov

sub read_stats {
  my($fh) = @_;
  my @cov;
  while (<$fh>) {
    chomp;
    my @x = split m/\t/;
    next unless @x >= 7;
    next unless $x[1] =~ m/^\d+$/;
    next unless $x[1] >= $minlen;
  #  my $cov = $x[6] + $x[8]; # add both short read channels
    my $cov = $x[5] + $x[7]; # add both short read channels
    next unless $cov > $mincov;
    push @cov, int($cov);
  }
  return \@cov;
}

# Mathematical "mode" is the most frequent observation (peak of histogram)

sub math_mode {
  my($array_ref) = @_;
  my %freq_of;
  map { $freq_of{$_}++ } @{$array_ref};
  my $mode=0;
  $freq_of{$mode} = 0; # sentinel
  for my $x (keys %freq_of) {
    $mode = $x if $freq_of{$x} > $freq_of{$mode};
  }
  return $mode; 
}

# Print a bar graph of an array of values

sub print_histogram {
  my($cov) = @_;
  my %freq_of;
  for my $c (@$cov) { 
    $freq_of{$c}++ 
  }
  my $max = max(values %freq_of);
  for my $c (sort { $a <=> $b } keys %freq_of) {
    my $f = $freq_of{$c};
    my $N = int( $f  * $width / $max );
    printf "%6d | %6d | %-${width}s\n", $c, $f, ('*'x$N);
  }
}


#----------------------------------------------------------------------
# Option setting routines

sub setOptions {
  use Getopt::Long;

  @Options = (
    {OPT=>"help",    VAR=>\&usage,             DESC=>"This help"},
    {OPT=>"v|verbose!",  VAR=>\$verbose, DEFAULT=>0, DESC=>"Verbose output"},
    {OPT=>"l|minlen=i",  VAR=>\$minlen, DEFAULT=>61, DESC=>"Minimum size contigs to include"},
    {OPT=>"c|mincov=i",  VAR=>\$mincov, DEFAULT=>0, DESC=>"Minimum kmer coverage to include"},
    {OPT=>"w|width=i",  VAR=>\$width, DEFAULT=>0, DESC=>"Width of output graph (0=auto)"},
  );

  (!@ARGV) && (usage());

  &GetOptions(map {$_->{OPT}, $_->{VAR}} @Options) || usage();

  # Now setup default values.
  foreach (@Options) {
    if (defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
      ${$_->{VAR}} = $_->{DEFAULT};
    }
  }
}

sub usage {
  print "Usage: $0 [options] stats.txt\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
