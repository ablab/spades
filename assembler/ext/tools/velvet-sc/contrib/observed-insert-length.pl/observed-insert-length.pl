#!/usr/bin/perl -w

#    observed-insert-length.pl             
#
#    Displays the insert lengths of the read pairs which happen to be on the same
#    contig of a velvetg Graph2 file.
#
#    Copyright (C) 2009 Torsten Seeman 
#    Modified by Daniel Zerbino,  August 17, 2009
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
################################################################################
use strict;
use List::Util qw(max);

my(@Options, $verbose, $library, $width);
setOptions();

# cmdline parameter sanity

$width = $width || $ENV{'COLUMNS'} || 80;
$width -= 20; # for labels
$library = $library || 1;
$library = $library * 2 - 1;

# MAIN

my $insert_lengths = read_stats($ARGV[0]);
print_histogram($insert_lengths);
my $ins_length_mode = math_mode($insert_lengths);
my $ins_length_median = math_median($insert_lengths);
my $ins_length_sd = math_sd($insert_lengths);
print "Observed median insert length: $ins_length_median\n";
print "Observed mode of insert length: $ins_length_mode\n";
print "Observed sample standard deviation: $ins_length_sd\n"; 
if ($library == 1) {
	$library = "";
}
if ($ins_length_sd == 0) {
	$ins_length_sd = 1;
}	
print "Suggested velvetg parameters: -ins_length$library $ins_length_median -ins_length${library}_sd $ins_length_sd\n";

# END

sub read_stats {
  my $directory = shift; 

  die "Directory $directory does not exist, please check your parameters.\n" unless (-e $directory);
  die "Sequence file $directory/Sequences does not exist, please check that you ran Velvetg correctly.\n" unless (-e "$directory/Sequences");
  die "Graph file $directory/Graph2 does not exist, please re-run velvetg with the -read_trkg parameter on.\n" unless (-e "$directory/Graph2");

  open IN, "$directory/Sequences"; 
  my %read_pair = ();
  my $unmatched_read = 0;

  while(<IN>) {
    chomp;
    my @x = split m/\t/;
    next unless @x > 1;
    my $read_no = $x[@x - 2];
    my $read_cat = $x[@x - 1];
    next unless $read_cat == $library;
    if ($unmatched_read == 0) {
      $unmatched_read = $read_no;
    } else {
      $read_pair{ $read_no } = $unmatched_read;
      $read_pair{ $unmatched_read } = $read_no;
      $unmatched_read = 0;
    }
  }

  die "No paired reads found for library ".($library/2 + 0.5)."\n" unless scalar(keys %read_pair) > 0;

  close IN;

  open IN, "$directory/Graph2";
  $_ = <IN>;
  chomp;
  my @line = split m/\t/;
  my $kmer = $line[2];
  my @insert_lengths;
  my %read_position = ();
  my %node_lengths = ();
  my $node = 0;
  while (<IN>) {
    chomp;
    my @x = split m/\t/;
    next unless defined $x[0];
    if ($x[0] =~ m/NODE/) {
	$node_lengths{ $x[1] } = $x[2];
	next;
    } elsif ($x[0] =~ m/NR/) {
	$node = int($x[1]);
	next;
    }

    next unless $node != 0;

    if ($node < 0) {
	my $read_no = $x[0];
	my $read_pos = $x[1];
	my $read_off = $x[2];
	if ($read_pos > 0) {$read_position{ $read_no } = [-$node, $node_lengths{ -$node } - $read_pos + $read_off];}
    } else {
      my $read_no = $x[0];
      my $read_pos = $x[1];
      my $read_off = $x[2];
	
      if ($read_pos > 0 && $read_pair{$read_no} && $read_position{$read_pair{$read_no}} && $read_position{ $read_pair{ $read_no }}->[0] == $node) {
	my $length = abs($read_pos - $read_off - $read_position{ $read_pair{ $read_no} }->[1]); 	
	push @insert_lengths, $length + $kmer - 1;
      }
    }
  }
  close IN;
  return \@insert_lengths;
}

# Mathematical mode is the most frequent observation (peak of histogram)

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

# Mathematical median is the middle obersvation (center of the histogram)

sub math_median {
  my($array_ref) = @_;
  my @sorted = sort { $a <=> $b } @{$array_ref};
  return $sorted[int(@sorted/2)]; 
}

# Standard error of dataset

sub math_sd {
  my($array_ref) = @_;
  my $mean = 0;
  ($mean += $_) for @{$array_ref};
  $mean /= @{$array_ref};
  my $total = 0;
  ($total += ($_ - $mean)**2) for @{$array_ref};
  return sqrt ($total / @{$array_ref});
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
    {OPT=>"w|width=i",  VAR=>\$width, DEFAULT=>0, DESC=>"Width of output graph (0=auto)"},
    {OPT=>"l|library=i", VAR=>\$library, DEFAULT=>1, DESC=>"Read library as given to Velvet"}
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
  print "Usage: $0 [options] Velvet_directory\n";
  foreach (@Options) {
    printf "  --%-13s %s%s.\n",$_->{OPT},$_->{DESC},
           defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
  }
  exit(1);
}
 
#----------------------------------------------------------------------
