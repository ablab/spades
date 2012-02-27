#!/usr/bin/env perl

############################################################################
# Title:          SummarizeContigs.pl
# Author:         Glenn Tesler
# Created:        2010
# Last modified:  01/24/2010
#
# Copyright (c) 2009-2010 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# TODO:
#
#  Add command-line options
#
#  Summaraize many things about assembly, not just contig lengths
#
#  Give statistics about reads
#    number of reads used
#    number of reads discarded
#    mean coverage of contigs by reads
#
#  Give further statistics on graph, including
#    number of components
#    number of pairs of components, and number of self-dual components
#  Would be easy to do in C++ instead of perl, with Mark's routines.
#
#  Determine useful routines in PrintGraphSummary.cpp, some of which
#  need to be fixed up, and use them
#
#  Stats restricted to high-coverage contigs (mean coverage >= 5)


sub get_contig_lengths {
		my $nc = 0;
		my $tot=0;
		my $len = 0;
		my @lens;
		while (<>) {
				if (/^>/) {
						if ($len > 0) {
								push @lens, $len;
						}
						$len = 0;
						next;
				}
				
				chomp;
				$len += length($_);
		}

		if ($len > 0) {
				push @lens, $len;
		}

		# sort lengths from largest to smallest
		@lens = sort { $b <=> $a } @lens;

		return \@lens;
}

sub calc_stats {
		my ($lens0,$min_size) = @_;

		my $lens = [ grep { $_ >= $min_size } @$lens0 ];

#		print join "\n", @$lens;
#		print "\n\n";

    # number of contigs
		my $nc = scalar(@$lens);

		my ($min,$max) = (0,0);
		if ($nc > 0) {
				$max = $lens->[0];
				$min = $lens->[$nc-1];
		}

    # sum of contig lengths and squared lengths
		my ($tot,$tot2) = (0,0);
		$tot += $_     for @$lens;
		$tot2 += $_*$_ for @$lens;

		# mean and standard deviation
		my ($mean, $sd) = (0,0);
		if ($nc > 0) {
				$mean = $tot / $nc;
				if ($nc > 1) {
						$sd = sqrt(($tot2 - $nc*$mean*$mean)/($nc-1));
				}
		}

		# median
		my $med = 0;
		if ($nc > 0) {
				if ($nc & 1) {
						# odd # entries, take middle one
						$med = $lens->[($nc-1)/2];
				} else {
						# even # entries, average the middle two
						$med = ($lens->[($nc/2)-1] + $lens->[$nc/2]) / 2;
				}
		}

		# N50, N75, N90
		# Nx is the size of the largest contig s.t. x% of all bases are in
		# contigs of size >= Nx
		my @n_pct = (50,75,90);  # in order
		my @n_cutoffs = map { int( ($_*$tot + 99) / 100); } @n_pct;
		my @n_values;

		if ($nc == 0) {
				@n_values = (0) x scalar(@n_pct);
		} else {
				my $cutoff = shift @n_cutoffs;
				my $cum_sum = 0;

			LEN_LOOP:
				foreach my $len (@$lens) {
						$cum_sum += $len;
						while ($cum_sum >= $cutoff) {
								push @n_values, $len;
								$cutoff = shift @n_cutoffs;
								last LEN_LOOP if (!defined $cutoff);
						}
				}
		}

		my $stats = {
				'min_size' => $min_size,  # minimum size considered

				'num_contigs' => $nc,
				'tot' => $tot,

				'min' => $min,
				'max' => $max,
				'mean' => $mean,
				'sd' => $sd,
				'median' => $med,
				'Nx_x' => \@n_pct,
				'Nx_len' => \@n_values,
		};

		return $stats;
}

sub print_stats {
		my ($stats) = @_;

		my $min_size = $stats->{'min_size'};
		my $desc;
		if ($min_size == 0) {
				$desc = "all contigs";
		} else {
				$desc = "contigs >= $min_size nt";
		}
		print <<"END_report";
Contig statistics [$desc]:
  Number of contigs:       $stats->{'num_contigs'}
  Number of bases:         $stats->{'tot'}
  Min length:              $stats->{'min'}
  Max length:              $stats->{'max'}
  Mean length:             $stats->{'mean'}
  Standard deviation:      $stats->{'sd'}
  Median:                  $stats->{'median'}
END_report

    my $num_Nx = scalar(@{$stats->{'Nx_x'}});
		for (my $i=0; $i < $num_Nx; $i++) {
				my $Nx = sprintf("%-5s","N" . $stats->{'Nx_x'}->[$i] . ":");
				my $len = $stats->{'Nx_len'}->[$i];
				print << "END_Nx";
  ${Nx}                    $len
END_Nx
		}
		print "\n";
}

sub main {
		if ($ARGV[0] =~ /^-/) {
				print_usage();
				exit(0);
		}

		my $lens = get_contig_lengths();

		if (scalar(@$lens) == 0) {
				print_usage();
				exit(0);
		}

		my $stats0 = calc_stats($lens,0);
		my $stats500 = calc_stats($lens,500);
		print_stats($stats0);
		print_stats($stats500);
}

sub print_usage {
		print <<END_usage
usage:     SummarizeContigs.pl < fasta_file
suggested: SummarizeContigs.pl < reads.fasta.contig > reads.fasta.summary
NOTE: This may change, to input reads, graph, and contig files
END_usage
}

main();
1;
