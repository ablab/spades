#!/usr/bin/perl
use strict;
use warnings;
use Math::Random qw(:all);
use Getopt::Long;

my $genomefile = "~/hammerdata/MG1655-K12.fasta";
my $bowtielogfile = "ecoli_lane1_cov/0.log.tinyhead";
my $outfileprefix = "graph";
my @resolutions = qw(1 10 100 250 500 1000);
my $color = "lightblue";

my $result = GetOptions ("genome=s" => \$genomefile,
		         "log=s"  => \$bowtielogfile,
		         "output=s"  => \$outfileprefix,
			 "color=s"  => \$color,
);

sub revcomp {
	my $res = shift;
	$res = reverse $res;
	$res =~ tr/ACGTacgt/TGCAtgca/;
	return $res;
}

sub sum_arr {
	my ($arr, $start, $end) = @_;
	my $res = 0;
	for (my $i = $start; $i < $end; ++$i) {
		$res += $arr->[$i];
	}
	return $res;
}

open GNM, $genomefile or die $!;
my $g = "";
while (<GNM>) {
	if ($_ =~ /^>/) { next; }
	if ($_ =~ /^([ACGT]+)$/) {
		$g = $g . $1;
	}
}
close GNM;

my $len = length($g);
my $rg = revcomp( $g );
print length($g) . "\n";
print length($rg) . "\n";

my @cnt;
for (my $i = 0; $i < $len; ++$i) { $cnt[$i] = 0; }

my @ok;
my @err;
for (my $i = 0; $i < 101; ++$i) { $ok[$i] = 0; $err[$i] = 0; }

$| = 1;

my $tmpletter;

my $read_cnt = 0;
open LOG, $bowtielogfile or die $!;
while (<LOG>) {
	++$read_cnt;
	if ($read_cnt % 100000 == 0) { print "."; }
	if ($read_cnt % 2000000 == 0) { print "\n"; }
	if ( $_ =~ /(.*)\t([+-])\t(.*)\t([0-9]*)\t([ACGT]*)\t([^\t]*)\t([0-9]*)\t([^\n]*)$/ ) {
		my $n = length($5);
		my $gs = substr( $g, $4, $n );
		#print "\t\t$gs\n";

		if ($5 eq $gs) { # if all is well, all is very simple
			for (my $i = 0; $i < $n; ++$i) {
				if ($cnt[$4 + $i] < ($n - $i)) {
					$cnt[$4 + $i] = $n - $i;
				}
			}
		} else { # otherwise, it's not :)
			my $min = 0;
			for (my $i=0; $i < $n; ++$i) {
				if (not (substr($gs, $i, 1) eq substr($5, $i, 1)) ) {
					for (my $j=$min; $j<$i; ++$j) {
						if ($cnt[$4 + $j] < ($j - $min)) {
							$cnt[$4 + $j] = $j - $min;
						}
					}
					$min = $i+1;
				}
			}
		}

		## count errors
		if ($2 eq "+") {
			for (my $i=0; $i < $n; ++$i) {
				$tmpletter = substr($5, $i, 1);
				if ($tmpletter eq 'N') { next; }
				if (substr($gs, $i, 1) eq $tmpletter) {
					$ok[$i]++;
				} else {
					$err[$i]++;
				}
			}
		} else {
			for (my $i=0; $i < $n; ++$i) {
				$tmpletter = substr($5, $i, 1);
				if ($tmpletter eq 'N') { next; }
				if (substr($gs, $i, 1) eq $tmpletter) {
					$ok[$n-1-$i]++;
				} else {
					$err[$n-1-$i]++;
				}
			}
		}
	}
}
close LOG;
print "\n";

for (my $i = 0; $i < 101; ++$i) { print (($ok[$i]+$err[$i])." "); } print "\n";
for (my $i = 0; $i < 101; ++$i) { print $err[$i] . " "; } print "\n";

open OUT, ">$outfileprefix.errors.tex" or die $!;
print OUT "\\addplot [sharp plot,color=$color,mark=+,mark size=2pt,mark repeat=5,smooth,line width=1.5pt] coordinates {\n";
for (my $i=0; $i<101; ++$i) {
	if ( $err[$i] + $ok[$i] > 0 ) {
		print OUT "($i, " . $err[$i] / ($err[$i] + $ok[$i]) . ")\n";
	} else {
		print OUT "($i, 0)\n";
	}
}
print OUT "};\n";
close OUT;

foreach my $resolution (@resolutions) {
	open OUT, ">$outfileprefix.$resolution.tex" or die $!;
	print OUT "\\addplot[const plot,fill=$color,draw=$color] coordinates {\n(-1, 0)\n";
	my $i = 0; my $sum;
	for ( ; $i < ($len / $resolution) - 1; ++$i) {
		$sum = sum_arr( \@cnt, $i * $resolution, ($i + 1) * $resolution);
		print OUT "($i, " . ($sum / $resolution) . ")\n";
	}
	$sum = sum_arr( \@cnt, $i * $resolution, $len);
	print OUT "($i, " . ($sum / ($len - $i * $resolution) ) . ")\n";
	print OUT "(".($i+1).", 0)\n};\n";
	close OUT;
}

die "OK!";


