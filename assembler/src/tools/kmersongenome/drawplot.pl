#!/usr/bin/perl
use strict;
use warnings;
use Math::Random qw(:all);
use Getopt::Long;

my $pref1 = "output";
my $pref2 = "res2";
my $pref3 = "res3";
my $pref4 = "res4";
my $pref5 = "res5";

my $kmin = 15;
my $kmax = 29;
my $begin = 1883063;
my $end = 1883076;
my $repeat = 20;
my $numinp = 1;

my @marks = qw(x + asterisk star o oplus oplus* otimes otimes* square square* triangle diamond);
my @color = qw(black!75 red!50 lightblue green blue darkgreen red orange black!50 red!80!blue black);

#print "\\addplot [const plot,color=$color[$ind],mark=$marks[$ind],mark size=2pt,mark repeat=$repeat,line width=1pt] coordinates {\n";


my $result = GetOptions ("kmin=i" => \$kmin,
			 "kmax=i" => \$kmax,
			 "begin=i" => \$begin,
			 "end=i" => \$end,
			 "repeat=i" => \$repeat,
			 "numprefs" => \$numinp,
			 "pref=s" => \$pref1,
			 "pref2=s" => \$pref2,
			 "pref3=s" => \$pref3,
			 "pref4=s" => \$pref4,
			 "pref5=s" => \$pref5,
);

my @inputs;
push @inputs, $pref1, $pref2, $pref3, $pref4, $pref5;

$begin++;
$end++;

my $ind = 0;

for ( my $i = 0; $i < $numinp; ++$i ) {
	my $prefix = $inputs[$i];
	for ( my $k = $kmin; $k <= $kmax; $k = $k + 2 ) {
		if ( ($k == 23) or ($k == 27) ) { next; }
		my $lines = "$begin,$end" . "p";
		my $plot = `sed -n '$lines' $prefix.$k.tex`;
		print "\%$prefix k=$k:\n";
		print "\\addplot \[const plot,color=" . $color[$ind] . ",fill=" . $color[$ind] . ",mark=none,mark size=2pt,mark repeat=$repeat,line width=1.5pt] coordinates {\n";
		print "(" . ($begin-2) . ", 0)\n" . $plot . "(" . ($end) . ", 0)\n};\n";
		++$ind;
	}
}

