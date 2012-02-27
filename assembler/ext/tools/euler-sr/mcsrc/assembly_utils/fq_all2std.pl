#!/usr/bin/env perl 

# Author: lh3

use strict;
use warnings;
use Getopt::Std;

my $usage = qq(
Usage:   fq_all2std.pl <command> <in.txt>

Command: scarf2std      Convert SCARF format to the standard/Sanger FASTQ
         fqint2std      Convert FASTQ-int format to the standard/Sanger FASTQ
         sol2std        Convert Solexa/Illumina FASTQ to the standard FASTQ
         fa2std         Convert FASTA to the standard FASTQ
         fq2fa          Convert various FASTQ-like format to FASTA
         instruction    Explanation to different format
         example        Show examples of various formats

Note:    Read/quality sequences MUST be presented in one line.
\n);

die($usage) if (@ARGV < 1);

# Solexa->Sanger quality conversion table
my @conv_table;
for (-64..64) {
  $conv_table[$_+64] = chr(int(33 + 10*log(1+10**($_/10.0))/log(10)+.499));
}

# parsing command line
my $cmd = shift;
my %cmd_hash = (scarf2std=>\&scarf2std, fqint2std=>\&fqint2std, sol2std=>\&sol2std, fa2std=>\&fa2std,
				fq2fa=>\&fq2fa, example=>\&example, instruction=>\&instruction);
if (defined($cmd_hash{$cmd})) {
  &{$cmd_hash{$cmd}};
} else {
  die("** Unrecognized command $cmd");
}

sub fa2std {
  my %opts = (q=>25);
  getopts('q:', \%opts);
  my $q = chr($opts{q} + 33);
  warn("-- The default quality is set to $opts{q}. Use '-q' at the command line to change the default.\n");
  while (<>) {
	if (/^>(\S+)/) {
	  print "\@$1\n";
	  $_ = <>;
	  print "$_+\n", $q x (length($_)-1), "\n";
	}
  }
}

sub fq2fa {
  while (<>) {
	if (/^@(\S+)/) {
	  print ">$1\n";
	  $_ = <>; print;
	  <>; <>;
	}
  }
}

sub scarf2std {
  while (<>) {
	my @t = split(':', $_);
	my $name = join('_', @t[0..4]);
	print "\@$name\n$t[5]\n+\n";
	my $qual = '';
	@t = split(/\s/, $t[6]);
	$qual .= $conv_table[$_+64] for (@t);
	print "$qual\n";
  }
}

sub fqint2std {
  while (<>) {
	if (/^@/) {
	  print;
	  $_ = <>; print; $_ = <>; $_ = <>;
	  my @t = split;
	  my $qual = '';
	  $qual .= $conv_table[$_+64] for (@t);
	  print "+\n$qual\n";
	}
  }
}

sub sol2std {
  my $max = 0;
  while (<>) {
	if (/^@/) {
	  print;
	  $_ = <>; print; $_ = <>; $_ = <>;
	  my @t = split('', $_);
	  my $qual = '';
	  $qual .= $conv_table[ord($_)] for (@t);
	  print "+\n$qual\n";
	}
  }
}

sub instruction {

  print "
FASTQ format is first used in the Sanger Institute, and therefore
we take the Sanger specification as the standard FASTQ. Although
Solexa/Illumina reads file looks pretty much like the standard
FASTQ, they are different in that the qualities are scaled
differently. In the quality string, if you can see a character
with its ASCII code higher than 90, probably your file is in the
Solexa/Illumina format.

Sometimes we also use an integer, instead of a single character,
to explicitly show the qualities. In that case, negative
qualities indicates that Solexa/Illumina qualities are used.

";

}

sub example {
  my $exam_scarf = '
USI-EAS50_1:4:2:710:120:GTCAAAGTAATAATAGGAGATTTGAGCTATTT:23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 19 23 23 23 18 23 23 23
USI-EAS50_1:4:2:690:87:GTTTTTTTTTTTCTTTCCATTAATTTCCCTTT:23 23 23 23 23 23 23 23 23 23 23 23 12 23 23 23 23 23 16 23 23 9 18 23 23 23 12 23 18 23 23 23
USI-EAS50_1:4:2:709:32:GAGAAGTCAAACCTGTGTTAGAAATTTTATAC:23 23 23 23 23 23 23 23 20 23 23 23 23 23 23 23 23 23 23 23 23 12 23 18 23 23 23 23 23 23 23 23
USI-EAS50_1:4:2:886:890:GCTTATTTAAAAATTTACTTGGGGTTGTCTTT:23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23
USI-EAS50_1:4:2:682:91:GGGTTTCTAGACTAAAGGGATTTAACAAGTTT:23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 20 23 23 23 23 23 23 23 23 23 23 23 18 23 23 23 23
USI-EAS50_1:4:2:663:928:GAATTTGTTTGAAGAGTGTCATGGTCAGATCT:23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23 23
';

  my $exam_fqint = '
@4_1_912_360
AAGGGGCTAGAGAAACACGTAATGAAGGGAGGACTC
+4_1_912_360
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 21 40 40 40 40 40 40 40 40 40 26 40 40 14 39 40 40
@4_1_54_483
TAATAAATGTGCTTCCTTGATGCATGTGCTATGATT
+4_1_54_483
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 16 40 40 40 28 40 40 40 40 40 40 16 40 40 5 40 40
@4_1_537_334
ATTGATGATGCTGTGCACCTAGCAAGAAGTTGCATA
+4_1_537_334
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 21 29 40 40 33 40 40 33 40 40 33 31 40 40 40 40 18 26 40 -2
@4_1_920_361
AACGGCACAATCCAGGTTGATGCCTACGGCGGGTAC
+4_1_920_361
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 39 40 40 40 40 40 40 40 40 31 40 40 40 40 40 40 15 5 -1 3
@4_1_784_155
AATGCATGCTTCGAATGGCATTCTCTTCAATCACGA
+4_1_784_155
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 31 40 40 40 40 40
@4_1_595_150
AAAGACGTGGCCAGATGGGTGGCCAAGTGCCCGACT
+4_1_595_150
40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 30 40 40 40 40 40 40 40 40 40 20 40 40 40 40 40 14 40 40
';

  my $exam_sol = '
@SLXA-B3_649_FC8437_R1_1_1_610_79
GATGTGCAATACCTTTGTAGAGGAA
+SLXA-B3_649_FC8437_R1_1_1_610_79
YYYYYYYYYYYYYYYYYYWYWYYSU
@SLXA-B3_649_FC8437_R1_1_1_397_389
GGTTTGAGAAAGAGAAATGAGATAA
+SLXA-B3_649_FC8437_R1_1_1_397_389
YYYYYYYYYWYYYYWWYYYWYWYWW
@SLXA-B3_649_FC8437_R1_1_1_850_123
GAGGGTGTTGATCATGATGATGGCG
+SLXA-B3_649_FC8437_R1_1_1_850_123
YYYYYYYYYYYYYWYYWYYSYYYSY
@SLXA-B3_649_FC8437_R1_1_1_362_549
GGAAACAAAGTTTTTCTCAACATAG
+SLXA-B3_649_FC8437_R1_1_1_362_549
YYYYYYYYYYYYYYYYYYWWWWYWY
@SLXA-B3_649_FC8437_R1_1_1_183_714
GTATTATTTAATGGCATACACTCAA
+SLXA-B3_649_FC8437_R1_1_1_183_714
YYYYYYYYYYWYYYYWYWWUWWWQQ
';

  print qq(
solexa
======
$exam_sol
scarf
=====
$exam_scarf
fqint
=====
$exam_fqint
);
}
