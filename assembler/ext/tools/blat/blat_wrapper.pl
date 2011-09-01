#!/usr/bin/perl
use warnings;
use strict;
use List::Util qw(max);

# alignes fasta file to reference & outputs alignment stats to STDOUT
# expectes fasta file (to align) and reference sequence (fasta)

die "\n\tusage: blat_wrapper.pl <to-align.fa> <ref.fa>\n\n" if (@ARGV != 2);
my ($to_align_f, $ref_seq_f) = @ARGV;
my $aln = "out1.psl";

######################### align ###########################

`./blat $ref_seq_f $to_align_f $aln 1>/dev/null`;

# read PSL output, parse, and print shorthand to screen
my $curr_score;
my $name = "";
my @scores;
open(PSL, "<", $aln) or die "could nop open PSL output. Quitting...\n";
while(my $l = <PSL>){
	if ($l =~/^psL/ or $l eq "\n"){ next;} # don't write  
	if($l =~ /^---/){ print $l; next; } # write separator line
	my @l_arr = split(/\t/,$l);
	if($l !~ /^[0-9]+/){
		print_blat_line(\@l_arr);
		next;
	}
	else{
		$name = print_blat_line(\@l_arr);
		$curr_score = get_aln_score(\@l_arr);
		push @scores, $curr_score;
	}
}
close PSL;

if(@scores == 0){
	print "\nBest Score: no alignment\n\n";	
}
else{
	print "\nBest Score $name: " . max(@scores) ."\n\n";
}

system("rm $aln");

###########################################################
#################### Subroutines ##########################
###########################################################

sub print_blat_line{
	my ($line) = @_;	
	for my $i (0 .. (scalar @{$line} - 1)){
			print $line->[$i] if ($i != 13 and $i != 14);
			print "\t" if ($i < (scalar @{$line} - 1) and $i != 13 and $i != 14);
	}
	return $line->[9];	
}

sub get_aln_score{
	my ($line) = @_;
	my $score = 0;
	if(scalar @{$line} < 20){ return 0; }
	else{
		# get alignment score 
		$score = $line->[0] + ($line->[2]>>1) - $line->[1] - $line->[4] -$line->[6];
	}
	return $score;
}
