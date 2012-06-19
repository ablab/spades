#!/usr/local/bin/perl -w
# Konrad Paszkiewicz, University of Exeter UK.
# k.h.paszkiewicz@exeter.ac.uk

use strict;

my $usage = "This script is designed to be used after filtering paired-end
Illumina reads. Such filtering can end up removing one member of a pair. To 
ensure compatibility with the shuffleSequences_fasta.pl script, the reads in
the forward and reverse file must be in the same order. 
 
This script will take the forward reads and reverse read files in FASTA format
and ensures the reads are presented in the correct order. Reads which appear in
the forward or reverse direction ONLY will be put into the singletons file. The
output files can then be run through the shuffleSequences_fasta.pl script. The
singletons can be added as an additional channel to aid in the velvet assembly
using the -short2 parameter in velveth.
 
USAGE: 

select_paired.pl <file containing forward reads in FASTA format> <name of output file for sorted forward reads> <name of file containing reverse reads in FASTA format> <name of output file for sorted reverse reads> <singletons output file (i.e. reads present in only one of the input files>
";

my $infile1 = shift or die $usage;
my $outfile1 = shift or die $usage;

our $infile2 = shift or die $usage;
my $outfile2 = shift or die $usage;

our $outfile3 = shift or die $usage;

my %hash1;
my %hash2;

open(FILE1, $infile1) or die "Cannot open $infile1\n";
open(FILE2, $infile2) or die "Cannot open $infile2\n";

open(OUTFILE1, ">$outfile1") or die "Cannot open $outfile1\n";
open(OUTFILE2, ">$outfile2") or die "Cannot open $outfile2\n";
open(OUTFILE3, ">$outfile3") or die "Cannot open $outfile3\n";

my $name1;
my $name2;

while(<FILE1>){
        if(/^(\>.*)\/\d$/){
                $hash1{$1}=1;
                $name1=$1;
        }else{
                $hash1{$name1}=$_;
        }
}
close(FILE1);

while(<FILE2>){
        if(/^(\>.*)\/\d$/){
                $hash2{$1}=1;
                $name2=$1;
        }else{
                $hash2{$name2} = $_;
        }
}
close(FILE2);

for $name1 ( keys %hash1 ) {
        if(exists $hash2{$name1}){
                print OUTFILE1 "$name1/1\n$hash1{$name1}";
                print OUTFILE2 "$name1/2\n$hash2{$name1}";
        }else{
                print OUTFILE3 "$name1/1\n$hash1{$name1}";
        }
}

for $name2 (keys %hash2){
        if(!exists $hash1{$name2}){
                print OUTFILE3 "$name2/2\n$hash2{$name2}";
        }
}

close(OUTFILE1);
close(OUTFILE2);
