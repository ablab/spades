#!/usr/bin/perl
use Getopt::Long;
use strict;
use POSIX;
use File::Basename;

## The script that runs all stages of Hammer 0.2.
## 
## Author: snikolenko

my $builddir = "build/hammer";
my $datadir = "data/hammer";
my $fastqfile = $datadir . "/smalltest.fastq";
my $prefix = "test.preproc";

my $qvoffset = 64;
my $tau = 2;
my $freqCutoff = -1;
my $mem = 2;
my $threads = 1;

GetOptions ('tau=i' => \$tau, 'threads=i' => \$threads, 'cutoff=f' => \$freqCutoff, 'qv_off=i' => \$qvoffset, 'mem=f' => \$mem, 'fastq=s' => \$fastqfile, 'prefix=s' => \$prefix) or usage();

my $meminfull = $mem . "000000000";
my $kmerfile = "$datadir/$prefix.kmers";
my $preprocfile = "$datadir/$prefix.fastq";

print "Removing temporary files...\n";
my $command = "rm $datadir/$prefix.*";
print $command . "\n";
`$command`;
my $command = "rm $datadir/reads.uf*";
print $command . "\n";
`$command`;

## PREPROC

my $command = "./" . $builddir . "/preproc $tau $qvoffset $fastqfile $preprocfile $kmerfile > log.preproc.txt";
print $command . "\n";
`$command`;


## SPLIT

my $numbytes = `ls -l $kmerfile.? | awk '{print \$5}' | awk '{ s += \$1 } END{print s}'`;
chomp $numbytes;
$numbytes = int($numbytes / (4* $threads )) + 1000;
$command = "cat $kmerfile.? | split -a 2 -C $numbytes -d  - $kmerfile\.dup\. > log.split.txt";
print $command . "\n";
system($command);


## CLUSTER

die if (!($tau >= 1));
my $threads2 = `ls -1 $kmerfile\.dup* | wc -l`;
chomp $threads2;
$command = "./$builddir/cluster $tau $kmerfile $kmerfile\.dup $meminfull 0 $threads2 $datadir > log.cluster.txt";
print $command . "\n";
system($command);


## CENTER

$command = "./$builddir/center $datadir/reads.uf $freqCutoff 1 > log.center.txt";
print $command . "\n";
system($command);


