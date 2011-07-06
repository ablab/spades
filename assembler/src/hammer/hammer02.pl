#!/usr/bin/perl
use Getopt::Long;
use strict;
use POSIX;
use File::Basename;

## The script that runs all stages of Hammer 0.2.
## 
## Author: snikolenko

my $builddir = "/home/snikolenko/algorithmic-biology/assembler/build/hammer";
my $datadir = "data/hammer";
my $logdir = "alllogs";
my $fastqfile = $datadir . "/smalltest.fastq";
my $prefix = "test.preproc";

my $qvoffset = 64;
my $tau = 2;
my $freqCutoff = -1;
my $mem = 2;
my $threads = 1;
my $iterations = 1;

GetOptions ('tau=i' => \$tau, 'threads=i' => \$threads, 'cutoff=f' => \$freqCutoff, 'qv_off=i' => \$qvoffset, 'mem=f' => \$mem, 'fastq=s' => \$fastqfile, 'prefix=s' => \$prefix, 'iterations=i' => \$iterations, 'logdir=s' => \$logdir, 'datadir=s' =>\$datadir) or usage();

my $meminfull = $mem . "000000000";
my $kmerfile = "$datadir/$prefix.kmers";
my $preprocfile = "$datadir/$prefix.fastq";

sub execute {
	my $command = shift;
	print $command . "\n";
	`$command`;
}


sub make_iteration {

	## REMOVE TEMP FILES
	#execute "rm $datadir/$prefix.*";
	#execute "rm $datadir/reads.uf*";

	## PREPROC
	execute "$builddir/preproc $tau $qvoffset $fastqfile $preprocfile $kmerfile $threads > $logdir/log.preproc.txt";
	execute "$builddir/preproc2 $tau $qvoffset $fastqfile $preprocfile $kmerfile $threads > $logdir/log.preproc2.txt";


	## SPLIT
	my $numbytes = `ls -l $kmerfile.? | awk '{print \$5}' | awk '{ s += \$1 } END{print s}'`;
	print "got numbytes=$numbytes\n";
	chomp $numbytes;
	$numbytes = int($numbytes / (4* $threads )) + 1000;
	execute "cat $kmerfile.? | split -a 2 -C $numbytes -d  - $kmerfile\.dup\. > $logdir/log.split.txt";

	## CLUSTER
	die if (!($tau >= 1));
	my $threads2 = `ls -1 $kmerfile\.dup* | wc -l`;
	chomp $threads2;
	execute "$builddir/cluster $tau $kmerfile $kmerfile\.dup $meminfull 0 $threads2 $datadir > $logdir/log.cluster.txt 2> $logdir/log.cluster.2.txt";
	#execute "rm -rf $datadir/test.preproc.kmers.dup.* $datadir/test.preproc.kmers.*";

	## CENTER
	execute "$builddir/center $datadir/reads.uf $freqCutoff $threads > $logdir/log.center.txt";
	my $command = "cat";
	for (my $i=0; $i<$threads; ++$i) { $command .= " reads.uf.corr.$i"; }
	execute $command . " > $datadir/reads.uf.corr";
	#execute "rm -rf reads.uf.corr.*";

	## RECONSTRUCT
	execute "$builddir/reconstruct $datadir/reads.uf.corr $fastqfile 1 $qvoffset > $logdir/log.reconstruct.txt";
}

sub add_iterno {
	my $iterno = shift;
	my $dirprefix = shift;
	foreach my $l (@_) {
		execute "mv " . $dirprefix . $l . " " . $dirprefix . $iterno . "." . $l;
	}
}

my $keepgoing = 1;
for (my $iter=0; $keepgoing > 0; ++$iter) {
	if ($iter >= $iterations - 1) { $keepgoing = 0; }
	print "   === ITERATION " . ($iter+1) . " ===\n";

	make_iteration();

	add_iterno $iter, "", "reads.processed";
	add_iterno $iter, "$datadir/", "reads.uf.corr", "test.preproc.kmers", "test.preproc.fastq", "reads.uf";
	add_iterno $iter, "$logdir/", "log.split.txt", "log.cluster.txt", "log.cluster.2.txt", "log.preproc.txt", "log.reconstruct.txt";

	my $diffresult = execute "diff $fastqfile $iter.reads.processed";
	if (length($diffresult) == 0) { 
		print "Nothing changed in the reads! Stopping iterations.\n";
		$keepgoing = 0;
	}

	$fastqfile="$iter.reads.processed";
	print "\n";
}

