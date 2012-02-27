#!/usr/bin/env perl
use POSIX;

$readsFile = shift @ARGV;
$nJobs     = shift @ARGV;


open(RF, $readsFile) or die "cannot read $readsFile\n";
$cwd = $ENV{"PWD"};
for ($job =0; $job < $nJobs; $job++ ) {
  open(CMD, ">verify.$job.csh") or die "cannot open $runfix.$job.csh\n";
	print CMD "cd $cwd ; ~/projects/mcsrc/assembly/VerifySeq.pl ../spneumoniae.fasta $readsFile.$job.blasttab $readsFile.$job.split.fixed > $readsFile.$job.split.verified \n";
}
     

