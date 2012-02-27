#!/usr/bin/env perl
use POSIX;

$readsFile = shift @ARGV;
$nJobs     = shift @ARGV;


open(RF, $readsFile) or die "cannot read $readsFile\n";
$cwd = $ENV{"PWD"};
for ($job =0; $job < $nJobs; $job++ ) {
  open(CMD, ">compfixed.$job.csh") or die "cannot open $runfix.$job.csh\n";
  print CMD "cd $cwd ;  blastall -pblastn -d ../db/spneumoniae.fasta -K 1 -e 1e-20 -a 2 -m 8 -i $readsFile.$job.split.fixed > $readsFile.$job.split.fixed.blasttab;  ~/projects/mcsrc/assembly/VerifySeq.pl ../spneumoniae.fasta $readsFile.$job.split.fixed.blasttab $readsFile.$job.split.fixed  > $readsFile.$job.split.fixed.verified\n";
}
     

