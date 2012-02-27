#!/usr/bin/env perl
use POSIX;

$readsFile = shift @ARGV;
$nJobs     = shift @ARGV;
$spectFile = shift @ARGV;

open(RF, $readsFile) or die "cannot read $readsFile\n";

$nreads = 0;
# count the # reads
while(<RF>) {
  $line = $_;
  if ($line =~ /^>/) {
    $nreads++;
  }
}
print "got $nreads reads\n";
$readsPerJob = POSIX::floor($nreads / ($nJobs - 1));
print "using $readsPerJob reads per job\n";
close(RF);

open(RF, $readsFile) or die "cannot read $readsFile\n";
$cwd = $ENV{"PWD"};
$job = 0;
open(OUT, ">$readsFile.$job.split") or die "cannot open $readsFile.$job\n";
open(CMD, ">runfix.$job.csh") or die "cannot open $runfix.$job.csh\n";
print CMD "cd $cwd ; ~/projects/mcsrc/assembly/i686/fixerrorssap $readsFile.$job.split $spectFile $readsFile.$job.split.fixed ";
print CMD " -tupleSize 16 -minMult 7 -maxScore 10 -discards $readsFile.$job.discards > $readsFile.$job.split.output \n";

for ($i =0; $i < $nJobs; $i++ ) {
  $readNum = 0;
  while(<RF>) {
    $line = $_;
    if ($line eq "") {
      last;
    }
    if ($line =~ /^>/) {
      if ($readNum >= $readsPerJob and $i < $nJobs -1) {
	$job++;
	close OUT;
	close CMD;
	print "opening $readsFile.$job.split and runfix.$job.csh\n";
	open(OUT, ">$readsFile.$job.split") or die "cannot open $readsFile.$job\n";
	open(CMD, ">runfix.$job.csh") or die "cannot open $runfix.$job.csh\n";
	print CMD "cd $cwd ; ~/projects/mcsrc/assembly/i686/fixerrorssap $readsFile.$job.split $spectFile $readsFile.$job.split.fixed ";
	print CMD " -tupleSize 16 -minMult 7 -maxScore 10 -discards $readsFile.$job.discards > $readsFile.$job.split.output \n";
	$readNum = 0;
	last;
      }
      ++$readNum;
    }
    print OUT $line;
  }
}
     

