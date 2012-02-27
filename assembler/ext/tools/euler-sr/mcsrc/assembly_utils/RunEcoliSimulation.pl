#!/usr/bin/env perl
if ($#ARGV < 4) {
		print "usage: RunEcoliSimulation.pl workingDir genome readlength npasses resultsDir\n";
		print "  This simluations npasses through the genome, readlength reads from a clone\n";
		print "  of a fixed size (300bases now). \n";
		exit(0);
}
$workingDir = shift @ARGV;
$genome  = shift @ARGV;
$readLength = shift @ARGV;
$nPasses = shift @ARGV;
$resultsDir = shift @ARGV;

$workingDir .= "/mchaisso.$$";
system("mkdir -p $workingDir");
$mach = $ENV{"MACHTYPE"};
chdir($workingDir);
$cloneLength = 300 - ($readLength * 2);
system("~/projects/mcsrc/assembly/$mach/readsimulator $genome $workingDir/reads.fa -rl $readLength -cloneLib  $cloneLength 30 100 -stratify 1 $nPasses");
system("~/projects/mcsrc/assembly_utils/AssembleToMates.pl reads.fa 24 5 $workingDir");
system("echo $readLength $nPasses > $resultsDir/ecolisim.$readLength.$nPasses.$$");
system("pcl $workingDir/MateTransformed/*.edge | N50 >> $resultsDir/ecolisim.$readLength.$nPasses.$$");
system("rm -rf $workingDir");


