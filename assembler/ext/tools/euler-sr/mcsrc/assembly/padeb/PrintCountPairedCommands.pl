#!/usr/bin/env perl
if ($#ARGV < 6) {
		print "usage: PrintCountPairedCommands.pl seq startTuple endTuple strideTuple startClone endClone strideClone\n";
		exit(0);
}

$seq = shift @ARGV;
$startTuple = shift @ARGV;
$endTuple   = shift @ARGV;
$strideTuple = shift @ARGV;

$startClone = shift @ARGV;
$endClone   = shift @ARGV;
$strideClone = shift @ARGV;

$commandNum = 0;
$EUSRC = $ENV{"EUSRC"};
$mach  = $ENV{"MACHTYPE"};
$curdir = $ENV{"PWD"};
for ($tuple = $startTuple; $tuple <= $endTuple; $tuple += $strideTuple) {
		for ($clone = $startClone; $clone <= $endClone; $clone += $strideClone) {
				print "cd $curdir; $EUSRC/assembly/padeb/$mach/countPaired $seq $tuple $clone -histogram $seq.$tuple.$clone.hist -contigLengths $seq.$tuple.$clone.contigs -printNumDuplicated $seq.$tuple.$clone.duplicated -printN50 $seq.$tuple.$clone.n50 -printNumContigs $seq.$tuple.$clone.nc -printUniqueDuplicated $seq.$tuple.$clone.ud \n";
#				if ($commandNum % 2 != 0 and $commandNum > 0) {
#						print "&\n";
#				}
#				else {
#						print "\n";
#				}
				$commandNum++;
		}
}
