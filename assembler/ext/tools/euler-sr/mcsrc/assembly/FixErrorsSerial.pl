#!/usr/bin/env perl
use POSIX;

print "This script is deprecated.  Please use FixErrorsSerialParam.pl instead\n";
print "It is still included with the distribution if you want to edit it \n";
print "and run the script yourself.\n";
exit(0);


if ($#ARGV != 2) {
  print "usage: FixErrorsSerial.pl reads coverage working_directory \n";
	print "  -grid   Use a grid environment to fix errors.\n";
  exit(1);
}
$readsFile = shift @ARGV;
$coverage  = shift @ARGV;
$workingDir = shift @ARGV;

$useGrid = 0;
while ($#ARGV >= 0) {
	$opt = shift @ARGV;
	if ($opt eq "-grid") {
		$useGrid = 1;
	}
}

$mach    = FwgLib::CrucialGetEnv("MACHTYPE");
$origDir = FwgLib::CrucialGetEnv("PWD");
$srcDir  = FwgLib::CrucialGetEnv("EUSRC");


chdir($workingDir);

$maxMult = POSIX::floor($coverage/2);
$minMult = POSIX::floor($coverage/20);
if ($minMult < 3) {
  $minMult = 3;
}
$steps  = 4;
$step  = POSIX::floor(($maxMult - $minMult)/$steps);
$mult = $maxMult;
$iter = 0;
`ln -s $readsFile $readsFile.$iter`;
while ($mult > $minMult) {
  $curReadsFile = "$readsFile.$iter";
  $cmd = "$srcDir/assembly/$machtype/countSpectrum $curReadsFile $curReadsFile.spect -tupleSize 20";
  system($cmd);
  $cmd = "$srcDir/assembly/$machtype/sortTupleList $curReadsFile.spect $curReadsFile.spect.sorted -minMult 2";
  system($cmd);
  $cmd = "$srcDir/assembly/$machtype/fixErrorsSAP $curReadsFile $curReadsFile.spect.sorted $curReadsFile.fixed -discardFile $curReadsFile.discards -minMult $mult -span 2 -edgeLimit 8 -maxTrim 4 -maxScore 7 -maxGap 3";
  print "running $cmd\n";
  system($cmd);
  $iter++;
  `cat $curReadsFile.fixed $curReadsFile.discards > $readsFile.$iter`;
  if ($mult == $minMult + 1) {
    #done
    $mult--;
  }
  else {
    if ($minMult + $step > $mult) {
      $mult = $minMult + 1;
    }
    else {
      $mult = $mult - $step;
    }
  }
}

`mv $curReadsFile.fixed $readsFile.fixed`;

  
    



