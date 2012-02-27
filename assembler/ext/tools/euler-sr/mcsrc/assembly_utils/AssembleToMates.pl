#!/usr/bin/env perl

$reads = shift @ARGV;
$vertexSize = shift @ARGV;
$rule = shift @ARGV;
if ($#ARGV >= 0) {
		$workingDir = shift @ARGV;
		chdir($workingDir);
		print "working in $workingDir\n";
		$simple = $workingDir . "/Simple";
		$mt     = $workingDir . "/MateTransformed";
		$wdreads = $workingDir . "/$reads";
		$wd = "$workingDir/";
} else {
		$simple = "Simple";
		$mt = "MateTransformed";
		$wdreads = $reads;
		$wd = "";
}

`mkdir -p $simple`;
`mkdir -p $mt`;
$mach = $ENV{"MACHTYPE"};
`~/projects/mcsrc/assembly/assemble.pl $wdreads -vertexSize $vertexSize`;
`~/projects/mcsrc/assembly/$mach/buildMateTable $wdreads ~/projects/mcsrc/assembly/readtitle.rules $wdreads.matetable`;
system("~/projects/mcsrc/assembly/$mach/simplifyGraph $wdreads $simple/$reads -removeSimpleBulges 75 -removeBulges 75");
system("~/projects/mcsrc/assembly/$mach/mateTransformGraph $simple/$reads $wdreads.matetable ~/projects/mcsrc/assembly/readtitle.rules $mt/$reads -findPaths -ruleType $rule -notStrict > $wd$reads.out");
