#!/usr/bin/env perl

if ($#ARGV < 1) {
		print "usage RunEulerET.pl readfile niter [-dir workingdir] [-w low cov]\n";
		exit(1);
}
$readsFile = shift @ARGV;
$niter = shift @ARGV;

$workingDir = "";
$lowCoverage = 3;
while($#ARGV >= 0) {
		$opt = shift @ARGV;
		if ($opt eq "-dir") {
				$workingDir = shift @ARGV;
		}
		if ($opt eq "-w" ) {
				$lowCoverage = shift @ARGV;
		}
}

if ($workingDir ne "" ) {
		chdir ($workingDir);
}


$cmd = "~/projects/mcsrc/euler/euler_et -s $readsFile -x 20 -o $readsFile.et.dot -w 3";
		$output = `$cmd`;
		if ($? != 0) {
				print "ERROR with command: $cmd\n";
				print "output: $output\n";
		}

$prev = $readsFile;
for ($i = 1; $i < $niter; $i++ ){
		`ln -s $readsFile $readsFile.$i`;
		
		`ln -s $prev.et.edge $readsFile.$i.edge`;
		`ln -s $prev.et.graph $readsFile.$i.graph`;
		`ln -s $prev.et.intv $readsFile.$i.intv`;
		$cmd = "~/projects/mcsrc/euler/euler_et -s $readsFile.$i -x 20 -o $readsFile.$i.dot";
		$output = `$cmd`;
		if ($? != 0) {
				print "ERROR with command: $cmd\n";
				print "output: $output\n";
		}
		$prev = "$readsFile.$i";
}




