#!/usr/bin/env perl

if ($#ARGV < 4) {
		print "usage RunEulerDB.pl readfile mates file #mates rulefile niter [-dir workingdir] [-w low cov -v vertexSize]\n";
		exit(1);
}
$readsFile = shift @ARGV;
$matesFile = shift @ARGV;
$nmates    = shift @ARGV;
$ruleFile  = shift @ARGV;
$niter = shift @ARGV;

$workingDir = "";
$lowCoverage = 3;
$vertexSize = 20;
while($#ARGV >= 0) {
		$opt = shift @ARGV;
		if ($opt eq "-dir") {
				$workingDir = shift @ARGV;
		}
		if ($opt eq "-w" ) {
				$lowCoverage = shift @ARGV;
		}
		if ($opt eq "-v" ) {
				$vertexSize = shift @ARGV;
		}
}


if ($workingDir ne "" ) {
		chdir ($workingDir);
}

$prev = $readsFile;

for ($m = 0; $m < $nmates; $m++) {
		for ($i = 0; $i < $niter; $i++ ){

				if ($i == 0 && $m == 0) {
						$cmd = "~/projects/mcsrc/euler/euler_db -i $readsFile -m $matesFile -r $ruleFile -x $vertexSize -o $readsFile.db.dot -w 2 -M $m";
						print "$cmd\n";
						$output = `$cmd`;
						if ($? != 0) {
								print "FAILED0: $cmd\n";
								exit(0);
						}
						$prev = $readsFile;
				}
				else {
						`ln -s $readsFile $readsFile.$m.$i`;
						
						`ln -s $prev.db.edge $readsFile.$m.$i.edge`;
						`ln -s $prev.db.graph $readsFile.$m.$i.graph`;
						`ln -s $prev.db.intv $readsFile.$m.$i.intv`;
						$cmd = "~/projects/mcsrc/euler/euler_db -i $readsFile.$m.$i -x $vertexSize -o $readsFile.$m.$i.dot -e $readsFile.$m.$i.edge -g $readsFile.$m.$i.graph -p $readsFile.$m.$i.intv -M $m -m $matesFile -r $ruleFile";
						if ($i == $niter-1) {
								$cmd .= " -X";
						}
						print "$cmd\n";
						$output = `$cmd`;
						if ($? != 0) {
								print "ERROR with command: $cmd\n";
								print "output: $output\n";
						}
						$prev = "$readsFile.$m.$i";
				}
		}
}


