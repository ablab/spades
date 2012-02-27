#!/usr/bin/env perl
if ($#ARGV < 1) {
  print "usage: $0 in reads [-vertexSize v] [-workingDir dir]\n";
  exit(0);
}
$in = shift @ARGV;
$reads = shift @ARGV;
$curdir = "";
$vertexSizeOpt = "";
while ($#ARGV >= 0) {
  $opt = shift @ARGV;
  if ($opt eq "-dir") {
    $curdir = shift @ARGV;
  }
	if ($opt eq "-vertexSize") {
			$v = shift @ARGV;
			$vertexSizeOpt = "-vertexSize $v";
	}
	if ($opt eq "-workingDir") {
			$workingDir = shift @ARGV;
			chdir($workingDir);
	}
}

$cwd = `pwd`;
print "running ReorderAndLink from $cwd\n";

if (! -e $reads) {
  print "Couldn't find $reads\n";
  exit(0);
}

$nr = `grep -c ">" $reads`;
chomp($nr);
$machtype=$ENV{"MACHTYPE"};
$src     =$ENV{"EUSRC"};
$cmd = "$src/assembly/$machtype/reorderIntervals $in $in.r $nr $vertexSizeOpt";
print "running: $cmd\n";
system($cmd);
system("ln -s $reads $in.r");
system("ln -s $in.edge $in.r.edge");
system("ln -s $in.graph $in.r.graph");

