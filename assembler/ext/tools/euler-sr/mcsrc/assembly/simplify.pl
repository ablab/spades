#!/usr/bin/env perl
use FwgLib;
if ($#ARGV < 0) {
	print "usage: simplify.pl infile \n";
	print "  There must exist infile.bgraph infile.edge infile.intv\n";
	print "  The advanced usage is: infile vertexSize minEdgeLength lowCoverage\n";
	exit(0);
}

$infile = shift @ARGV;
$vertexSize = 20;
$minEdgeLength = 100;
$bulgeSize = 100;
$lowCoverage = 3;
if ($#ARGV >= 0) {
  $vertexSize = shift @ARGV;
  $minEdgeLength = shift @ARGV;
  $bulgeSize = shift @ARGV;
  $lowCoverage = shift @ARGV;
}

$infile =~ /(.*)\.(.*)/;
$mach    = FwgLib::CrucialGetEnv("MACHTYPE");
$mcsrc   = FwgLib::CrucialGetEnv("EUSRC");
$outfile = "$infile.simple";
$cmd ="$mcsrc/assembly/$mach/simplifyGraph $infile $outfile -vertexSize $vertexSize -minEdgeLength $minEdgeLength -removeLowCoverage $lowCoverage 2 -removeBulges $bulgeSize -removeSimpleBulges $bulgeSize";
print "$cmd\n";
system($cmd);

