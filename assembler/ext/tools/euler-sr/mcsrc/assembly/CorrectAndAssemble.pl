#!/usr/bin/env perl
use POSIX;
use FwgLib;

if ($#ARGV < 1 ) {
  print "usage: CorrectAndAssemble readsFile workingDir\n";
  print "            [-matesFile matesFile] [-vertexSize v] [-exeDir e]\n";
	print "            [-illumina]\n";
	print "            [-454]\n";
	print "  Selecting Illumina or 454 will adjust some parameters. The default \n";
	print "  is to assemble 454.\n";
  exit(0);
}
$vertexSize = 20;
$machtype = $ENV{"MACHTYPE"};
$srcDir = $ENV{"EUSRC"};
$exeDir = "$srcDir/assembly/$machtype";

$readsFile = shift @ARGV;
$workingDir = shift @ARGV;
$matesFile = "";
$readType = "454";
while ($#ARGV >= 0) {
  $option = shift @ARGV;
  if ($option eq "-matesFile") {
    $matesFile = shift @ARGV;
  }
  if ($option eq "-vertexSize") {
    $vertexSize = shift @ARGV;
  }
	if ($option eq "-illumina") {
			$readType= "illumina";
	}
	if ($option eq "-454") {
			$readType = "454";
	}
}
FwgLib::RunCommand("mkdir $workingDir");
FwgLib::RunCommand("ln -s ../$readsFile $workingDir/$readsFile");
FwgLib::RunCommand("$srcDir/assembly/$machtype/integralCountSpectrum $workingDir/$readsFile $vertexSize $workingDir/$readsFile.ispect -binary -printCount");
FwgLib::RunCommand("$srcDir/assembly/$machtype/sortIntegralTupleList $workingDir/$readsFile.ispect -printCount");

$expErrCmd = "$srcDir/assembly/$machtype/estErrorDist $workingDir/$readsFile.ispect -binSpect";
print "$expErrCmd\n";
$expErr = `$expErrCmd`;
$expErr /= 2;
if ($expErr == 0) {
		$expErr = 2;
}
print "Estimated the lower bound on correct tuples to be $expErr\n";
$exp2 = POSIX::floor($expErr * 1.5);

if ($readType eq "454" ) {
		
		FwgLib::RunCommand("$srcDir/assembly/$machtype/fixErrorsISAP $workingDir/$readsFile $workingDir/$readsFile.ispect $vertexSize $workingDir/$readsFile.fixed -discardFile $workingDir/$readsFile.discards -maxGap 3 -maxScore 5 -maxTrim 4 -edgeLimit 3 -minMult $expErr -startScore 2 -stepScore 2 -maxScore 6");
}
else {
		# Require at least 2x redundancy, but grab reads with low multiplicity.
		$catchall = POSIX::floor($exp2 / 2);
		if ($carchall == 1) {
				$catchall = 2;
		}
		FwgLib::RunCommand("$srcDir/assembly/$machtype/fixErrorsIVoting $workingDir/$readsFile $workingDir/$readsFile.ispect $vertexSize $workingDir/$readsFile.fixed -discardsFile $workingDir/$readsFile.discards -minMult $expErr -minVotes 2 -maxTrim 5 -catchall $catchall");
}

$cwd = $ENV{"PWD"};

$origReadsFile = $readsFile;
$fixedFile = "$readsFile.fixed";
FwgLib::RunCommand("$srcDir/assembly/assemble.pl $workingDir/$fixedFile -vertexSize $vertexSize");

#$vertexSize = 20;
$minEdgeLength = $vertexSize * 2;
if ($readType eq "illumina") {
		$bulgeSize = $vertexSize * 3;
}
else {
		$bulgeSize = $vertexSize * 5;
}

#
# Simplify errors in the graph.
#

FwgLib::RunCommand("mkdir $workingDir/Simple");
FwgLib::RunCommand("cd $workingDir/Simple; ln -s ../$fixedFile .");
FwgLib::RunCommand("$srcDir/assembly/$machtype/simplifyGraph $workingDir/$fixedFile $workingDir/Simple/$fixedFile.simple -minEdgeLength $minEdgeLength -removeBulges $bulgeSize -removeLowCoverage 5 3 -vertexSize $vertexSize");

@sources = `$srcDir/assembly/$machtype/printGraphSummary $workingDir/Simple/$fixedFile.simple -sources`;
$nSources = $#sources;
if ($nSources < 1000) {
		FwgLib::RunCommand("$srcDir/assembly/$machtype/joinsas $workingDir/Simple/$fixedFile.simple 10 $workingDir/Simple/$fixedFile.simple.j -vertexSize $vertexSize");
		$joinedFile = "$fixedFile.simple.j";
}
else {
		print "The graph is too fragmented.  Not attempting to de-fragment it.\n";
		$joinedFile = "$fixedFile.simple";
}

#
# Transform paths in the graph.
#

FwgLib::RunCommand("mkdir $workingDir/Transformed");
if ($readType eq "454") {
		FwgLib::RunCommand("$srcDir/assembly/ReorderAndLink.pl $joinedFile $fixedFile -vertexSize $vertexSize -workingDir $workingDir/Simple");
		if ($readType ne "illumina") {

				FwgLib::RunCommand("cd $workingDir/Transformed ; ln -s ../Simple/$joinedFile.r* .");
				
				$eqtransCmd = "$srcDir/euler/euler_et -s $workingDir/Transformed/$joinedFile.r -x $vertexSize -o $workingDir/Transformed/$joinedFile.et.dot -S";

				if ($matesFile ne "") {
						# Run equivalent transformation that does not split the graph.
						#
						FwgLib::RunCommand($eqtransCmd);
						FwgLib::RunCommand("$srcDir/euler/euler_db -i $workingDir/Transformed/$joinedFile.r -x $vertexSize -o $workingDir/Transformed/$joinedFile.db.dot -X -M 0 -m $matesFile");
						system("mv $workingDir/Transformed/$joinedFile.r.db.contig $origReadsFile.contig");
				}
				else {
						# No mate-paired transformations will take place, split the graph
						# to resolve paths that pass through partially resolved repeats.
						FwgLib::RunCommand("$eqtransCmd -X");
						FwgLib::RunCommand("mv $workingDir/Transformed/$joinedFile.r.et.contig $origReadsFile.contig");
						system("rm $workingDir/Transformed/*.ace");
				}
		}
}
else {
		FwgLib::RunCommand("$srcDir/assembly/$machtype/transformGraph $workingDir/Simple/$joinedFile $vertexSize $workingDir/Transformed/$joinedFile.t -notStrict");
		if (-e "$workingDir/Transformed/$joinedFile.t.edge") {
				FwgLib::RunCommand("$srcDir/assembly/$machtype/printContigs $workingDir/Transformed/$joinedFile.t");
				FwgLib::RunCommand("mv $workingDir/Transformed/$joinedFile.t.contig $readsFile.contig");
		}
		else {
			FwgLib::RunCommand("$srcDir/assembly/$machtype/printContigs $workingDir/Simple/$joinedFile.simple");
			FwgLib::RunCommand("mv $workingDir/Simple/$joinedFile.simple.contig $readsFile.contig");
		}
}

# clean up some files
system("rm $workingDir/$fixedFile*.edge");
system("rm $workingDir/$fixedFile*.graph");
system("rm $workingDir/$fixedFile*.path");
system("rm $workingDir/$fixedFile*.bgraph");
system("rm $workingDir/$origReadsFile.ispect");
system("rm $workingDir/$fixedFile*.intv");

