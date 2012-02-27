#!/usr/bin/env perl
use FwgLib;
if ($#ARGV < 1) {
		print "usage: CountSpectrum.pl readsFile vertexSize spectrumFile\n";
}
$readsFile = shift @ARGV;
$vertexSize = shift @ARGV;
$spectrumFile = shift @ARGV;

$mcsrc = FwgLib::CrucialGetEnv("EUSRC");
$machtype = FwgLib::CrucialGetEnv("MACHTYPE");
$exeDir = "$mcsrc/assembly/$machtype";
if ($vertexSize > 30) {
		FwgLib::RunCommand("$exeDir/countSpectrum $readsFile $spectrumFile -tupleSize $vertexSize");
		FwgLib::RunCommand("$exeDir/sortVertexList $spectrumFile $readsFile $vertexSize $spectrumFile.sorted");
		FwgLib::RunCommand("mv $spectrumFile.sorted $spectrumFile");
}
else {
		# prints the list of tuples that are present.
		$cmd1 = "$exeDir/integralCountSpectrum $readsFile $vertexSize $spectrumFile -binary";
		$cmd2 = "$exeDir/sortIntegralTupleList $spectrumFile";
		$cmd3 = "$exeDir/countIntegralTuples $spectrumFile $vertexSize $readsFile";
		$cmd4 = "$exeDir/binSpectToAscii $spectrumFile $vertexSize $spectrumFile.ascii -printCount";
		$cmd5 = "mv $spectrumFile.ascii $spectrumFile";
		$cmd6 = "$exeDir/sortTupleList $spectrumFile $spectrumFile.sorted";
		$cmd7 = "mv $spectrumFile.sorted $spectrumFile";
		$cwd = FwgLib::CrucialGetEnv("PWD");
		$logFile = "$cwd/countspectrum.log";
		FwgLib::RunLoggedCommand($cmd1, $logFile);
		FwgLib::RunLoggedCommand($cmd2, $logFile);
		FwgLib::RunLoggedCommand($cmd3, $logFile);
		FwgLib::RunLoggedCommand($cmd4, $logFile);
		FwgLib::RunLoggedCommand($cmd5, $logFile);
		FwgLib::RunLoggedCommand($cmd6, $logFile);
		FwgLib::RunLoggedCommand($cmd7, $logFile);
}

