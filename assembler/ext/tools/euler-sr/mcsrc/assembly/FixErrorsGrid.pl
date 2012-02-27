#!/usr/bin/env perl

use FwgLib;

if ($#ARGV < 1) {
		PrintUsage();
		exit(0);
}

$commandStr = join(" ", @ARGV);

$readsFile    = shift @ARGV;
$workingDir   = shift @ARGV;

# The results are should be stored in $readsFile.fixed, and $readsFile.discards

#
# Configure the working directory.  This is necessary because
# each job that is spawned needs to chdir to $workingDir because
# they are re-started in the same directory as where the parent
# script started.
#
$readsPerFile = 30000;
$discardFile  = "";
$smp           = 0;

$nJobs = -1;
$nProc = 2;
$vertexSize = 20;
$curDir = "";
$minMult = 5;

$minVotes = -1;
$minVotesOpt = "";
$maxGapOpt = "";
$startScoreOpt = "";
$stepScoreOpt  = "";
$maxScoreOpt   = "";
$maxTrimOpt    = "";
$spanOpt       = "";
$edgeLimitOpt  = "";
$searchOpt     = "";
$logFileOpt    = "";
$logFile       = "";
$printMap      = "";
$countSpectrum = 1;
$printMap      = 0;
$csTrimFrontOpt = "";
$csTrimEndOpt = "";
while ($#ARGV >= 0) {
		$opt = shift @ARGV;
		if($opt eq "-curDir" ) {
				$curDir = shift @ARGV;
		}
		elsif($opt eq "-readsPerFile") {
				$readsPerFile = shift @ARGV;
		}
		elsif($opt eq "-vertexSize") {
				$vertexSize = shift @ARGV;
		}
		elsif($opt eq "-trimFront") {
			$trimFront = shift @ARGV;
			$csTrimFrontOpt = "-trimFront $trimFront";
		}
		elsif($opt eq "-trimEnd") {
			$trimEnd = shift @ARGV;
			$csTrimEndOpt = "-trimEnd $trimEnd";
		}
		elsif($opt eq "-discardFile") { 
				$discardFile = shift @ARGV;
		}
		elsif($opt eq "-prog" ) {
				$fixProg = shift @ARGV;
		}
		elsif($opt eq "-search" ){
				$search = shift @ARGV;
				$searchOpt = "-search $search";
		}
		elsif($opt eq "-nJobs" ){
				$nJobs = shift @ARGV;
		}
		elsif($opt eq "-edgeLimit") {
				$edgeLimit = shift @ARGV;
				$edgeLimitOpt = "-edgeLimit $edgeLimit";
		}
		elsif($opt eq "-trim") {
				$trim  = 1;
			  $trimOpt = "-trim";
		}
		elsif($opt eq "-maxTrim" ){
				$maxTrim    = shift @ARGV;
				$maxTrimOpt = "-maxTrim $maxTrim";
		}		
		elsif($opt eq "-trim") {
			$trimOpt = "-trim";
		}
		elsif($opt eq "-maxGap") {
				$maxGap = shift @ARGV;
				$maxGapOpt = "-maxGap $maxGap";
		}
		elsif($opt eq "-stepScore") {
				$stepScore = shift @ARGV;
				$stepScoreOpt = "-stepScore $stepScore";
		}
		elsif($opt eq "-maxScore") {
				$maxScore = shift @ARGV;
				$maxScoreOpt = "-maxScore $maxScore";
		}
		elsif($opt eq "-startScore") {
				$startScore = shift @ARGV;
				$startScoreOpt = "-startScore $startScore";
		}
		elsif($opt eq "-minVotes") {
				$minVotes = shift @ARGV;
				$minVotesOpt = "-minVotes $minVotes";
		}
		elsif($opt eq "-nProc" ){ 
				$nProc = shift @ARGV;
				$nProcOpt = "-nProc $nProc";
		}
		elsif($opt eq "-minMult") {
				$minMult = shift @ARGV;
		}
		elsif($opt eq "-smp") {
				$smpOpt = "-smp";
				$smp = 1;
		}
		elsif($opt eq "-span") {
				$span = shift @ARGV;
				$spanOpt = "-span $span";
		}
		elsif($opt eq "-log" ){ 
				$logFile = shift @ARGV;
				open(LOG, ">>$logFile") or die "cannot open $logFile\n";
				$logFileOpt = "-log $logFile";
		}
		elsif($opt eq "-spectrum") {
				$spectrumFile = shift @ARGV;
				$countSpectrum = 0;
		}
		elsif($opt eq "-map" ) {
			$printMap = 1;
			$mapFile = shift @ARGV;
		}
		else {
				PrintUsage();
				print "bad option: $opt\n";
				exit(0);
		}
}

FwgLib::PrintLog($logFile, $commandStr);

# Move to the correct starting directory if necessary
$curDirOpt = "";
if ($curDir ne "" ) {
		chdir($curDir);
		$curDirOpt = "-curDir $curDir";
}
else {
	$curDir = FwgLib::CrucialGetEnv("PWD");
}

if ($workingDir eq "") {
		print "ERROR, must specify a working directory\n";
		exit(0);
}

# move to current directory so this may be called from
# a program that is moving around.

print "fixerrorsgrid setting up working dir $workingDir\n";
$localReadsFile = FwgLib::SetupWorkingDir($workingDir, "$readsFile");

# Configure environment variables
$srcDir = FwgLib::CrucialGetEnv("EUSRC");
$mach   = FwgLib::CrucialGetEnv("MACHTYPE");
print "src: $srcDir mach: $mach\n";

# Configure fast access to executables
$asm      = "$srcDir/assembly/$mach";
$asmScript= "$srcDir/assembly/";
$utils    = "$srcDir/UTILS";

$splitExt = $$;
$splitBase= "split";

# Divide the input file into more manageable files

$splitCmd = "$srcDir/UTILS/SplitReads.pl $localReadsFile $readsPerFile $splitBase $splitExt";
print "running: $splitCmd\n";
$output   = `$splitCmd`;
print "res: $output";

# Store the names of split files and number of split files
@filesToFix = glob("$splitBase.*.$splitExt");
$nFilesToFix = scalar @filesToFix;
if ($nJobs == -1) {
		$nJobs = $nFilesToFix;
}


if ($countSpectrum) {
# Compute the spectrum of all these files, then combine them together into one.
# use the 'alreadySplit' option for computing the spectrum so that the 
# reads aren't divided again.
		$computeSpectCommand = "$srcDir/assembly/CountSpectrumGrid.pl $localReadsFile " .
		" $workingDir $curDirOpt " . 
		" -minMult $minMult " .
		" -outputFile $localReadsFile.spect -vertexSize $vertexSize" . 
		" $smpOpt $nProcOpt $logFileOpt" .
		" $csTrimFrontOpt $csTrimEndOpt" . 
		" -alreadySplit $splitBase $splitExt";
# actually run it
		print "going to compute spectrum with: $computeSpectCommand\n";
		FwgLib::PrintLog($logFile, "compuing spectrum of $localReadsFile\n");
		FwgLib::PrintLog($logFile, $computeSpectCommand);
		system("$computeSpectCommand");
		if ($? != 0) {
				print "error running $computeSpectCommand\n";
				exit(0);
		}
		$spectrumFile = "$localReadsFile.spect";
}

# 3. Fix errors 

# 3.1 create the error fixing commands

# configure the type of error fixing that will happen

if ($fixProg eq "vote") {
		$fixExe = "$asm/fixErrorsVoting ";
}
elsif ($fixProg eq "sap") {
		$fixExe = "$asm/fixErrorsSAP ";
}
else {
		print "Cannot find program $fixProg\n";
		exit(1);
}



@fixCommands = ();
foreach $file (@filesToFix) {
		if ($discardFile eq "") {
				$discardOpt = "";
		}
		else {
				$discardOpt = " -discardFile $file.discards ";
		}
		$fixedFile = "$file.fixed";
		if ($printMap) {
			$mapFile    = "$file.map";
			$mapFileOpt = "-map $mapFile";
		}
		push @fixCommands, "cd $workingDir; " .
				" $fixExe $file $spectrumFile $vertexSize $file.fixed ".
		    " $minMultOpt $discardOpt " .
				" $startScoreOpt $stepScoreOpt $maxScoreOpt $spanOpt " .
				" $edgeLimitOpt $maxTrimOpt  $maxGapOpt" . 
				" $minVotesOpt $searchOpt $trimOpt $mapFileOpt";
}


# Print some information about the commands to run
$numFixCommands = scalar @fixCommands;
#print "Fixing with commands $numFixCommands: \n";
#foreach $fc (@fixCommands) {
#		print "  $fc\n";
#}

FwgLib::SubmitCommandListAndWait($workingDir, $nJobs, $nProc, $smp, @fixCommands);

print "Split error correction done, collecting results\n";
# now collect the results
@fixed = ();
@discards = ();
foreach $file (@filesToFix) {
		push @fixed, "$file.fixed";
		push @discards, "$file.discards";
		push @map, "$file.map";
}
print "cur dir is: $curDir\n";
$catFixedCmd = "cd $workingDir; cat @fixed > $curDir/$readsFile.fixed";
if ($discardFile ne "") {
		$catAllCmd   = "cd $workingDir; cat @discards > $curDir/$readsFile.discards";
}
if ($printMap) {
	$catMapCmd = "cd $workingDir; cat @map > $curDir/$readsFile.map";
	`$catMapCmd`;
}

$res = `$catFixedCmd`;
$res = `$catAllCmd`;
`rm @fixed`;
`rm @discards`;
`rm @filesToFix`;
if ($printMap) {
#	`rm @map`;
}

#`rm @filesToFix`;
print "Done.\n";

sub PrintUsage() {
		print "usage: FixErrorsGrid.pl reads.fasta workingDir\n";
		print "options: \n";
		print "  -readsPerFile R   Make each sub-problem contain at most R reads.\n";
		print "  -vertexSize   V   Fix on vertex size V\n";
		print "  -discardFile file Print discards to 'discardFile'.  If this is not\n";
		print "                      used, all reads will be printed to the original file.\n";
		print "  -curDir dir       Start running in 'curDir'.  This is necessary when the script\n";
		print "                      is called from inside another script.\n";
		print "  -prog [vote|sap] options \n";
		print "                    Use either voting, or dp, and supply with options.  This\n";
		print "                      should be in quotes.  Include the discard file in the option\n";
		print "                      above, since it is processed differently.\n";
		print "  -nJobs N          Submit N jobs to the grid.\n";
		print "  -nProc P          Each job submits to P processors.\n";
		print "  -smp              Run in symmetric multiprocessor mode.  This doesn't submit\n";
		print "                      to a grid, but runs a few jobs at the same time to the \n";
		print "                      local machine.  Use this for fixing erros on a multi-core \n";
		print "                      machine, like a mac pro.\n";
		print "  -map mapFile      Print a map of the aligned reads to 'mapFile'\n";
		print "  -trim             Trim parts of reads that are not solid\n";
		print "  -maxTrim M        If a read must be trimmed by more than 'M' bases on each end \n";
		print "                      to be solid, it is discarded.\n";
		print "  -minVotes V       When performing error correction by voting, require at least\n";
		print "                      'V' votes to make a rea solid.\n";
		print "  -minMult  M       Require tuples of multiplicity at least 'M'\n";
}
