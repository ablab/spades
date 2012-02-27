#!/usr/bin/env perl
use FwgLib;

if ($#ARGV < 1) {
  print "usage: CountSpectrumGrid.pl seqFile workingDir [-readsPerFile r] [-nJobs n] [-vertexSize v] [-alreadySplit splitBase splitExt] \n";
	print " -lockFile file Have worker jobs wait on 'file' to read data.\n";
  print "       -smp    Run in smp mode.\n";
  print "    -minMult M Only save tuples of multiplicity >= M\n";
  exit(1);
}

$commandStr = join(" ", @ARGV);
$inFile     = shift @ARGV;
$workingDir = shift @ARGV;

# split this amongst nJobs. Default is nreads / readsPerFile jobs
$nJobs        = -1; 
# Split file into fils of 50000 reads
$readsPerFile = 50000;
# Count on vertex size
$vertexSize   = 20;
$outputFile = "$inFile.spect";
$lockFile = "";
$curDir = "";
$doSplit = 1;
$splitBase = "split";
$splitExt  = $$;
$smp = 0;
$nProc = 0;

$trimEndOpt = "";
$trimFrontOpt = "";
$lockFileOpt = "";
$minMult = -1;

while ($#ARGV >= 0) {
  $opt = shift @ARGV;
  if ($opt eq "-readsPerFile") {
    $readsPerFile = shift @ARGV;
  }
  elsif ($opt eq "-nJobs") {
    $nJobs = shift @ARGV;
  }
  elsif ($opt eq "-nProc") {
    $nProc = shift @ARGV;
  }
  elsif ($opt eq "-minMult") {
    $minMult = shift @ARGV;
  }
  elsif ($opt eq "-vertexSize") {
    $vertexSize = shift @ARGV;
  }
  elsif ($opt eq "-outputFile" ) {
    $outputFile = shift @ARGV;
  }
  elsif ($opt eq "-curDir") {
    $curDir = shift @ARGV;
		print "starting in $curDir\n";
  }
	elsif ($opt eq "-trimFront") {
		$trimFront = shift @ARGV;
		$trimFrontOpt = "-trimFront $trimFront";
	}
	elsif ($opt eq "-trimEnd") {
		$trimEnd = shift @ARGV;
		$trimEndOpt = "-trimEnd $trimEnd";
	}
	elsif ($opt eq "-log" ){
			$logFile =shift @ARGV;
			open(LOG, ">>$logFile") or die "cannot write to $logFile\n";
	}
  elsif ($opt eq "-alreadySplit") {
    $splitBase = shift @ARGV;
    $splitExt  = shift @ARGV;
    $doSplit   = 0;
  }
  elsif ($opt eq "-smp") {
			$smp = 1;
  }
	elsif ($opt eq "-lockFile") {
		$doLock = 1;
		$lockFile = shift @ARGV;
		$lockFileOpt = "-lockFile $lockFile";
	}
  else {
    print "ERROR: bad command line argument: $opt\n";
    exit(1);
  }
}

FwgLib::PrintLog($logFile, "$0 $commandStr");
		
if ($curDir ne "") {
  chdir($curDir);
}
else {
	$curDir = $env{"PWD"};
}

$localFile = FwgLib::RemovePath($inFile);
$cwd = $ENV{"PWD"};

FwgLib::SetupWorkingDir($workingDir, $inFile);
print "done setting up working dir\n";

$EUSRC = $ENV{"EUSRC"};
$arch  = $ENV{"MACHTYPE"};

# split the input file .
if ($doSplit != 0) {
  $splitCommand = "$EUSRC/UTILS/splitfasta.pl $localFile $readsPerFile $splitBase $splitExt";
  $result = `$splitCommand`;
  print "for: $splitCommand\n";
  print "result: $result\n";
}

# find out what was split
@splitFiles = glob("$splitBase.*.$splitExt");
if ($nJobs == -1) {
  $nJobs = scalar @splitFiles;
}
# create the commands to run the jobs on the grid.
@commands = ();
for ($c = 0; $c <= $#splitFiles; $c++) {
		push @commands, "cd $workingDir; $EUSRC/assembly/CountSpectrum.pl $splitFiles[$c] $vertexSize $splitFiles[$c].spect";

#$EUSRC/assembly/$arch/sortTupleList $splitFiles[$c].spect $splitFiles[$c].spect.sorted $trimFrontOpt $trimEndOpt $lockFileOpt";
}

# count tuple frequency of each sub file
print "going to sumit commands\n";
print "workingDir $workingDir\n";
print "njobs: $nJobs\n";
print "nProc: $nProc\n";
print "smp: $smp\n";
FwgLib::SubmitCommandListAndWait($workingDir, $nJobs, $nProc, $smp, @commands);

# collect frequencies into one file

if ($minMult != -1) {
  $minMultOpt = "-minMult $minMult";
}
else {
  $minMultOpt = "";
}

$cwd = `pwd`;
chomp $cwd;
@collate[0] = "cd $cwd; $EUSRC/assembly/$arch/combineTupleLists $workingDir/$outputFile $workingDir/*.sorted  $minMultOpt -bufferOutput";
print "Collating to $cwd  $workingDir/$outputFile. $collate[0]\n";
FwgLib::SubmitCommandListAndWait($workingDir, 1, 1, $smp, @collate);

#remove the split files
foreach $sp (@splitFiles) {
  if (! $alreadySplit) {
	 	`rm $workingDir/$sp`;
  }
  `rm $workingDir/$sp.spect`;
  `rm $workingDir/$sp.spect.sorted`;
}
