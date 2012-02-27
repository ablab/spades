#!/usr/bin/env perl
use POSIX;
use FwgLib;

# 
# Helper script for fixing errors.  This takes as input a reads file, 
# a compressed suffix array of the original genome (this 
# assumes errors are being fixed when the reference sequence is known), 
# a directory to perform all the work in (it will be called directory.pid),
# a directory to store the results


if ($#ARGV < 3) {
		if ($#ARGV == -1) {
				PrintUsage();
				exit(1);
		}
		$opt = shift @ARGV;
		if ($opt eq "-verbose") {
				PrintUsageVerbose();
				exit(0);
		}
		else {
				PrintUsage();
				exit(0);
		}
}

$argvStr = join(" ", @ARGV);

$readsFile  = shift @ARGV;
$workingDir = shift @ARGV;
$resultsDir = shift @ARGV;

@minMult    = ();
@maxScore   = ();
@startScore = ();
@step       = ();
$grid          = 0;
$readsPerFile = 200000;
$vertexSize   = 15;
$discardFile  =  "";
$fixProg      =  "";
$nJobs        = -1;
$nProc        = -1;
$smp          = 0;
$maxTrim      = 1;
$maxGap       = 4;
$minVotes			= -1;
$minVotesOpt  = "";
$prog            = "sap";
$progOpt         = "-prog sap";
$readsPerFileOpt = "-readsPerFile $readsPerFile";
$vertexSizeOpt   = "-vertexSize $vertexSize";
$discardFileOpt  = "";
$maxScoreOpt     = "";
$maxGapOpt       = "";
$searchOpt       = "";
$logOpt          = "";
$maxScoreOpt     = "";
$curDir  = $ENV{"PWD"};
while ($#ARGV >= 0) {
  $opt = shift @ARGV;
	if ($opt eq "-itMult" ) {
		push @minMult, shift @ARGV;
		# configure optional arguments
		if ($#ARGV >= 0) {
			if ($ARGV[0] eq "-itOpt") {
				shift @ARGV;
				push @itOptions, shift @ARGV;
			}
			else {
				push @itOptions, "";
			}
		}
  }
  elsif($opt eq "-itStart") {
		push @startScore, shift @ARGV;
  }
  elsif($opt eq "-itStep") {
		push @step, shift @ARGV;
  }
  elsif($opt eq "-itMax") {
		push @maxScore, shift @ARGV;
	}
  elsif($opt eq "-curDir" ) {
    $curDir = shift @ARGV;
  }
  elsif($opt eq "-maxGap") {
		$maxGap = shift @ARGV;
		$maxGapOpt = "-maxGap $maxGap";
 	}
	elsif($opt eq "-maxScore" ) {
		$maxScore = shift @ARGV;
		$maxScoreOpt = "-maxScore $maxScore";
	}
	elsif($opt eq "-minVotes" ) {
		$minVotes = shift @ARGV;
		$minVotesOpt = "-minVotes $minVotes";
	}
	elsif($opt eq "-search") {
		$search = shift @ARGV;
		$searchOpt = "-search $search";
	}
	elsif($opt eq "-csa" ) {
		$csaFile    = shift @ARGV;
	}
  elsif($opt eq "-grid") {
	  $grid = 1;
  }
  elsif($opt eq "-readsPerFile") {
    $readsPerFile = shift @ARGV;
    $readsPerFileOpt = "-readsPerFile $readsPerFile";
  }
  elsif($opt eq "-vertexSize") {
    $vertexSize = shift @ARGV;
	  $vertexSizeOpt = "-vertexSize	$vertexSize";
  }
  elsif($opt eq "-discardFile") { 
    $discardFile = shift @ARGV;
	  $discardFileOpt = "-discardFile $discardFile";
  }
  elsif($opt eq "-prog" ) {
    $prog = shift @ARGV;
		$progOpt = "-prog $prog";
  }
  elsif($opt eq "-nJobs" ){
    $nJobs = shift @ARGV;
    $nJobsOpt = "-nJobs $nJobs";
  }
  elsif($opt eq "-nProc" ){ 
    $nProc = shift @ARGV;
 	  $nProcOpt = "-nProc $nProc";
  }
  elsif($opt eq "-smp") {
			$smp = 1;
			$smpOpt = "-smp";
  }
  elsif($opt eq "-span") {
	 $span = shift @ARGV;
   $spanOpt = "-span $span";
  } 	
	elsif($opt eq "-edgeLimit") {
		$edgeLimit = shift @ARGV;
    $edgeLimitOpt = "-edgeLimit $edgeLimit";
  }
  elsif($opt eq "-maxTrim") {
	 	$maxTrim = shift @ARGV; 
		$maxTrimOpt = "-maxTrim $maxTrim ";
  }
	elsif($opt eq "-verbose") {
			PrintUsageVerbose();
			exit(0);
	}
	else {
		PrintUsage();
		print "Bad option: $opt\n";
		exit(1);
	}

}

CheckCommandSanity();


# configure the parallel grid options
  # Create the logging and results files
FwgLib::InitDir($resultsDir);

$logFile = "$resultsDir/fixerror.$$.log";
$resultFile= "$resultsDir/fixerror.$$.result";

open(LOG, ">$logFile") or 
	die "cannot open $logFile\n";
open(RESULT, ">$resultFile") or 
	die "cannot open $resultFile\n";

#perform all operations in $workingDir

print LOG "Creating working dir $workingDir.\n";
$localReadsFile = FwgLib::SetupWorkingDir($workingDir, $readsFile);
print LOG "Using reads file $readsFile, reads in working dir: $localReadsFile\n";
print LOG "Currently at: $curdir .\n";
print LOG "Running: $argvStr .\n";
$logOpt = "-log $logFile";

#link the input reads file so that it works on the first iteration
$iter = 0;
#print "linking $readsFile to $localReadsFile.$iter\n";
`ln -sf $localReadsFile ./$localReadsFile.$iter`;

# log some interesting things
print LOG "local reads: $localReadsFile\n";

$paramStr  = "";
$nparamSets = scalar @minMult;
print LOG "Iterating $nparamSets param sets\n";
$resultStr = "";
$mach    = FwgLib::CrucialGetEnv("MACHTYPE");
$origDir = FwgLib::CrucialGetEnv("PWD");
$srcDir  = FwgLib::CrucialGetEnv("EUSRC");

# Iteratively fix the errors.  
# Each iteration is processed in the directory "workingdir.$i"


if ($prog eq "sap") {
  $fixExe = "$srcDir/assembly/$mach/fixErrorsISAP ";
}
elsif ($prog eq "vote") {
	$fixExe = "$srcDir/assembly/$mach/fixErrorsIVoting ";
}


for ($i = 0; $i <= $#minMult; $i++) {
	#
	# Each iteration fixes the file $localReadsFile.$iter
	# The result should be $localReadsFile.$iter.fixed,
	# and $localReadsFile.$iter.discards
	
  # Configure options that all error correction routines use
  $curReadsFile = "$localReadsFile.$iter";
  $minMultOpt = "-minMult $minMult[$i]";
	$optStr     = "$itOptions[$i]";
  # Initialize options that only some ec routines use
  $startScoreOpt = "";
  $stepScoreOpt  = "";
  $maxScoreOpt   = "";
  $edgeLimitOpt  = "";
  if ($prog eq "sap") {
    $startScoreOpt = "-startScore $startScore[$i]";
		$stepScoreOpt  = "-stepScore $step[$i]";
    $maxScoreOpt   = "-maxScore $maxScore[$i]";
  }

	$discardFile    = "$curReadsFile.discards";
	$discardFileOpt = "-discardFile $discardFile";
  
	if ($grid == 0 && $smp == 0) {
    # Do all the processing in just one machine.  Eventually this will all
	  # be handled by grid scripts that are aware they are on just one machine
	  # (in particular one machine that's multi-core .
		print LOG "Iteration: $iter\n";	
		# Count tuple frequencies
		$cmd = "$srcDir/assembly/CountSpectrum.pl $curReadsFile " .
	         "  $vertexSize $curReadsFile.spect";
		RunCmd($cmd);
		# Make the frequency list query-able
#		$cmd = "$srcDir/assembly/$mach/sortTupleList $curReadsFile.spect " .
#	         " $curReadsFile.spect.sorted -minMult 2";
		RunCmd($cmd);
    # Fix errors for this iteration
  	$cmd = "$fixExe " .
			# The first 4 arguments are required
	    " $curReadsFile $curReadsFile.spect  $vertexSize $curReadsFile.fixed " .
			# Split the input file
		  " -discardFile $curReadsFile.discards ".
			# Options for both SA and voting
			" $minMultOpt " . 
			# Options that are configured for spectral alignment
	    "  $startScoreOpt $stepScoreOpt $maxScoreOpt " .
			# Random options that may be changed for each step
			" $optStr " . 
			# Options for spectral alignment.
			" $edgeLimitOpt $maxTrimOpt $maxGapOpt $spanOpt "  .
			# Options that are configured for voting 
			" $minVotesOpt $searchOpt ";
		
	  RunCmd($cmd);
	$dirContents = `ls`;
	print LOG "after fixing, dir:: $dirContents\n";
		
	}
	else {
    # Run the command on the grid. 
    $gridCmd = "$srcDir/assembly/FixErrorsGrid.pl $curReadsFile " .
	  " $workingDir/grid -vertexSize $vertexSize $discardFileOpt " .
	  " -curDir $workingDir $nJobsOpt $nProcOpt $smpOpt " .
		" $progOpt $minVotesOpt $searchOpt $minMultOpt $readsPerFileOpt " .
	 	" $startScoreOpt $stepScoreOpt $maxScoreOpt $edgeLimitOpt " .
		" $maxGapOpt $maxTrimOpt " . 
		" $logOpt " . # optional file for logging
		" $itOpt "; # optional custom arguments
    print "running grid: $gridCmd\n";
		print LOG "running grid cmd: $gridCmd\n";
	  system($gridCmd);
    #RunCmd($gridCmd);
	}

	# After correcting errors, an iteration is finished.  Prepare for 	
	# the beginnign of this loop.
  $iter++;
	$dirContents = `ls`;
	print LOG "dir contents: $dirContents\n";
  `cat $curReadsFile.fixed $curReadsFile.discards > $localReadsFile.$iter`;

	# Evaluate the correction, if possible
	if ($csaFile != "") {
		$cmd = "$srcDir/bbbwt/$mach/countSeqMatches $localReadsFile.$iter " .
					 " $csaFile > $localReadsFile.$iter.nmatches";  
		RunCmd($cmd);
		#system($cmd);
	
		$iterRes = `cat $localReadsFile.$iter.nmatches`;
		chomp $iterRes;
		$resultStr .= "$iterRes, ";
	}
  $paramStr .= " $minMultOpt $sstartScoreOpt $stepScoreOpt $maxScoreOpt ";
}


# All done processing.  Move the data to where it was expected
`mv $curReadsFile.fixed $origDir/$readsFile.fixed`;
`mv $curReadsFile.discards $origDir/$readsFile.discards`;
print LOG "Results: \n";
print LOG "  $origDir/$readsFile.fixed\n";
print LOG "  $origDir/$readsFile.discards\n";

if ($csaFile ne "") {
	$cmd ="$srcDir/bbbwt/$mach/countSeqMatches $localReadsFile.fixed $csaFile > $localReadsFile.nmatches";
	print "running $cmd\n";
	RunCmd($cmd);
	$res = `cat $localReadsFile.nmatches`;
	print "ran counter $cmd\n";
	chomp $res;
	$resultStr .= " $res";
	$resultStr .= ", $paramStr";
	print RESULT "$resultStr\n";
}
close RESULT;
#`gzip $localReadsFile.fixed`;
#`cp $localReadsFile.fixed $resultsDir/$localReadsFile.fixed.$resultsTag`;

close LOG;

sub RunCmd {
  my ($cmd) = @_;
	print LOG "running $cmd\n";
  $output = `$cmd`;
  $status = $?;
  print LOG "$status $cmd\n";
  print LOG "output:\n";
  print LOG "$output";
  if ($status != 0) {
    $host = $ENV{"HOST"};
    $dir  = $workDir;
    print LOG "$jobID FAILED, $status at $host, $dir, command $cmd : $argstr\n";
    close LOG;
    close RESULT;
    exit(0);
  }
}

sub PrintUsage() {
  print "usage: FixErrorsSerialParam.pl reads workingDir logDir \n";
	print "                           -itMult mult -itStart start -itStep step -itMax max \n";
  print "    -itMult minMultiplicity to consider a tuple solid\n";
  print "    -itStart start Multiplicity for iterative error correction\n";
  print "    -itStep  progress from 'start' to 'max' at intervals of 'step'\n";
  print "    -itMax  max score to consider or a fix\n";
  print "      Each parameter set specifies a 4-set: minMult, startScore, step, \n";
	print "      max score\n";
  print "      ResultsTag should differentiate the results of this run from other runs.\n";
	print "	     It will be stored in 'resultsTag.result'\n";
	print "\n";
	print "  -csa CSAFile.csa  - Compare fixed reads with the original genome stored in \n";
	print "                      the compressed suffix array 'CSA'\n";
	print "  -grid               Use a grid to divide the workload.  The parameters that \n";
	print "                      control the grid programs are: \n";
  print "    -readsPerFile R   Make each sub-problem contain at most R reads.\n";
  print "    -vertexSize   V   Fix on vertex size V\n";
  print "                      used, all reads will be printed to the original file.\n";
  print "    -prog [vote|sap] options \n";
  print "                      Use either voting, or dp, and supply with options.  This\n";
  print "                      should be in quotes.  Include the discard file in the \n";
	print "                      option above, since it is processed differently.\n";
  print "    -nJobs N          Submit N jobs to the grid.\n";
  print "    -nProc P          Each job submits to P processors.\n";
  print "    -smp              Run in symmetric multiprocessor mode.  This doesn't submit\n";
  print "                      to a grid, but submits a few jobs at the same time to the \n";
  print "                      local machine.  Use this for fixing erros on a multi-core \n";
  print "                      machine, like a mac pro.\n";
	print "  -minVotes V         When correcting using voting, require at least 'V' votes\n";
	print "                      in order to call a correction.\n";
	print "  -verbose            Print verbose usage, including examples.\n";
}


sub PrintUsageVerbose() { print <<ENDVERBOSE; 


FixErrorsSerialParam.pl is the front-end script to the error
correction routines.  Error correction is currently an extremely
parameter-dependent process, and so this helps to handle some
of the complexity.  

This script iteratively fixes errors in a set of reads.  The basic
idea of 'serial parameters' is to take a set of reads and
several sets of error-correction parameters, and to fix a subset of reads
with the first parameter set, a subset of the unfixed reads with the
following parameter set, and so-on.  We define accuracy as the
number of corrected reads that are actually correct, and selectivity
as the number of unfixed reads.  The first parameter set is highly
accurate, but also highly selective, and the final parameter set
should be less selective, but will sacrifice some accuracy.  By
serially applying several parameter sets to error correction, some of
the dependencies on parameters are mitigated.

Because the methods are applied in serial, a lot of files are
generated.  To keep your read directories clean, you must supply a
'working' directory, and a 'results' directory for the script.  All
temporary files are stored in the 'working' directory.  For now the
full paths to the working directory and results directory must be supplied.


Two different programs are available for fixing reads: the
spectral-alignment method for fixing errors, and a voting method.
The spectral-alignment method is good at fixing errors in reads that
are somewhat long (70+ bases), and have many indels.  It does have the
drawback that it must find one short contiguous high-quality section
of a read.  When reads are short, and may have errors in the middle
such as Solexa reads, the voting procedure is better. We are still
tuning the voting procedure, so it should be used with caution.

Only one error correction routine may be used at a time.  The routine
is specified on the command line using 

   -prog  "vote|sap"  For selecting between voting or spectral
                      alignment.  Each program has a parameter set 
											that is required which will be described later.


Both methods work by taking tuples (or k-mers) that are considered to
be incorrect, and fix reads by making the tuples in the read correct.
Correct tuples are found by counting the frequency of all tuples
within a set of reads, and picking a threshold for which all tuples
appearing with a frequency above that threshold are correct.  Typical
values for tuple-size are 15 for BACs and 20 for bacterial genomes.


Since error correction is an embarrassingly parallel process, this
script also can break down the error correction into parallel
processes, and can either submit them to a grid using the Sun Grid
Engine (SGE) platform, or can run them in parallel on a multi-core
processor (called Symmetric Multiprocessor mode, or SMP mode).




The grid options are as follows:
  -readsPerFile  N    Divide a reads file of R reads into R/N files,
                      and process each file separately.
  
  -nJobs N            When running on a grid, N jobs are submitted to
                    	the grid.  This defaults to 1 in SMP mode.
  -nProc N            Run on N processors.  This option is necessary
                    	when running in SMP mode.
  -smp                Use SMP mode.
  -grid               Submit jobs to the grid.



When fixing errors using spectral alignment, there are several
parameters that are constant in each iteration, and several that are
changed on each iteration.  

Consider a (generic) read shown by 'N's, with '*'s marking the
positions of errors.  Assume that any tuple that overlaps a * is
considered erroneous (has a low frequency). 


0         10        20        30        40        50  54
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
   *      *                   *         *      *     * 


-span   N   Attempt to fix reads where there is a stretch of N tuples
            that are all solid.  If the vertex size is 15, and span is
            6, this read cannot be fixed, but if the span were 1, the
            read may be fixed.

-edgeLimit N Spectral alignment is bad at fixing errors close to the 
            ends of reads.  If an error is detected within 'N' of the
            end of a read, the read is considered unfixable.  In the 
						read above, if 'edgeLimit' is 3, the read cannot be fixed.

-maxTrim    One way to get around discarding reads that have errors at
            the end of the read is to trim starting at the position of
						the error until the end of the read.  The maxTrim option
						determines the maximum amount that may be trimmed off the 
						end of a read.  This leaves higher quality edges of reads,
            but if the coverage is low in the sequencing project, a
            high value of 'maxTrim' may often trim all ends of reads
            that are low coverage, thus creating 'hard stops' in the
            read coverage.  This is often set to 1 when fixing errors
						in 454 reads.

-maxGap     The maximum number of indels to consider when trying to
						fix a read.  Increasing 'maxGap' increases runtime.
						Typically, even for high homopolymer content reads, no
						more than 4 indels are necessary.
						
The options that are used in the iterative steps for spectral alignment
						are:

-itMult     Consider tuples with frequency 'itMult' to be correct.

To save time, the number of modifications to each read are limited.
The spectral alignment method first tries to fix reads with 'itStart'
modifications, then 'itStart + itStep', ... until 'itMax' changes are
tried to fix a read.  Many reads may be fixed with just one or two
modifications, and searching below about 4 is very fast.  Each of
'itStart, itStep, and itMax' must be supplied on each error correction
iteration.  Increasing itMax decreases selectivity at an expense of
some accuracy, and fixes very low quality reads.  Sometimes the
sequencing reaction is just bad, so it is often better to not try and
correct very low quality reads.

So, for each iteration of spectral alignment error correction, the
four parameters must be given:

-itMult mult -itStart start -itStep step -itmax max

For each quartet (mult, start, step, max), error correction is ran once.

A decent set of parameters for a sequencing project that produced 20x
coverage on e-coli with 454 GS20 reads is:
15 2 2 10  
10 2 2 10  
8 2 2 10  
6 2 2 6  
3 2 2 10 


Fixing errors with voting is a more simple process that (currently)
takes a long time.  Each read may sometimes be fixed several different
ways to make it correct.  Each candidate change is a
(position,mutation) pair.  The number of 'votes' cast by a change is
the number of tuples that are made correct by a change.  

The options controlling voting are:

-itMult  M   Consider tuples with frequency above 'M' to be correct
             (same as above).

-minVotes   V   Require at least 'V' votes to allow a position to be
 fixed.

-search  S   Search for 'S' positions to fix in order to fix a read.
		         This runs in O(n^S) time, where n is the length of the
		         read, so S > 2 is VERY SLOW!



Here are some examples of error correction 


FixErrorsSerialParam.pl test.reads /data/shortassembly/test/test_errcorr/work \\
		/data/shortassembly/test/test_errcorr/result -smp \\
		-readsPerFile 2 -nJobs 10 -nProc 3 -prog sap -itMult 10 -itStart 2 \\
		-itStep 1 -itMax 4  -itMult 5 -itStart 2 -itStep 1 -itMax 5 \\
		-edgeLimit 4 

This runs spectral alignment error correction, on a SMP machine with 3
		processors.  


FixErrorsSerialParam.pl test.reads \\
		/data/shortassembly/test/test_errcorr/work \\
		/data/shortassembly/test/test_errcorr/result -readsPerFile 2\\
		-nJobs 10 -nProc 3 -prog vote -itMult 10 -itMult 3 -minVotes 3 \\
		-smp

This runs the voting procedure twice, requiring 3 votes per run.





ENDVERBOSE

exit(0);

}

sub CheckCommandSanity() {
		if ($smp != 0 and $nProc == 0) {
				print "ERROR, when running in smp mode, a number of processors must\n";
				print "       be supplied using the -nProc option\n";
		}
		if ($smp != 0) {
				$nJobs = 1;
				$nJobsOpt = "-nJobs 1";
		}
}
