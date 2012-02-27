#!/usr/bin/env perl

############################################################################
# Title:          Assemble.pl
# Author:         Mark Chaisson, Glenn Tesler
# Created:        2008
# Last modified:  03/04/2010
# 
# Copyright (c) 2008-2010 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
############################################################################


BEGIN {
		unshift(@INC, $ENV{"EUSRC"} . "/assembly/");
}

use strict;
use POSIX;
use RunCmd;

my $VERSION = "2.0.1-beta";
my $RELEASE_DATE = "August 28, 2010";


my $machtype = $ENV{"MACHTYPE"};
my $srcDir = $ENV{"EUSRC"};
my $exeDir = "$srcDir/assembly/$machtype";

if (!defined $srcDir || !defined $machtype) {
		print "ERROR: Need to define environment variables EUSRC and MACHTYPE\n";
		print_usage();
		exit(1);
}

$Usage::exeDir = $exeDir;

if ($#ARGV < 1) {
		print_usage();
		exit(1);
}


my $cmdLine = join " ", ($0, @ARGV);
my $readsFile  = shift @ARGV;
my $vertexSize = shift @ARGV;
my $edgeSize = $vertexSize + 1;
my $ruleFile = "";
my $minMult = 0;
my $opt;
my $fixErrors = 1;
my $onlyFixErrors = 0;
my $numJobs = 1;
my $numJobsOpt = "";
my $verbose = 0;
my $verboseOpt = "";
my $earlyTrim = 0;
my $earlyTrimOpt = "";
#my $minCov = 0;
#my $minCovOpt = "";
my $forceLarge = 0;
my $useAltEdges = 0;
my $script = 0;
my $useDebug = 0;
my $highCov = -1; # -1: not yet specified; 0: low; 1: high

if ($#ARGV >= 0 && $ARGV[0] !~ /^\-/) {
		# Allow rules file as an optional 3rd command line argument
		$ruleFile = shift @ARGV;
}

while ($#ARGV >= 0) {
		$opt   = shift @ARGV;
		if ($opt =~ /^\-/) {
				if ($opt eq "-minMult") {
						$minMult = shift @ARGV;
				}
				elsif ($opt eq "-noFixErrors") {
						print "skipping fix errors\n";
						$fixErrors = 0;
				}
				elsif ($opt eq "-onlyFixErrors") {
						$onlyFixErrors = 1;
				}
				elsif ($opt eq "-ruleFile") {
						$opt = shift @ARGV;
						$ruleFile = $opt;
				}
				elsif ($opt eq "-numJobs") {
						$numJobs = shift @ARGV;
						$numJobsOpt = " -numJobs $numJobs ";
				}
				elsif ($opt eq "-verbose") {
						$verbose = 1;
						$verboseOpt = " -verbose";
				}
				elsif ($opt eq "-script") {
						$script = 1;
				}
				elsif ($opt eq "-earlyTrim") {
						$earlyTrim = 1;
						$earlyTrimOpt = " -earlyTrim";
				}
				elsif ($opt eq "-highCov") {
						$highCov = 1;
				}
				elsif ($opt eq "-noHighCov") {
						$highCov = 0;
				}
#				elsif ($opt eq "-minCov") {
#						$minCov = shift @ARGV;
#						$minCovOpt = " -minCov $minCov";
#				}
				elsif ($opt eq "-useAltEdges") {
						$useAltEdges = 1;
				}
				elsif ($opt eq "-noUseAltEdges") {
						$useAltEdges = 0;
				}
				elsif ($opt eq "-forceLarge") {
						$forceLarge = 1;
				}
				elsif ($opt eq "-debug") {
						$useDebug = 1;
				}
				else {
					print "bad option: $opt\n";
					print_usage();
					exit(1);
				}
		}
}

# -highCov: determine what to do if it wasn't specified
if ($highCov == -1) {
		my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,
		 $atime,$mtime,$ctime,$blksize,$blocks)
				= stat($readsFile);

		# TODO: determine appropriate threshold or other criteria to
		# decide when to use -highCov
		if ($size >= 0x40000000) {
				$highCov = 1;
		} else {
				$highCov = 0;
		}
}

#
# If the debug code is specified, change the directory for the executables.
#

if ($useDebug != 0) {
		$exeDir = $exeDir . "_debug";
		$Usage::exeDir = $exeDir;
}

###############################################################################
# Validate vertex size
###############################################################################

# Get maximum K allowed
my $maxK = 0;
my $runMaxK = "$exeDir/maxK";
if (defined $exeDir && -e $runMaxK && -x $runMaxK) {
    # TODO: probably should verify that all executables are available

    # WARNING:
    # Use RunCmd::RunCommandLog for all other system commands unless
    # there's a good reason.
		$maxK = int(`$runMaxK`);
}
if ($maxK == 0) {
		print <<"END_error";
ERROR: executable $runMaxK missing or not working.
Need to define environment variables EUSRC and MACHTYPE correctly,
and need to build EULER-SR before running this script.

END_error
		print_usage();
		exit(1);
}
if ($vertexSize > $maxK
		&& !$forceLarge) {
		print <<"END_error";
ERROR: EULER-SR has been compiled with a maximum vertex size $maxK.
To use a larger vertex size, we strongly recommend setting compilation
options for a larger maximum.  However, you may also run this script
with the '-forceLarge' option, which uses an alternate, less effficient,
set of code.

END_error
		print_usage();
		exit(1);
}

###############################################################################
# Set up various files
###############################################################################



# Option to generate a script instead of executing commands.
# Most steps will be the same, except:
# * Computing minMult with estErrorDist
# * running assemblesec.pl
if ($script) {
		$verbose = -1;
		my $cmd = $cmdLine;
		$cmd =~ s/\s-script\s/ /g;
		$cmd =~ s/\s-script$/ /g;
		print <<"START_script1";
#!/bin/sh
#
START_script1

		print_copyright($script);

		print <<"START_script2";
# Commands generated by Assemble.pl with these options:
# $cmd

START_script2
}

if (!$script) {
		print_copyright($script);
}

#
# Report file in case of fatal errors
#
$RunCmd::reportfile = "${readsFile}.report";

###############################################################################
# Generate commands to run the different parts of the assembly
###############################################################################


#
# Cleanup in case of re-running in same directory
#
my $cleanUpReports = "rm -f ${readsFile}.report fixed/${readsFile}.report simple/${readsFile}.report transformed/${readsFile}.report matetransformed/${readsFile}.report";
RunCmd::RunCommandLog($cleanUpReports, "Clean up residual of any previous run", $verbose);

#
# Initial log entries
if (!$script) {
		RunCmd::OnlyLogCommand($cmdLine, "Assembly parameters", 1, 1);
	}



#
# Commands to correct errors.
#

my $fixErrorsCmd = "";
if ($fixErrors) {

#		my $countSpectCmd =
#		"mkdir -p fixed; $exeDir/integralCountSpectrum $readsFile $vertexSize fixed/$readsFile.spect -printCount -binary";
		my $countSpectCmd = 
				"mkdir -p fixed; $exeDir/readsToSpectrum $readsFile $edgeSize fixed/$readsFile.spect -printCount -minMult $minMult; $exeDir/sortIntegralTupleList fixed/$readsFile.spect -printCount";

		my $estErrDistCmd = 	
				"$exeDir/estErrorDist fixed/$readsFile.spect -binary";

		RunCmd::RunCommandLog($countSpectCmd, "Count spectrum", $verbose);  # 1
		if ($minMult == 0) {

				if (!$script) {
						$minMult = RunCmd::RunCommandLog($estErrDistCmd, "Estimate minimum multiplicity", $verbose, 1);
						chomp $minMult;
						if ($minMult <= 1) {
								$minMult = 2;
						}
						print "Using minimum multiplicity: $minMult\n";
				} else {
						# Generating a command-line script.
						# Need to mimic the above steps through a script
						print <<"END_emulate";

# When -minMult is not specified for Assemble.pl, it normally would
# invoke estErrorDist to get the value of $minMult and then
# plug it into all further commands as needed.  We emulate that here:
minMult=\`$estErrDistCmd\`
if \\\[ \$minMult -le 1 \\\] ; then minMult=2 ; fi

END_emulate

						$minMult = '$minMult';
				}
		}


#### TEMP DISABLE MINMULT (sitmm0)
		my $sortSpectCmd =
				"$exeDir/sortIntegralTupleList fixed/$readsFile.spect -printCount -minMult $minMult";
#		my $sortSpectCmd =
#				"$exeDir/sortIntegralTupleList fixed/$readsFile.spect -printCount";

		RunCmd::RunCommandLog($sortSpectCmd, "Sort spectrum", $verbose);   # 2
		$fixErrorsCmd  = "$exeDir/fixErrors $readsFile fixed/$readsFile.spect $edgeSize fixed/$readsFile -minMult $minMult -maxScore 3 -startScore 2 -stepScore 1 -minVotes 2 -edgeLimit 3 -replaceN $numJobsOpt";
		
		
}
else {
		$fixErrorsCmd = "mkdir -p fixed; cd fixed; ln -sf ../$readsFile ./$readsFile";
}

RunCmd::RunCommandLog("$fixErrorsCmd","Fix errors",$verbose); # 4


if ($onlyFixErrors) {
		exit(0);
}

my $aMinMult = $minMult;
if ($aMinMult eq '$minMult') {
		# for -script
		$aMinMult = '\$minMult';
}
my $assembleCmd =
		"$srcDir/assembly/assemblesec.pl fixed/$readsFile -vertexSize $vertexSize $verboseOpt $earlyTrimOpt";
if ($highCov) {
		$assembleCmd .= "  -minCov $aMinMult -minKmerCount $aMinMult";
}
if ($useAltEdges) {
		$assembleCmd .= " -skipIntervals";
}
if ($forceLarge) {
		$assembleCmd .= " -forceLarge";
}
if ($useDebug) {
		$assembleCmd .= " -debug";
}

if (!$script) {
		RunCmd::RunCommandLog($assembleCmd,"Perform assembly",1);   #5  # verbose for this command even if not in general
	} else {
			# WARNING:
			# Use RunCmd::RunCommandLog for all other system commands unless
			# there's a good reason.
			my @lines = `$assembleCmd -script`;
			my $assembleScript = join "", @lines;
			print << "END_assembly";

############################################################
$assembleScript
############################################################
# Back to commands generated by Assemble.pl

END_assembly
	}


my $minEdgeLength = $vertexSize * 4;
my $bulgeLength   = $vertexSize * 4;

my $simplifyCmd;

if (!$useAltEdges) {
		# without AltEdges
		$simplifyCmd = 
				"$exeDir/simplifyGraph fixed/$readsFile simple/$readsFile -minEdgeLength $minEdgeLength -removeBulges $bulgeLength -removeLowCoverage 5 3";
} else {
		# with AltEdges
		$simplifyCmd = 
				"$exeDir/simplifyGraph fixed/$readsFile simple/$readsFile -minEdgeLength $minEdgeLength -removeBulges $bulgeLength -skipIntervals";
}

my $cpReadsToSimple =
		"mkdir -p simple; ln -sf ../$readsFile simple/$readsFile";

RunCmd::RunCommandLog($cpReadsToSimple,"Link reads to simple/ directory",$verbose);

my $printSimplifyGVZCmd =
		"$exeDir/printComponentsToGVZ simple/$readsFile";

RunCmd::RunCommandLog($simplifyCmd,"Simplify graph",$verbose);   #6

# temporarily disable
# RunCmd::RunCommandLog($printSimplifyGVZCmd,"Printing simplified graph to GraphViz",$verbose);

if ($useAltEdges) {
		my $buildEdgeOvpCmd = 
				"$exeDir/integralEdgesToOverlapList simple/$readsFile.edge $vertexSize simple/$readsFile.iovp";
		my $buildAltEdgeOvpCmd = 
				"$exeDir/integralEdgesToOverlapList simple/$readsFile.altEdges $vertexSize simple/$readsFile.altEdges.iovp";
		my $printReadIntervalsCmd =
				"$exeDir/integralPrintReadIntervals simple/$readsFile -altEdgeOvp simple/$readsFile.altEdges.iovp";

		RunCmd::RunCommandLog($buildEdgeOvpCmd,"Build edges",$verbose);
		RunCmd::RunCommandLog($buildAltEdgeOvpCmd,"Build alternative edges",$verbose);
		RunCmd::RunCommandLog($printReadIntervalsCmd,"Print read intervals",$verbose);
}


#
# Customize read transformation.
# 
my $erodePathsOpt = "";
if ($ruleFile eq "") {
  $erodePathsOpt = " -erodeShortPaths $vertexSize ";
}

my $transformCmd = 
		"mkdir -p transformed; $exeDir/transformGraph simple/$readsFile transformed/$readsFile -minPathCount 3 -notStrict $erodePathsOpt";

RunCmd::RunCommandLog($transformCmd,"Transformation based on read paths",$verbose);  #7

if ($ruleFile ne "") {
		my $buildMateListCmd = 
				"$exeDir/buildMateTable fixed/$readsFile $ruleFile $readsFile.mates";
		my $mateTransformCmd = 
				"mkdir -p matetransformed; $exeDir/mateTransformGraph transformed/$readsFile $readsFile.mates $ruleFile matetransformed/$readsFile -notStrict -minMatePairCount 5 -onlyMatePaths -scaffold";
		my $printMateContigsCmd =
				"$exeDir/printContigs matetransformed/$readsFile; cp matetransformed/$readsFile.contig .";
		
		RunCmd::RunCommandLog($buildMateListCmd,"Build list of mate pairs",$verbose); #8 
		RunCmd::RunCommandLog($mateTransformCmd,"Transformation based on mate pairs",$verbose); #9
		RunCmd::RunCommandLog($printMateContigsCmd,"Print contigs",$verbose);
}
else {
		my $printTransContigsCmd = 
				"$exeDir/printContigs transformed/$readsFile; cp transformed/$readsFile.contig .";
		RunCmd::RunCommandLog($printTransContigsCmd,"Print contigs",$verbose);
}

my $summarizeContigsCmd =
		"$srcDir/assembly_utils/SummarizeContigs.pl $readsFile.contig > $readsFile.summary";
if (!$script) {
		RunCmd::OnlyLogCommand($summarizeContigsCmd, "Summarize contigs", 1, 0);
}
RunCmd::RunCommandLog($summarizeContigsCmd,"Summarize contigs",$verbose);


sub print_copyright {
		my ($script) = @_;
		
		my $copyright = <<"END_copyright";
EULER-SR
Copyright (c) 2007-2010 The Regents of the University of California
All Rights Reserved

VERSION: $VERSION
RELEASE DATE: $RELEASE_DATE

END_copyright

		if ($script) {
				$copyright = "# ${copyright}";
				$copyright =~ s/\n/\n\# /gs;
				$copyright =~ s/\# $//;
		}
		print $copyright;

}


sub print_usage {

		my $exeDir = $Usage::exeDir;
		my $runMaxK = "$exeDir/maxK";
		my $maxK = 0;
		if (defined $exeDir && -e $runMaxK && -x $runMaxK) {
				# WARNING:
				# Use RunCmd::RunCommandLog for all other system commands unless
				# there's a good reason.
				$maxK = int(`$runMaxK`);
		}
		if (!$maxK) {
				$maxK = "[NOT AVAILABLE]";
		}


		print <<"END_usage";
usage: Assemble.pl readsFile vertexSize [-ruleFile file] [options]
   [-ruleFile file]  File with regexps for mate paired reads

   [-minMult m]      Minimum multiplicity to keep a k-mer (vertex)
                     or (k+1)-mer (edge), depending on the stage of EULER.

   [-noFixErrors]    Skip error correction
   [-onlyFixErrors]  Only do error correction, no assembly

   [-numJobs n]      For parallel error correction.

   [-useAltEdges]    Use new Alternative Edges method (IN DEVELOPMENT)
   [-noUseAltEdges]  Use original interval-based method

   [-earlyTrim]      Trim short edges of coverage 1 during initial
                     construction of deBruijn graph.

   [-highCov]        Very high coverage reads.
                     Simplify graph should remove low coverage edges
                     that were not fixed by error correction.

   [-noHighCov]      Simplify graph should not remove such edges.
                     Default (subject to change):
                     use -highCov for readsFile size >= 1 GB

TROUBLESHOOTING OPTIONS:

   [-script]         Do not execute commands.  Instead, output a
                     script that would execute the commands.

   [-verbose]        Show output from subprocesses.
                     Output is suppressed without this option.

   [-debug]          Run the debug version of the code, compiled by
                     'make debug'.

   [-forceLarge]     EULER has been compiled to handle vertex size <= ${maxK}.
                     It can be configured at compile-time to handle larger
                     vertex size, and we strongly recommend doing this
                     if you need a larger size.

                     Previously, EULER used different programs for
                     small vertex size vs. large vertex size.  To run
                     the old programs for large vertex size, use '-forceLarge'.
END_usage

#   [-minCov c]       Minimum coverage to keep an edge.
#                     Used in small vertex de Bruijn construction.
#
}
