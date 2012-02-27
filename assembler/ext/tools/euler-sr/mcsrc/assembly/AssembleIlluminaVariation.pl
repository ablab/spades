#!/usr/bin/env perl

use POSIX;
use FwgLib;
use Getopt::Long;
use strict;

sub help {
    print "usage: AssembleIllumina readsFile workingDirectory \n\n";
    print "  -vertexSize v   Assemble on graph with vertices of size 'v'. This\n";
    print "                  must be less than the read length.\n\n";
    print "  -readLength r   ******Disabled for now******\n";
		print "                  Specify the read length that is used.  When this is set,\n";
    print "                  a simple heuristic will be used to attempt to solve repeats of\n";
    print "                  length   v < repeat length < r.\n\n";
    print "                  When read length is greater than 40, we will attempt to recover\n";
    print "                  full length reads that were trimmed during error correction via our.\n";
    print "                  threading method.\n\n";
		print "  -qualityFilter  Run a couple of Illumina-specific quality filters on reads prior\n";
		print "                  to assembling.  This trims reads that are low quality, and removes\n";
		print "                  poly-A, and poly-C sequences.  This typically doesn't change assembly\n";
		print "                  qualities too much, but it speeds it up.\n\n";
		print "  -maxTrim n      Trim reads by up to 'n' bases to make them solid.\n\n";
		print "  -skipEnd n      Skip the trailing 'n' bases when counting spectra.\n\n";
    print "  -tryBeta        Try running some experimental code, after producing a normal assembly.  This may: \n";
    print "                      1. Crash.\n";
    print "                      2. Improve assembly just a little.\n";
    print "                      3. Improve assembly a fair amount.\n";
    print "                      So far, 1 and 2 are more likely, but 3 happens every once in a while.\n";
    exit(0);
}

sub DetermineFileType {
		my ($fileName) = @_;

		open(IN, $fileName) or die "cannot open $fileName\n";
		my $firstLine = <IN>;
		if ($firstLine =~ /^>/) {
				return 1;
		}
		elsif ($firstLine =~ /^@/) {
				return 2;
		}
		else {
				return 0;
		}
}

my $vertexSize = 20;
my $machtype   = $ENV{"MACHTYPE"};
my $srcDir     = $ENV{"EUSRC"};
my $exeDir     = "$srcDir/assembly/$machtype";

my $matesFile = "";
my $readLength = -1;
my $tryBeta = 0;
my $qualityFilter = 0;
my $maxTrim = 0;
my $trimEnd = 0;
&help() if ($#ARGV < 1 );

GetOptions(
    'matesFile:s'  => \$matesFile,
    'vertexSize:s' => \$vertexSize,
    'readLength:s' => \$readLength,
    'tryBeta'    => \$tryBeta,
		'qualityFilter' => \$qualityFilter,
		'maxTrim:s'=> \$maxTrim,
		'trimEnd:s'=> \$trimEnd,
    'h|help:s'     => sub { &help() },
   );

my $readsFile  = shift @ARGV;
my $workingDir = shift @ARGV;

my $origReadsFile = $readsFile;

my $fastA = 1;
my $fastQ = 2;

my $fileType = &DetermineFileType($readsFile);

if ($fileType == 0) {
		print "ERROR, unable to determine the type of input file used.\n";
		print "This should be FASTA (beginning with '>') or FASTQ (beginning \n";
		print "with '@').\n";
		exit(0);
}

FwgLib::RunCommand("mkdir $workingDir");

my $maxTrimOpt = ""; # default is no trimming.
if ($maxTrim != 0) {
		$maxTrimOpt = " -maxTrim $maxTrim";
}

my $trimEndOpt = ""; #default count full reads.
if ($trimEnd != 0) {
		$trimEndOpt = " -trimEnd $trimEnd";
}

my $altReadFileCreated = 0;
if ($fileType == $fastQ) {
		FwgLib::RunCommand("$srcDir/assembly/$machtype/qualityTrimmer -fastQ $readsFile -outFasta $workingDir/$readsFile.fasta -minQual 15 -span 3");
		$readsFile = "$workingDir/$readsFile.fasta";
		$altReadFileCreated = 1;
}

if ($qualityFilter != 0) {
		FwgLib::RunCommand("$srcDir/assembly_utils/FilterReads.pl $readsFile $srcDir/assembly/readtitle.rules $readsFile.filtered -maxA 0.85 -maxC 0.85");
		$readsFile = "$readsFile.filtered";
		$altReadFileCreated = 1;
}

if ($altReadFileCreated == 0) {
		FwgLib::RunCommand("ln -s ../$readsFile $workingDir/$readsFile");
		$readsFile = "$workingDir/$readsFile";
}


FwgLib::RunCommand("mkdir $workingDir/Simple");
FwgLib::RunCommand("$srcDir/assembly/$machtype/integralCountSpectrum $readsFile $vertexSize $readsFile.ispect -binary -printCount $trimEndOpt");

my $cmd = "$srcDir/assembly/$machtype/estErrorDist $readsFile.ispect -binSpect";
print "$cmd\n";
my $errEqCorr = `$cmd`;
print "erreqcorr: $errEqCorr\n";
chomp $errEqCorr;
if ($errEqCorr > 4) {
		# nudge the estimate a bit to get low coverage regions.
		$errEqCorr -= 2;
}
my $fixedFile = "$readsFile.fixed";
my $catchall = POSIX::ceil($errEqCorr / 2);

FwgLib::RunCommand("$srcDir/assembly/$machtype/sortIntegralTupleList $readsFile.ispect -printCount -minMult 2");

FwgLib::RunCommand("$srcDir/assembly/$machtype/fixErrorsIVoting $readsFile $readsFile.ispect $vertexSize $fixedFile -discardFile $readsFile.discards -minVotes 2 -minMult $errEqCorr $maxTrimOpt  -catchall $catchall ");

system("$srcDir/assembly/assemble.pl $fixedFile -vertexSize $vertexSize");
#system("$srcDir/assembly/CleanUpIntermadiates.pl $fixedFile");

# correct erroneous edges
my $minEdgeLength = POSIX::floor($vertexSize * 2.5);
my $bulgeSize     = $vertexSize * 3;
FwgLib::RunCommand("$srcDir/assembly/$machtype/simplifyGraph $fixedFile $workingDir/Simple/$origReadsFile.fixed.simple -minEdgeLength $minEdgeLength -vertexSize $vertexSize -removeLowCoverage 5.0 2");

if ($tryBeta) {
		system("mkdir $workingDir/Transformed");
		FwgLib::RunCommand("$srcDir/assembly/$machtype/transformGraph $workingDir/Simple/$origReadsFile.fixed.simple $vertexSize $workingDir/Transformed/$origReadsFile.fixed.simple.t -notStrict");
		if (-e "$workingDir/Transformed/$origReadsFile.fixed.simple.t.edge") {
				FwgLib::RunCommand("$srcDir/assembly/$machtype/printContigs $workingDir/Transformed/$origReadsFile.fixed.simple.t");
				FwgLib::RunCommand("mv $workingDir/Transformed/$origReadsFile.fixed.simple.t.contig $origReadsFile.contig");
		}
		else {
			FwgLib::RunCommand("$srcDir/assembly/$machtype/printContigs $workingDir/Simple/$origReadsFile.fixed.simple");
			FwgLib::RunCommand("mv $workingDir/Simple/$origReadsFile.fixed.simple.contig $origReadsFile.contig");
		}
}
else {
		FwgLib::RunCommand("$srcDir/assembly/$machtype/printContigs $workingDir/Simple/$origReadsFile.fixed.simple");
		FwgLib::RunCommand("mv $workingDir/Simple/$origReadsFile.fixed.simple.contig $origReadsFile.contig");
}



