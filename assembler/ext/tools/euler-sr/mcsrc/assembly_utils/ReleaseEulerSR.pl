#!/usr/bin/env perl

############################################################################
# Title:          ReleaseEulerSR.pl
# Author:         Mark Chaisson, Glenn Tesler
# Created:        2008
# Last modified:  03/03/2010
# 
# Copyright (c) 2008-2010 The Regents of the University of California
# All Rights Reserved
# See file LICENSE for details.
############################################################################


if ($#ARGV < 1) {
		print "usage: ReleaseEulerSR.pl srcDir releaseDir\n";
		exit(0);
}

$srcDir = shift @ARGV;
$releaseDir = shift @ARGV;

$EUSRC = $srcDir;


`mkdir $releaseDir`;
`cp $EUSRC/Makefile $releaseDir/`;


`cp $EUSRC/LICENSE $releaseDir/`;
`cp $EUSRC/common.mak $releaseDir/`;
#`perl -pi -e "s/CPPOPTS = -g/CPPOPTS = -O3 -DNDEBUG/g" $releaseDir/common.mak`;
`cp $EUSRC/common.mak $releaseDir/`;
`cp $EUSRC/README.eulersr $releaseDir/`;
`cp $EUSRC/make.rules $releaseDir/`;
`cp $EUSRC/mkdep.pl $releaseDir/`;
`cp $EUSRC/mkfiles.pl $releaseDir/`;
`cp $EUSRC/CreateExecBuildCommands.pl $releaseDir/`;
`mkdir $releaseDir/assembly_utils/`;
`cp $EUSRC/assembly_utils/*.pl $releaseDir/assembly_utils/`;
`mkdir $releaseDir/lib`;
`mkdir $releaseDir/lib/align`;
`mkdir $releaseDir/lib/graph/`;
`mkdir $releaseDir/lib/hash/`;
`cp $EUSRC/lib/*.cpp $EUSRC/lib/*.h $EUSRC/lib/Makefile $releaseDir/lib`;
`cp $EUSRC/lib/graph/*.h $releaseDir/lib/graph/`;
`cp $EUSRC/lib/hash/Makefile $releaseDir/lib/hash/`;
`cp $EUSRC/lib/hash/*.cpp $releaseDir/lib/hash/`;
`cp $EUSRC/lib/hash/*.h $releaseDir/lib/hash/`;
`cp $EUSRC/lib/align/*.cpp  $releaseDir/lib/align/`;
`cp $EUSRC/lib/align/*.h $EUSRC/lib/align/Makefile $releaseDir/lib/align/`;
`cd $releaseDir; find . -type d -name  ".svn" -exec rm -rf {} \\;`;

# finally copy over the $releaseDir stuff
$dir = "$releaseDir/assembly";
$src = "$EUSRC/assembly";
`mkdir $dir`;

# copy make files and auxiliary make files.
`cp $src/Makefile $dir/`;
`cp $src/ExecList.txt $dir/`;

# copy lib files

`cp $src/DeBruijnGraph.cpp $src/ReadIntervals.cpp $src/IntervalGraph.cpp $src/GraphReader.cpp $src/Spectrum.cpp $src/ReadPaths.cpp $src/ReadMap.cpp $src/ReadPos.cpp $src/SortedTupleList.cpp $src/MapContigs.cpp $src/ContigMap.cpp $src/Tuple.cpp $src/ListSpectrum.cpp $src/StringTuple.cpp $src/HashedSpectrum.cpp $src/BitSpectrum.cpp $src/NumericTuple.cpp  $src/IntegralTuple.cpp $src/DeBruijnGraph.h $src/ReadIntervals.h $src/IntervalGraph.h $src/GraphReader.h $src/Spectrum.h $src/ReadPaths.h $src/ReadPos.h $src/SortedTupleList.h $src/MapContigs.h $src/Tuple.h $src/ContigMap.h $src/ListSpectrum.h $src/StringTuple.h $src/HashedSpectrum.h $src/BitSpectrum.h  $src/NumericTuple.h $src/IntegralTuple.h $src/IntegralTupleStatic.h $src/BEdge.h $src/BVertex.h $src/BVertex.cpp $src/ThreadPath.h $src/ThreadUtils.cpp $src/ThreadUtils.h $src/ReadMap.h $src/Vertex.h $src/Edge.h $src/MultTuple.h $src/NumericHashedSpectrum.h $src/StringMultTuple.h  $src/MateTable.h $src/MateLibrary.cpp $src/MateLibrary.h $src/PathLib.h $src/PathLib.cpp $src/PathBranch.h  $src/PrintGraphSummary.cpp  $src/RepeatSearch.cpp $src/RepeatSearch.h $src/IntegralEdgesToOverlapList.cpp $src/IntegralPrintReadIntervals.cpp $src/PathBranch.cpp $src/CompareAssemblies.cpp $src/PrintGraph.cpp $src/PrintGraph.h $src/RuleList.h $src/RuleList.cpp $src/Voting.cpp $src/Voting.h $src/SAP.cpp $src/SAP.h $src/Scaffold.cpp $src/Scaffold.h $src/Trace.h $src/FixErrorsStats.h $src/VoteUtils.h $src/VoteUtils.cpp $src/VectorHashedSpectrum.h $src/ReadsToSpectrum.cpp $src/ElandToFastq.cpp  $src/CreateMateScaffold.cpp $src/AlternativeEdge.h $src/AlternativeEdge.cpp $dir/`; 

# copy utilities
`cp $src/assemblesec.pl $src/Assemble.pl $src/FixErrors.pl $src/RunCmd.pm $dir/`;

# copy data
#`cp $src/readtitle.rules $dir`;

# copy executables
`cp $src/CountSpectrum.cpp $dir/`;
`cp $src/SortVertexList.cpp $dir/`;
`cp $src/DeBruijn.cpp $dir/`;
`cp $src/IntegralCountSpectrum.cpp $dir/`;
`cp $src/SortIntegralTupleList.cpp $dir/`;
`cp $src/SmallVertexDeBruijn.cpp $dir/`;
`cp $src/EdgesToOverlapList.cpp $dir/`;
`cp $src/PrintReadIntervals.cpp $dir/`;
`cp $src/PrintContigs.cpp $dir/`;
`cp $src/EstimateErrorDistribution.cpp $dir/`;
`cp $src/FixErrorsVoting.cpp $dir/`;
`cp $src/FixErrorsSAP.cpp $dir/`;
`cp $src/CountIntegralTuples.cpp $dir/`;
`cp $src/BinSpectToAscii.cpp $dir/`;
`cp $src/SortTupleList.cpp $dir/`;
`cp $src/SimplifyGraph.cpp $dir/`;
`cp $src/ReorderIntervals.cpp $dir`;
`cp $src/JoinSourcesAndSinks.cpp $dir`;
`cp $src/QualityTrimmer.cpp $dir`;
`cp $src/SplitLinkedClones.cpp $dir`;
`cp $src/FilterIlluminaReads.cpp $dir`;
`cp $src/TransformGraph.cpp $dir`;
`cp $src/FixErrorsI.cpp $dir`;
`cp $src/SFF2Fasta.cpp $dir`;
`cp $src/FilterFailedEndReads.cpp $dir`;
`cp $src/PrintVariants.cpp $dir`;
#`cp $src/ThreadReads.cpp $dir`;
`cp $src/ThreadReads2.cpp $dir`;
`cp $src/PrintMateLengthDistribution.cpp $dir`;
`cp $src/CleanGraphWithMates.cpp $dir`;
`cp $src/GraphExplorer.cpp $dir`;
`cp $src/MateTransformGraph.cpp $dir`;
`cp $src/FilterValidMatePaths.cpp $dir`;
`cp $src/BuildMateTable.cpp $dir`;
`cp $src/LastChanceReads.cpp $dir`;
`cp $src/PrintComponentsToGVZ.cpp $dir`;
`cp $src/MaxK.cpp $dir`;


# 
# Copy over some test data.
# 

`cp $srcDir/reads.fasta $releaseDir/`;
`cp $srcDir/reads.variants.fasta $releaseDir/`;
`cp $srcDir/readtitle.rules $releaseDir/`;





#
# No more old euler
#

`tar zcvf $releaseDir.tgz $releaseDir`;

