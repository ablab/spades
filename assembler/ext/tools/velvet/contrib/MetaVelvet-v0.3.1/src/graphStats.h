/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
#ifndef _GRAPHSTATS_H_
#define _GRAPHSTATS_H_

// This header file and the corresponding code file contain a load of 
// miscellaneous functions, many of which coded quickly and used only once 
// for reasons long forgotten since.
// Sorry for the mess ;-)

///////////////////////////////////////////////////////////////////
// Useful functions
///////////////////////////////////////////////////////////////////

// Original
double getNodeDensity(Node * node);
int * makeDummySubgraphMask(Graph * graph);
int estimated_cov_multi(Graph * graph, int * subgraphMask, double expCovMulti[100]);
void resolveRepeatOfAllSubgraphs(Graph * graph, ReadSet * reads, double expCovMulti[100],
				 boolean * dubious, boolean force_jumps, int pebbleRounds,
				 double rateChimericSubgraph, boolean discardChimericSubgraph,
				 double repeatNodeCovSD);
void resetUniqueness(Graph * graph);
// Original

void displayGraphStatistics(Graph * graph);

void displayGeneralStatistics(Graph * graph, char *filename, ReadSet * reads);

void exportLongNodeSequences(char *filename, Graph * graph,
			     Coordinate minLength);

void exportMediumNodeSequences(char *filename, Graph * graph,
			       Coordinate minLength);

IDnum readStarts(Node * node);

Coordinate readCoverage(Node * node);

IDnum strainMarkerCount(Node * node, IDnum firstStrain);

IDnum nodeMultiplicity(Node * node);

Coordinate n50(Graph * graph);

double estimated_cov(Graph * graph, char * directory);

Coordinate maxLength(Graph * graph);

boolean *removeLowCoverageNodesAndDenounceDubiousReads(Graph * graph,
						       double minCov);

void exportAMOSContigs(char *filename, Graph * graph,
		       Coordinate cutoff_length, ReadSet * reads);

IDnum usedReads(Graph * graph, Coordinate minContigLength); 

Coordinate totalAssemblyLength(Graph * graph);

void logFinalStats(Graph * graph, Coordinate minContigKmerLength, char *directory);

void exportUnusedReads(Graph* graph, ReadSet * reads, Coordinate minContigKmerLength, char* filename);

///////////////////////////////////////////////////////////////////
// Dodgy functions
///////////////////////////////////////////////////////////////////

IDnum countSinksAndSources(Graph * graph);

IDnum countTangles(Graph * graph);

IDnum countRepeats(Graph * graph);

IDnum countSNPs(Graph * graph, IDnum firstStrain, int WORDLENGTH);

void displayGraphStatisticsSelective(Graph * graph, IDnum first);

void grossErrorRemoval(Graph * graph, IDnum firstStrain);

Coordinate countCommonLength(Graph * graph, IDnum firstStrain);

IDnum countBreakpoints(Graph * graph, IDnum firstStrain);

IDnum countStrainOnlyNodes(Graph * graph, IDnum firstStrain);

Coordinate countStrainOnlyBp(Graph * graph, IDnum firstStrain);

void displayStrainOnlySequences(Graph * graph, IDnum firstStrain,
				char *inputFilename, char *filename,
				int WORDLENGTH);

void displayStrainOnlyDescriptors(Graph * graph, IDnum firstStrain);

void chainSawCorrection(Graph * graph, int minMult);

void displayBreakpoints(Graph * graph, IDnum firstStrain);

void destroyStrainSpecificIslands(Graph * graph, IDnum firstStrain);

void spotIrregularReads(Graph * graph, IDnum firstStrain,
			char *sequenceFile, char *outputFile);

void displayAlignmentToReference(Graph * graph, IDnum seqID,
				 IDnum firstStrain,
				 TightString ** sequences, int WORDLENGTH,
				 char *filename);

void removeReferenceMarkers(Graph * graph, IDnum firstStrain);

void testForBizarreMarkers(Graph * graph);

void surveyPaths(Graph * graph);

void destroyMixedReads(Graph * graph, IDnum minCoverage);

void destroySinglePoolNodes(Graph * graph);
void destroySinglePoolNodesStrict(Graph * graph);
void destroyShortTips(Graph * graph);

void destroyDisconnectedElements(Graph * graph);
void measureTangleSizes(Graph * graph, Coordinate maxLength);

void destroyEmptyNodes(Graph * graph);

void removeShortReads(Graph * graph);

Coordinate totalGraphLength(Graph * graph);

void contigStats(Node ** node, IDnum readCount);

void exportContigs(Node ** contigs, ReadSet * reads, char *filename,
		   int WORDLENGTH, int pairedReadsCount);

void removeLowCoverageNodes(Graph * graph, double minCov);
void removeHighCoverageNodes(Graph * graph, double maxCov);

void removeMissingStrain(Graph * graph, Category cat);

boolean isNatural(Graph * graph);

void searchForHallidayJunction(Graph * graph);

#endif
