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


void displayGeneralStatistics(Graph * graph, char *filename, ReadSet * reads);

void exportLongNodeSequences(char *filename, Graph * graph,
			     Coordinate minLength, ReadSet * reads, ShortLength * readLengths, IDnum minCov);

Coordinate readCoverage(Node * node);

Coordinate n50(Graph * graph);

double estimated_cov(Graph * graph, char * directory);

Coordinate maxLength(Graph * graph);

boolean *removeLowCoverageNodesAndDenounceDubiousReads(Graph * graph,
						       double minCov,
						       ReadSet * reads,
						       boolean export,
						       Coordinate minLength,
						       char *filename);

void removeLowLongCoverageNodesAndDenounceDubiousReads(Graph * graph,
						       double minCov,
						       ReadSet * reads,
						       boolean * dubious,
						       boolean export,
						       Coordinate minLength,
						       char *filename);

void removeLowCoverageReferenceNodes(Graph * graph, double minCov, double minLongCov, ReadSet * reads);

void exportAMOSContigs(char *filename, Graph * graph,
		       Coordinate cutoff_length, ReadSet * reads);

IDnum usedReads(Graph * graph, Coordinate minContigLength); 

Coordinate totalAssemblyLength(Graph * graph);

void logFinalStats(Graph * graph, Coordinate minContigKmerLength, char *directory);

void exportUnusedReads(Graph* graph, ReadSet * reads, Coordinate minContigKmerLength, char* filename);

void exportLongNodeMappings(char *filename, Graph * graph, ReadSet * reads,
			     Coordinate minLength, SequencesReader *seqReadInfo);

void removeHighCoverageNodes(Graph * graph, double maxCov, boolean export, Coordinate minLength, char * filename);

void removeLowArcs(Graph * graph, double cutoff);
#endif
