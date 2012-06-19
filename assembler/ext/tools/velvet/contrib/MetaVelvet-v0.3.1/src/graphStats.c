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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "globals.h"
#include "graph.h"
#include "readSet.h"
#include "tightString.h"
#include "passageMarker.h"
#include "concatenatedGraph.h"
#include "readCoherentGraph.h"
#include "fibHeap.h"
#include "utility.h"
// Original
#include "shortReadPairs.h"
// Original

// Original
#define LEN_HISTO_X 10000
// Original

static void surveyPath(PassageMarker * marker)
{
	Coordinate length = 0;
	Coordinate realLength = 0;
	PassageMarker *current = marker;

	if (passageMarkerDirection(current) < 0)
		current = getTwinMarker(current);

	for (; current != NULL; current = getNextInSequence(current)) {
		length +=
		    getNodeLength(getNode(current)) -
		    getStartOffset(current) - getFinishOffset(current);
		if (getPassageMarkerFinish(current) > 0)
			realLength = getPassageMarkerFinish(current);
	}

	printf("SURVEY %ld %lld %lld\n", (long) getAbsolutePassMarkerSeqID(marker),
	       (long long) realLength, (long long) length);
}

void surveyPaths(Graph * graph)
{
	IDnum ID;
	PassageMarker *marker;

	for (ID = 1; ID <= nodeCount(graph); ID++)
		for (marker = getMarker(getNodeInGraph(graph, ID));
		     marker != NULL; marker = getNextInNode(marker))
			if ((passageMarkerDirection(marker) > 0
			     && isInitial(marker))
			    || (passageMarkerDirection(marker) < 0
				&& isTerminal(marker)))
				surveyPath(marker);
}

static PassageMarkerList *copyMarkers(Node * node)
{
	PassageMarkerList *list = NULL;
	PassageMarkerList *new;
	PassageMarker *currentMarker;

	for (currentMarker = getMarker(node); currentMarker != NULL;
	     currentMarker = getNextInNode(currentMarker)) {
		new = newPassageMarkerList(currentMarker, list);
		list = new;
	}

	return list;
}

static boolean removeDead(PassageMarkerList ** list)
{
	PassageMarkerList *current, *next;
	boolean removed = false;

	if (*list == NULL)
		return false;

	current = *list;

	while (current->next != NULL) {
		next = current->next;

		if (isTerminal(next->marker)) {
			removed = true;
			current->next = next->next;
			deallocatePassageMarkerList(next);
		} else
			current = current->next;
	}

	current = *list;
	if (isTerminal(current->marker)) {
		removed = true;
		*list = current->next;
		deallocatePassageMarkerList(current);
	}

	return removed;
}

static Node *chooseDestination(PassageMarkerList * list)
{
	PassageMarkerList *current = list;
	Node *destination;

	destination = getNode(getNextInSequence(current->marker));
	while (current != NULL) {
		if (getNode(getNextInSequence(current->marker)) !=
		    destination)
			return NULL;
		current = current->next;
	}

	return destination;
}

static void destroyPassageMarkerList(PassageMarkerList ** list)
{
	PassageMarkerList *ptr;

	while (*list != NULL) {
		ptr = *list;
		*list = ptr->next;
		deallocatePassageMarkerList(ptr);
	}
}

static void updateMarkers(PassageMarkerList * list)
{
	PassageMarkerList *current;

	for (current = list; current != NULL; current = current->next)
		current->marker = getNextInSequence(current->marker);
}

Coordinate computeSubsequentNodesLength(Node * node)
{
	PassageMarkerList *list;
	Node *nextNode;
	Coordinate totalLength = 0;
	boolean uncertain = false;

	list = copyMarkers(node);

	while (true) {
		if (removeDead(&list))
			uncertain = true;

		if (uncertain && simpleArcCount(node) > 1) {
			destroyPassageMarkerList(&list);
			return totalLength;
		}

		if (list == NULL)
			return totalLength;

		nextNode = chooseDestination(list);
		if (nextNode == NULL) {
			destroyPassageMarkerList(&list);
			return totalLength;
		}

		totalLength += getNodeLength(nextNode);

		updateMarkers(list);
	}

	// Impossible instruction
	return -1;
}

Coordinate computeVirtualNodeLength(Node * node)
{
	Coordinate virtualLength;

	if (node == NULL)
		return 0;

	virtualLength = getNodeLength(node);

	virtualLength += computeSubsequentNodesLength(node);
	virtualLength += computeSubsequentNodesLength(getTwinNode(node));

	return virtualLength;
}

void testForBizarreMarkers(Graph * graph)
{
	IDnum index;
	Node *node;
	PassageMarker *marker;

	for (index = 1; index < nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL)
			continue;

		for (marker = getMarker(node); marker != NULL;
		     marker = getNextInNode(marker)) {
			if (getTwinMarker(marker) == NULL) 
				exitErrorf(EXIT_FAILURE, false, "Bizarre marker %s",
				       readPassageMarker(marker));
		}
	}
}

// COunts how many nodes are dead-ends
IDnum countSinksAndSources(Graph * graph)
{
	IDnum nodeIndex;
	IDnum result = 0;
	Node *node;

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex);
		if (arcCount(node) == 0
		    || arcCount(getTwinNode(node)) == 0)
			result++;
	}

	return result;
}

// Counts how many nodes have several arcs either going in or coming out
IDnum countTangles(Graph * graph)
{
	IDnum nodeIndex;
	IDnum result = 0;
	Node *node;

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex);
		if (arcCount(node) > 1 || arcCount(getTwinNode(node)) > 1)
			result++;
	}
	return result;
}

// Counts nodes with exactly one incoming and one outgoing arc
IDnum countRepeats(Graph * graph)
{
	IDnum nodeIndex;
	IDnum result = 0;
	Node *node;

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex);
		if (arcCount(node) == 1
		    && arcCount(getTwinNode(getDestination(getArc(node))))
		    == 1)
			if (getNodeID
			    (getTwinNode(getDestination(getArc(node)))) < 0
			    ||
			    getNodeID(getTwinNode
				      (getDestination(getArc(node)))) >
			    getNodeID(node))
				result++;
	}
	return result;

}

// Counts the number of markers for one node
int nodeGenomicMultiplicity(Node * node, IDnum firstStrain)
{
	int counter = 0;
	PassageMarker *marker;

	if (node == NULL)
		return 0;

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) < firstStrain)
			counter++;

	return counter;
}

// Counts the number of markers for one node
IDnum nodeMultiplicity(Node * node)
{
	int counter = 0;
	PassageMarker *marker;

	if (node == NULL)
		return 0;

	marker = getMarker(node);
	while (marker != NULL) {
		counter++;
		marker = getNextInNode(marker);
	}
	return counter;
}

// Prints out a set of predefined statistics for one node
char *nodeStatistics(Node * node)
{
	char *s = mallocOrExit(100, char);
	sprintf(s, "NODE %ld\t%lld\t%i\t%i\t%ld", (long) getNodeID(node),
		(long long) getNodeLength(node), simpleArcCount(node),
		simpleArcCount(getTwinNode(node)), (long) nodeMultiplicity(node));
	return s;
}

// Prints out a table of statistics for all the nodes of the graph
void displayGraphStatistics(Graph * graph)
{
	IDnum nodeIndex;
	printf("NODE ID\tlgth\tFwd\tBck\tMult\n");
	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++)
		printf("%s\n",
		       nodeStatistics(getNodeInGraph(graph, nodeIndex)));
}

void displayNodeStatisticsSelective(Node * node, IDnum first)
{
	PassageMarker *marker;
	boolean originalGenome;
	boolean strain;

	if (node == NULL)
		return;

	marker = getMarker(node);

	originalGenome = false;
	strain = false;
	while (marker != NULL) {
		if (getAbsolutePassMarkerSeqID(marker) < first)
			originalGenome = true;

		if (getAbsolutePassMarkerSeqID(marker) >= first)
			strain = true;

		marker = getNextInNode(marker);
	}

	printf("%s", nodeStatistics(node));
	if (originalGenome && !strain)
		printf("\tTRUE");
	else
		printf("\tFALSE");

	if (originalGenome && strain)
		printf("\tTRUE");
	else
		printf("\tFALSE");


	if (strain && !originalGenome)
		puts("\tTRUE");
	else
		puts("\tFALSE");

}

void displayGraphStatisticsSelective(Graph * graph, IDnum first)
{
	IDnum index;

	for (index = 1; index <= nodeCount(graph); index++)
		displayNodeStatisticsSelective(getNodeInGraph
					       (graph, index), first);

}

boolean isOnlyGenome(Node * node, IDnum firstStrain)
{
	PassageMarker *marker;

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) >= firstStrain)
			return false;

	return true;
}

boolean isOnlyStrain(Node * node, IDnum firstStrain)
{
	PassageMarker *marker;

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) < firstStrain)
			return false;

	return true;
}

boolean isSNP(Node * node, IDnum firstStrain, int WORDLENGTH)
{
	IDnum sequence;
	Coordinate position;

	if (getNodeLength(node) != WORDLENGTH)
		return false;

	if (getMarker(node) == NULL)
		return false;

	if (getAbsolutePassMarkerSeqID(getMarker(node)) >= firstStrain)
		return false;

	if (getNextInNode(getMarker(node)) != NULL)
		return false;

	if (arcCount(node) != 1)
		return false;

	if (arcCount(getTwinNode(node)) != 1)
		return false;

	if (isOnlyGenome(getDestination(getArc(node)), firstStrain))
		return false;

	if (isOnlyGenome
	    (getDestination(getArc(getTwinNode(node))), firstStrain))
		return false;

	sequence = getPassageMarkerSequenceID(getMarker(node));

	if (sequence >= 0)
		position = getPassageMarkerStart(getMarker(node));
	else {
		sequence = -sequence;
		position = getPassageMarkerFinish(getMarker(node));
	}

	printf("SNP\t%lld\t%ld\n", (long long) position, (long) sequence);

	return true;
}

IDnum strainMarkerCount(Node * node, IDnum firstStrain)
{
	PassageMarker *marker;
	IDnum counter = 0;

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) >= firstStrain)
			counter++;

	return counter;
}

boolean isError(Node * node, IDnum firstStrain)
{
	return (strainMarkerCount(node, firstStrain) < 5);
}

void removeStrainMarkers(Node * node, IDnum firstStrain)
{
	PassageMarker *marker;
	PassageMarker *tmp = NULL;

	marker = getMarker(node);
	while (marker != NULL) {
		tmp = getNextInNode(marker);

		if (getAbsolutePassMarkerSeqID(marker) >= firstStrain)
			destroyPassageMarker(marker);
		marker = tmp;
	}

}

void chainSawCorrection(Graph * graph, int minMult)
{
	IDnum nodeIndex;
	IDnum removed = 0;

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		if (markerCount(getNodeInGraph(graph, nodeIndex)) <
		    minMult) {
			destroyNode(getNodeInGraph(graph, nodeIndex),
				    graph);
			removed++;
		}
	}

	printf("%d dubious nodes removed\n", removed);
	concatenateGraph(graph);
	printf("%d node in the end\n", nodeCount(graph));
}

void grossErrorRemoval(Graph * graph, IDnum firstStrain)
{
	IDnum nodeIndex;
	IDnum removed = 0;

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		if (isError(getNodeInGraph(graph, nodeIndex), firstStrain)) {
			if (isOnlyStrain
			    (getNodeInGraph(graph, nodeIndex),
			     firstStrain)) {
				destroyNode(getNodeInGraph
					    (graph, nodeIndex), graph);
				removed++;
			} else
				removeStrainMarkers(getNodeInGraph
						    (graph, nodeIndex),
						    firstStrain);
		}
	}

	printf("%d dubious nodes removed\n", removed);
	concatenateGraph(graph);
	printf("%d node in the end\n", nodeCount(graph));
}

IDnum countSNPs(Graph * graph, IDnum firstStrain, int WORDLENGTH)
{
	IDnum index;
	IDnum counter = 0;

	for (index = 1; index < nodeCount(graph); index++)
		if (isSNP
		    (getNodeInGraph(graph, index), firstStrain,
		     WORDLENGTH))
			counter++;

	return counter;
}

Coordinate commonLength(Node * node, IDnum firstStrain)
{
	PassageMarker *marker = getMarker(node);
	int orig = 0;
	int strain = 0;

	while (marker != NULL) {
		if (getAbsolutePassMarkerSeqID(marker) < firstStrain)
			orig++;
		else
			strain++;
		marker = getNextInNode(marker);
	}

	if (orig == 0 || strain == 0)
		return 0;

	return (Coordinate) orig *getNodeLength(node);
}

Coordinate countCommonLength(Graph * graph, IDnum firstStrain)
{
	IDnum index;
	Coordinate res = 0;

	for (index = 1; index <= nodeCount(graph); index++)
		res +=
		    commonLength(getNodeInGraph(graph, index),
				 firstStrain);

	return res;
}

boolean isMixed(Node * node, IDnum firstStrain)
{
	return !isOnlyStrain(node, firstStrain)
	    && !isOnlyGenome(node, firstStrain);
}

int countLocalBreakpoints(PassageMarker * marker, IDnum firstStrain)
{
	PassageMarker *localMarker;
	IDnum sequenceID = getAbsolutePassMarkerSeqID(marker);
	IDnum localSeqID;
	Coordinate start = getPassageMarkerStart(marker);
	Node *localNode = getNode(marker);
	Node *destination;
	Arc *arc;
	int arcCount = 0;
	int arcIndex;
	boolean *arcStatus;
	int counter = 0;

	if (!isMixed(localNode, firstStrain))
		return 0;

	// Count arcs
	for (arc = getArc(localNode); arc != NULL; arc = getNextArc(arc))
		arcCount++;
	arcStatus = callocOrExit(arcCount, boolean);
	// Check for other genomic markers in node
	for (localMarker = getMarker(localNode); localMarker != NULL;
	     localMarker = getNextInNode(localMarker)) {
		localSeqID = getAbsolutePassMarkerSeqID(localMarker);
		if (localSeqID >= firstStrain)
			continue;

		if (localSeqID < sequenceID)
			return 0;

		if (localSeqID == sequenceID
		    && getPassageMarkerStart(localMarker) < start)
			return 0;

		destination = getNode(getNextInSequence(localMarker));

		// Enter into table:
		arcIndex = 0;
		for (arc = getArc(localNode);
		     getDestination(arc) != destination;
		     arc = getNextArc(arc))
			arcIndex++;
		arcStatus[arcIndex] = true;
	}

	// Check other nodes
	arcIndex = 0;
	for (arc = getArc(localNode); arc != NULL; arc = getNextArc(arc)) {
		if (!arcStatus[arcIndex]
		    && isMixed(getDestination(arc), firstStrain))
			counter++;
		arcIndex++;
	}

	free(arcStatus);
	return counter;
}

IDnum countBreakpoints(Graph * graph, IDnum firstStrain)
{
	PassageMarker *marker;
	IDnum seqIndex;
	IDnum total = 0;

	for (seqIndex = 1; seqIndex < firstStrain; seqIndex++) {
		marker = getMarker(getNodeInGraph(graph, seqIndex));
		while (marker != NULL) {
			total +=
			    countLocalBreakpoints(marker, firstStrain);
			marker = getNextInSequence(marker);
		}
	}

	return total;
}

IDnum countStrainOnlyNodes(Graph * graph, IDnum firstStrain)
{
	IDnum index;
	IDnum total = 0;

	for (index = 1; index <= nodeCount(graph); index++)
		if (isOnlyStrain
		    (getNodeInGraph(graph, index), firstStrain))
			total++;

	return total;
}

Coordinate countStrainOnlyBp(Graph * graph, IDnum firstStrain)
{
	IDnum index;
	Coordinate total = 0;
	Node *node;
	Arc *arc;
	Coordinate local;

	for (index = 1; index <= nodeCount(graph); index++) {
		if (isOnlyStrain
		    (getNodeInGraph(graph, index), firstStrain)) {
			node = getNodeInGraph(graph, index);
			local = getNodeLength(node);

			for (arc = getArc(node); arc != NULL;
			     arc = getNextArc(arc)) {
				if (!isOnlyStrain
				    (getDestination(arc), firstStrain)) {
					local -= 24;
					break;
				}
			}

			for (arc = getArc(getTwinNode(node)); arc != NULL;
			     arc = getNextArc(arc)) {
				if (!isOnlyStrain
				    (getDestination(arc), firstStrain)) {
					local -= 24;
					break;
				}
			}

			if (local < 0)
				local = 1;

			total += local;
		}
	}

	return total;
}

IDnum genomeMarkerCount(Node * node, IDnum firstStrain)
{
	PassageMarker *marker;
	IDnum counter = 0;

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) < firstStrain)
			counter++;

	return counter;
}

Coordinate readCoverage(Node * node)
{
	PassageMarker *marker;
	Coordinate sum = 0;

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker)) {
		if (getTwinMarker(marker) == NULL) {
			printf("Node %d screwed up\n", getNodeID(node));
			printf("Sequence %d\n",
			       getPassageMarkerSequenceID(marker));
			abort();
		}
		sum += getPassageMarkerLength(marker);
	}

	return sum;
}

Coordinate refReadCoverage(Node * node, IDnum firstStrain)
{
	PassageMarker *marker;
	Coordinate sum = 0;

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) < firstStrain)
			sum += getPassageMarkerLength(marker);

	return sum;
}

Coordinate newReadCoverage(Node * node, IDnum firstStrain)
{
	PassageMarker *marker;
	Coordinate sum = 0;

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) >= firstStrain) {
			sum += getPassageMarkerLength(marker);
			if (getPassageMarkerLength(marker) < 0)
				printf("Bizarre marker %d at node %d\n",
				       getPassageMarkerSequenceID(marker),
				       getNodeID(node));
		}

	return sum;
}

IDnum readStarts(Node * node)
{
	PassageMarker *marker;
	IDnum sum = 0;

	if (node == NULL)
		return 0;

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker)) {
		if (getPassageMarkerSequenceID(marker) > 0
		    && isInitial(marker))
			sum++;
		else if (getPassageMarkerSequenceID(marker) < 0
			 && isTerminal(marker))
			sum++;
	}

	return sum;
}

static void printShortCounts(FILE * outfile, Node * node, Graph * graph, ReadSet * reads) {
	IDnum counts[CATEGORIES];
	Category cat;
	IDnum shortReadIndex;
	IDnum readID;
	IDnum shortReadCount;
	ShortReadMarker *array;
	ShortReadMarker *marker;

	if (!readStartsAreActivated(graph)) {
		for (cat = 0; cat < CATEGORIES; cat++)
			fprintf(outfile, "\tN/A");
		return;
	}

	shortReadCount = getNodeReadCount(node, graph);
	array = getNodeReads(node, graph);

	for (cat = 0; cat < CATEGORIES; cat++)
		counts[cat] = 0;

	for (shortReadIndex = 0; shortReadIndex < shortReadCount; shortReadIndex++) {
		marker = getShortReadMarkerAtIndex(array, shortReadIndex);
		readID = getShortReadMarkerID(marker);
		cat = reads->categories[readID - 1] / 2;
		counts[cat]++;
	}

	for (cat = 0; cat < CATEGORIES; cat++)
		fprintf(outfile, "\t%li", (long) counts[cat]);
}

void displayGeneralStatistics(Graph * graph, char *filename, ReadSet * reads)
{
	IDnum nodeIndex;
	Node *node;
	Category cat;
	FILE *outfile;

	outfile = fopen(filename, "w");
	if (outfile == NULL) {
		printf("Couldn't open file %s, sorry\n", filename);
		return;
	} else
		printf("Writing into stats file %s...\n", filename);

	fprintf(outfile, "ID\tlgth\tout\tin\tlong_cov");

	for (cat = 0; cat < CATEGORIES; cat++) {
		fprintf(outfile, "\tshort%i_cov", (int) (cat + 1));
		fprintf(outfile, "\tshort%i_Ocov", (int) (cat + 1));
	}

	fprintf(outfile, "\tlong_nb");
	for (cat = 0; cat < CATEGORIES; cat++) {
		fprintf(outfile, "\tshort%i_nb", (int) (cat + 1));
	}

	fprintf(outfile, "\n");

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex);
		if (node == NULL)
			continue;
		fprintf
		    (outfile, "%ld\t%lld\t%i\t%i",
		     (long) getNodeID(node), (long long) getNodeLength(node), arcCount(node),
		     arcCount(getTwinNode(node)));

		if (getNodeLength(node) > 0) {
			fprintf(outfile, "\t%f",
				readCoverage(node) /
				(double) getNodeLength(node));
			for (cat = 0; cat < CATEGORIES; cat++) {
				fprintf(outfile, "\t%f",
					getVirtualCoverage(node,
							   cat) /
					(double) getNodeLength(node));
				fprintf(outfile, "\t%f",
					getOriginalVirtualCoverage(node,
								   cat) /
					(double) getNodeLength(node));
			}
		} else {
			fprintf(outfile, "\tInf");
			for (cat = 0; cat < CATEGORIES; cat++)
				fprintf(outfile, "\tInf\tInf");
		}

		fprintf(outfile, "\t%li", (long) markerCount(node));
		printShortCounts(outfile, node, graph, reads); 

		fprintf(outfile, "\n");
	}

	fclose(outfile);
}

void destroyStrainSpecificIslands(Graph * graph, IDnum firstStrain)
{
	IDnum index;
	Arc *arc;
	boolean isModified = true;
	Node *node;
	IDnum counter = 0;

	resetNodeStatus(graph);

	puts("Destroying disconnected strain specific sub-graphs");

	// Mark all genomic nodes 
	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (!isOnlyStrain(node, firstStrain))
			setNodeStatus(node, true);
	}

	// Mark nodes connected to genomic nodes
	while (isModified) {
		isModified = false;

		for (index = 1; index <= nodeCount(graph); index++) {
			node = getNodeInGraph(graph, index);

			if (getNodeStatus(node))
				continue;

			for (arc = getArc(node); arc != NULL;
			     arc = getNextArc(arc)) {
				if (getNodeStatus(getDestination(arc))) {
					isModified = true;
					setNodeStatus(node, true);
				}
			}

			for (arc = getArc(getTwinNode(node)); arc != NULL;
			     arc = getNextArc(arc)) {
				if (getNodeStatus(getDestination(arc))) {
					isModified = true;
					setNodeStatus(node, true);
				}
			}
		}
	}

	// Remove all unmarked nodes
	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (!getNodeStatus(node)) {
			destroyNode(node, graph);
			counter++;
		}
	}

	// Renumber graph nodes
	printf("Removed %d nodes \n", counter);
	renumberNodes(graph);
}

void displayStrainOnlySequences(Graph * graph, IDnum firstStrain,
				char *inputFilename, char *filename,
				int WORDLENGTH)
{
	IDnum nodeIndex;
	Node *node;
	FILE *outfile = fopen(filename, "w");
	Coordinate start, finish;
	char str[100];
	TightString *tString;
	IDnum readID;
	Coordinate readCoord;

	if (outfile == NULL) {
		printf("Could not write into %s, sorry\n", filename);
		return;
	}

	destroyStrainSpecificIslands(graph, firstStrain);

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex);
		if (isOnlyStrain(node, firstStrain)
		    && getNodeLength(node) > 500) {
			tString = expandNode(node, WORDLENGTH);
			readID =
			    getPassageMarkerSequenceID(getMarker(node));
			readCoord = getPassageMarkerStart(getMarker(node));
			fprintf(outfile, "> UNIQUE SEQUENCE %ld; %lld\n",
				(long) readID, (long long) readCoord);

			start = 0;
			while (start <= getLength(tString)) {
				finish = start + 60;
				readTightStringFragment(tString, start,
							finish, str);
				fprintf(outfile, "%s\n", str);
				start = finish;
			}
		}
	}

	fclose(outfile);
}

void displayStrainOnlyDescriptors(Graph * graph, IDnum firstStrain)
{
	IDnum nodeIndex;
	Node *node;
	char *str;

	destroyStrainSpecificIslands(graph, firstStrain);

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		printf("node %d from %d\n", nodeIndex, nodeCount(graph));
		node = getNodeInGraph(graph, nodeIndex);

		if (isOnlyStrain(node, firstStrain)) {
			str = readNode(node);
			printf("> UNIQUE SEQUENCE %s\n", str);
			free(str);
		}
	}
}

void displayLocalBreakpoint(PassageMarker * strainMarker,
			    IDnum firstStrain,
			    PassageMarker * genomeMarker,
			    Node ** genomeDestination,
			    Node ** strainDestination, IDnum * counter,
			    IDnum nodeCount)
{
	boolean isTranslocation;
	PassageMarker *marker;
	Node *destination, *destinationA;
	Node *destination2, *destination2A;
	Node *node1, *node2;
	IDnum localID = getNodeID(getNode(strainMarker));

	// Eliminate genomic markers
	if (strainMarker == genomeMarker)
		return;

	destinationA = getNode(getNextInSequence(strainMarker));

	if (destinationA == NULL)
		return;

	// Eliminate those that follow some local strain
	if (isDestinationToMarker(genomeMarker, destinationA)) {
//              puts("Parallel paths");
		return;
	}

	destination2A = getNode(getNextInSequence(genomeMarker));

	if (destination2A == NULL)
		return;

	printf("Lengths %lld %lld\n", (long long) getNodeLength(destinationA),
	       (long long) getNodeLength(destination2A));

	// Hop to another genomic node
//      if (getNodeLength(destinationA) > 24) {
	//printf("wrong length %d %d\n", getNodeLength(destination) , getNodeID(destination));
//              return;
//      }

	destination =
	    getNode(getNextInSequence(getNextInSequence(strainMarker)));

	if (destination == NULL)
		return;

	// Eliminate those that point to uniquely strain sequences 
	if (nodeGenomicMultiplicity(destination, firstStrain) != 1) {
//              puts("Multiple genome reads");
		return;
	}
	// Hop to another genomic node
//      if (getNodeLength(destination2A) != 24) {
	//puts("wrong length 2");
//              return;
//      }

	destination2 =
	    getNode(getNextInSequence(getNextInSequence(genomeMarker)));

	if (destination2 == NULL)
		return;


	if (destination == destination2)
		return;

	// Eliminate those that point to uniquely strain sequences 
	if (isOnlyGenome(destination2, firstStrain))
		return;

	setSingleNodeStatus(getNode(strainMarker), true);
	strainDestination[localID + nodeCount] = destination;
	genomeDestination[localID + nodeCount] = destination2;

//      printf("Assigning %p and %p to %d\n", destination, destination2, localID);
	printf("lengths %lld\t%lld\n", (long long) getNodeLength(destinationA),
	       (long long) getNodeLength(destination2A));

	// Detect translocation
	isTranslocation = true;
	for (marker = getMarker(destination); marker != NULL;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) ==
		    getAbsolutePassMarkerSeqID(genomeMarker)) {
			isTranslocation = false;
			break;
		}

	if (isTranslocation) {
		printf("BREAK TRANS\t%ld\t%lld\t%lld\t%lld\n",
		       (long) getAbsolutePassMarkerSeqID(genomeMarker),
		       (long long) getPassageMarkerStart(genomeMarker),
		       (long long) getNodeLength(destinationA),
		       (long long) getNodeLength(destination2A));
		counter[2]++;
		return;
	}
	// Detect breakpoint
	printf("BREAK INTRA\t%ld\t%lld\t%lld\t%lld\n",
	       (long) getAbsolutePassMarkerSeqID(genomeMarker),
	       (long long) getPassageMarkerStart(genomeMarker),
	       (long long) getNodeLength(destinationA), (long long) getNodeLength(destination2A));
	counter[1]++;

	// Check for inversion
	if (getPassageMarkerSequenceID(marker) !=
	    -getPassageMarkerSequenceID(genomeMarker))
		return;

//      puts("potential!!");

	node1 = getTwinNode(destination);

	if (getNodeStatus(node1)) {
		node2 =
		    getTwinNode(genomeDestination
				[getNodeID(node1) + nodeCount]);
		if (getNodeStatus(node2))
			if (strainDestination[getNodeID(node2) + nodeCount]
			    == destination2) {
//                                      puts("Safe");
				counter[1] -= 4;
				counter[0]++;
			} else;
//                              puts("stopped 3");
		else;
//                      puts("stopped 2");
	} else;
//              puts("stopped 1");
}

void displayBreakpoints(Graph * graph, IDnum firstStrain)
{
	IDnum nodeIndex;
	Node *node;
	PassageMarker *strainMarker, *genomeMarker;
	Node **genomeDestination =
	    callocOrExit(2 * nodeCount(graph) + 1, Node *);
	Node **strainDestination =
	    callocOrExit(2 * nodeCount(graph) + 1, Node *);
	IDnum counters[3];

	counters[0] = 0;
	counters[1] = 0;
	counters[2] = 0;

	resetNodeStatus(graph);

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex);

		if (arcCount(node) <= 1
		    && arcCount(getTwinNode(node)) <= 1) {
			continue;
		}

		if (nodeGenomicMultiplicity(node, firstStrain) != 1) {
			continue;
		}

		if (isOnlyGenome(node, firstStrain)) {
			continue;
		}

		for (genomeMarker = getMarker(node); genomeMarker != NULL;
		     genomeMarker = getNextInNode(genomeMarker))
			if (getAbsolutePassMarkerSeqID(genomeMarker) <
			    firstStrain)
				break;

		// Go through all strain passage marker
		for (strainMarker = getMarker(node); strainMarker != NULL;
		     strainMarker = getNextInNode(strainMarker)) {
			displayLocalBreakpoint(strainMarker, firstStrain,
					       genomeMarker,
					       genomeDestination,
					       strainDestination, counters,
					       nodeCount(graph));
			displayLocalBreakpoint(getTwinMarker(strainMarker),
					       firstStrain,
					       getTwinMarker(genomeMarker),
					       genomeDestination,
					       strainDestination, counters,
					       nodeCount(graph));
		}
	}


	printf("%d\t%d\t%d\n", counters[0], counters[1], counters[2]);
	free(strainDestination);
	free(genomeDestination);
}

PassageMarker *genomeMarker(Node * node, IDnum firstStrain)
{
	PassageMarker *marker;

	if (genomeMarkerCount(node, firstStrain) != 1)
		return NULL;

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) < firstStrain)
			return marker;

	return NULL;
}

void exportArcSequence(Arc * arc, FILE * outfile, int WORDLENGTH,
		       TightString ** sequences)
{
	char *str;
	TightString *output =
	    newTightString(getNodeLength(getOrigin(arc)) +
			   getNodeLength(getDestination(arc)));
	appendNodeSequence(getOrigin(arc), output, 0);
	appendNodeSequence(getDestination(arc), output,
			   getNodeLength(getOrigin(arc)));
	str = readTightString(output);
	fprintf(outfile, "> ARC from NODE %d", getNodeID(getOrigin(arc)));
	fprintf(outfile, "%s\n", str);
	destroyTightString(output);
	free(str);
}

// Produce sequences necessary to recreate graph elsewhere...
void projectGraphToFile(Graph * graph, char *filename, int WORDLENGTH,
			TightString ** sequences)
{
	FILE *outfile = fopen(filename, "w");
	IDnum index;
	Node *currentNode;
	Arc *arc;

	if (outfile == NULL) {
		printf("Could not open %s, sorry\n", filename);
		return;
	}

	for (index = 1; index < nodeCount(graph); index++) {
		currentNode = getNodeInGraph(graph, index);
		for (arc = getArc(currentNode); arc != NULL;
		     arc = getNextArc(arc))
			exportArcSequence(arc, outfile, WORDLENGTH,
					  sequences);

		for (arc = getArc(getTwinNode(currentNode)); arc != NULL;
		     arc = getNextArc(arc))
			exportArcSequence(arc, outfile, WORDLENGTH,
					  sequences);
	}

	fclose(outfile);
}

void removeReferenceMarkers(Graph * graph, IDnum firstStrain)
{
	IDnum ID;
	Node *node;
	PassageMarker *marker, *oldMarker;

	for (ID = 1; ID <= nodeCount(graph); ID++) {
		node = getNodeInGraph(graph, ID);
		marker = getMarker(node);
		while (marker != NULL) {
			if (getAbsolutePassMarkerSeqID(marker) <
			    firstStrain) {
				if (!isInitial(marker))
					changeMultiplicity
					    (getArcBetweenNodes
					     (getNode
					      (getPreviousInSequence
					       (marker)), node, graph),
					     -1);
				if (!isTerminal(marker))
					changeMultiplicity
					    (getArcBetweenNodes
					     (node,
					      getNode(getNextInSequence
						      (marker)), graph),
					     -1);
				oldMarker = marker;
				marker = getNextInNode(marker);
				destroyPassageMarker(oldMarker);
			} else
				marker = getNextInNode(marker);
		}

		if (getMarker(node) == NULL)
			destroyNode(node, graph);
	}

	concatenateGraph(graph);
}

void exportLongNodeSequences(char *filename, Graph * graph,
			     Coordinate minLength)
{
	FILE *outfile = fopen(filename, "w");
	IDnum nodeIndex;
	TightString *tString;
	Coordinate position;
	char nucleotide;
	Node *node;
	int WORDLENGTH = getWordLength(graph);
	GapMarker *gap;
	//double sensitivity, specificity;

	if (outfile == NULL) {
		printf("Could not write into %s, sorry\n", filename);
		return;
	} else {
		printf("Writing contigs into %s...\n", filename);
	}

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex);

		if (node == NULL || getNodeLength(node) < minLength)
			continue;

		tString = expandNode(node, WORDLENGTH);
		fprintf(outfile, ">NODE_%ld_length_%lld_cov_%f\n",
			(long) nodeIndex, (long long) getNodeLength(node),
			(getVirtualCoverage(node, 0)
			 + getVirtualCoverage(node, 1)
			 + readCoverage(node)) /
			(float) getNodeLength(node));

		gap = getGap(node, graph);
		for (position = 0; position < WORDLENGTH; position++) {
			if (gap && position >= getGapFinish(gap))
				gap = getNextGap(gap);

			if (gap == NULL || position < getGapStart(gap)) {
				nucleotide =
				    getNucleotideChar(position, tString);
				fprintf(outfile, "%c", nucleotide);
			} else
				fprintf(outfile, "N");
		}

		gap = getGap(node, graph);
		for (; position < getLength(tString); position++) {
			if (position % 60 == 0)
				fprintf(outfile, "\n");

			if (gap
			    && position - WORDLENGTH + 1 >=
			    getGapFinish(gap))
				gap = getNextGap(gap);

			if (gap == NULL
			    || position - WORDLENGTH + 1 <
			    getGapStart(gap)) {
				nucleotide =
				    getNucleotideChar(position, tString);
				fprintf(outfile, "%c", nucleotide);
			} else
				fprintf(outfile, "N");
		}
		fprintf(outfile, "\n");

		destroyTightString(tString);
	}

	fclose(outfile);
}

/*
void exportMediumNodeSequences(char* filename, Graph * graph, Coordinate minLength)
{
	IDnum dummy;
	ReadSet *readSet = readFastAFile(sequenceFile);
	char **reads = readSet->sequences;
	TightString **sequences =
	    newTightStringArrayFromStringArray(reads, dummy);
	FILE *outfile = fopen(filename, "w");
	char str[100];
	IDnum nodeIndex;
	TightString *tString;
	Coordinate start, finish;
	Node *node;
	double sensitivity, specificity;

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex);
		if (getNodeLength(node) < minLength
		    || getNodeLength(node) >= maxLength)
			continue;

		tString = expandNode(node);
		compareSequences(tString, sequences[0], &sensitivity,
				 &specificity, WORDLENGTH);

		fprintf(outfile,
			"> MEDIUM NODE %d, Sensitivity = %f, Specificity = %f\n",
			nodeIndex, sensitivity, specificity);
		printf
		    ("> MEDIUM NODE %d, Sensitivity = %f, Specificity = %f\n",
		     nodeIndex, sensitivity, specificity);

		start = 0;
		while (start <= getLength(tString)) {
			finish = start + 60;
			readTightStringFragment(tString, start,
						finish, str);
			fprintf(outfile, "%s\n", str);
			start = finish;
		}

		destroyTightString(tString);
	}
}
*/

Coordinate maxLength(Graph * graph)
{
	IDnum index;
	Node *node;
	Coordinate max = 0;

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node != NULL && getNodeLength(node) > max)
			max = getNodeLength(node);
	}

	return max;
}

Coordinate n50(Graph * graph)
{
	FibHeap *heap = newFibHeap();
	IDnum index;
	Coordinate totalLength = 0;
	Coordinate sumLength = 0;
	Node *node;

	if (nodeCount(graph) == 0) {
		puts("EMPTY GRAPH");
		return 0;
	}

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL)
			continue;
		insertNodeIntoHeap(heap, getNodeLength(node), node);
		totalLength += getNodeLength(node);
	}
	totalLength /= 2;

	node = removeNextNodeFromHeap(heap);
	while (node != NULL) {
		sumLength += getNodeLength(node);
		if (sumLength >= totalLength)
			break;
		node = removeNextNodeFromHeap(heap);
	}

	destroyHeap(heap);
	return getNodeLength(node);
}

static Coordinate getTotalCoverage(Node * node)
{
	Category cat;
	Coordinate coverage = 0;

	for (cat = 0; cat < CATEGORIES; cat++)
		coverage += getVirtualCoverage(node, cat);

	return coverage;
}

int compareNodeCovs(const void * A, const void * B) {
	Node * nodeA = *((Node **) A);
	Node * nodeB = *((Node **) B);
	double covA;
	double covB;
	
	if (getNodeLength(nodeA) == 0)
		nodeA = NULL;

	if (getNodeLength(nodeB) == 0)
		nodeB = NULL;

	// Null nodes considered to have infinite coverage
	if (nodeA == NULL && nodeB == NULL)
		return 0;
	if (nodeA == NULL)
		return 1;
	if (nodeB == NULL)
		return -1;

	// Deal with real coverage numbers:
	covA = getTotalCoverage(nodeA) / (double) getNodeLength(nodeA);	
	covB = getTotalCoverage(nodeB) / (double) getNodeLength(nodeB);	

	if (covA > covB)
		return 1;
	if (covA == covB)
		return 0;
	return -1;
}

double estimated_cov(Graph * graph, char * directory)
{
	Node ** nodeArray = callocOrExit(nodeCount(graph), Node*); 
	IDnum index;
	Coordinate halfTotalLength = 0;
	Coordinate sumLength = 0;
	Node *node;
	char *logFilename =
	    mallocOrExit(strlen(directory) + 100, char);
	char *statsLine = 
	    mallocOrExit(5000, char);
	FILE *logFile;

	strcpy(logFilename, directory);
	strcat(logFilename, "/Log");
	logFile = fopen(logFilename, "a");

	if (logFile == NULL)
		exitErrorf(EXIT_FAILURE, true, "Could not write to %s",
		       logFilename);

	puts("Measuring median coverage depth...");

	if (nodeCount(graph) == 0) {
		puts("EMPTY GRAPH");
		return 0;
	}

	// Write nodes into array and compute total assembly length
	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		nodeArray[index - 1] = node;
		if (node == NULL)
			continue;
		halfTotalLength += getNodeLength(node);
	}
	halfTotalLength /= 2;

	// Sort nodes
	qsort(nodeArray, nodeCount(graph), sizeof(Node *), compareNodeCovs);

	// Compute the length weighted median node coverage
	for (index = 0; index < nodeCount(graph); index++) {
		node = nodeArray[index];
		sumLength += getNodeLength(node);
		if (sumLength >= halfTotalLength) {
			printf("Median coverage depth = %f\n", getTotalCoverage(node) / (double) getNodeLength(node));
			fprintf(logFile, "Median coverage depth = %f\n", getTotalCoverage(node) / (double) getNodeLength(node));
			free(nodeArray);
			fclose(logFile);
			free(logFilename);
			free(statsLine);
			return getTotalCoverage(node) / (double) getNodeLength(node);
		}
	}

	// In case something went wrong...
	free(nodeArray);
	fclose(logFile);
	free(logFilename);
	free(statsLine);

	return -1;
}	

static void destroyMixedNode(Node * node)
{
	PassageMarker *marker = getMarker(node);
	PassageMarker *current;

	while (marker != NULL) {
		while (!isInitial(marker))
			marker = getPreviousInSequence(marker);

		while (marker != NULL) {
			current = marker;
			marker = getNextInSequence(marker);
			destroyPassageMarker(current);
		}

		marker = getMarker(node);
	}
}

void destroyMixedReads(Graph * graph, IDnum minCoverage)
{
	IDnum index;
	Node *node;

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL)
			continue;

		if (markerCount(node) < minCoverage)
			destroyMixedNode(node);
	}

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL)
			continue;

		if (getMarker(node) == NULL)
			destroyNode(node, graph);
	}

	concatenateGraph(graph);
}

static boolean isConnectedRead(PassageMarker * marker)
{
	PassageMarker *current;

	for (current = marker; getNodeStatus(getNode(current));
	     current = getNextInSequence(current))
		if (isTerminal(current))
			return false;

	for (current = getTwinMarker(marker);
	     getNodeStatus(getNode(current));
	     current = getNextInSequence(current))
		if (isTerminal(current))
			return false;

	return true;
}

static void destroyWholeRead(PassageMarker * marker)
{
	PassageMarker *current = marker;
	PassageMarker *next;

	while (!isInitial(current))
		current = getPreviousInSequence(current);

	for (; current != NULL; current = next) {
		next = getNextInSequence(current);
		destroyPassageMarker(current);
	}
}

static void cleanUpNode(Node * node, Graph * graph)
{
	Category cat;
	Node *twin = getTwinNode(node);
	PassageMarker *marker, *twinMarker;

	for (cat = 0; cat < CATEGORIES; cat++)
		setVirtualCoverage(node, cat, 0);

	while (getArc(node) != NULL)
		destroyArc(getArc(node), graph);
	while (getArc(twin) != NULL)
		destroyArc(getArc(twin), graph);

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker)) {
		twinMarker = getTwinMarker(marker);

		if (getNode(getNextInSequence(marker)) != twin
		    || getPassageMarkerSequenceID(marker) > 0)
			createArc(node,
				  getNode(getNextInSequence(marker)),
				  graph);
		if (getNode(getNextInSequence(twinMarker)) != node
		    || getPassageMarkerSequenceID(twinMarker) > 0)
			createArc(twin,
				  getNode(getNextInSequence(twinMarker)),
				  graph);
	}
}

void destroySinglePoolNodes(Graph * graph)
{
	IDnum index;
	Node *node;
	PassageMarker *marker, *next;

	puts("Destroying single pool nodes");
	resetNodeStatus(graph);

	// Remove empty, single pool nodes, mark other single pool nodes
	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL)
			continue;

		if (getVirtualCoverage(node, 0) != 0
		    && getVirtualCoverage(node, 1) != 0)
			continue;

		if (getMarker(node) == NULL)
			destroyNode(node, graph);
		else
			setNodeStatus(node, true);
	}

	// Remove disconnected reads
	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL || !getNodeStatus(node))
			continue;

		for (marker = getMarker(node); marker != NULL;
		     marker = next) {
			if (isConnectedRead(marker))
				next = getNextInNode(marker);
			else {
				destroyWholeRead(marker);
				next = getMarker(node);
			}
		}
	}

	// Remove empty, single pool nodes, review coverage of the others
	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL || !getNodeStatus(node))
			continue;

		if (getMarker(node) == NULL)
			destroyNode(node, graph);
		else
			cleanUpNode(node, graph);
	}

	puts("Done");

	concatenateGraph(graph);
}

void destroyShortTips(Graph * graph)
{
	IDnum index;
	Node *node;
	boolean modified = true;

	puts("Removing short tips");

	while (modified) {
		modified = false;
		for (index = 1; index <= nodeCount(graph); index++) {
			node = getNodeInGraph(graph, index);
			if (node == NULL)
				continue;

			if (getArc(node) == NULL
			    || getArc(getTwinNode(node)) == NULL) {
				if (getNodeLength(node) < 500) {
					modified = true;
					destroyNode(node, graph);
				}
			}
		}
	}

	puts("Done");

	concatenateGraph(graph);
}

static Coordinate connectDomain(Node * node)
{
	Coordinate result = getNodeLength(node);
	Arc *arc;

	if (getNodeStatus(node))
		return 0;
	setNodeStatus(node, true);

	for (arc = getArc(node); arc != NULL; arc = getNextArc(arc))
		result += connectDomain(getDestination(arc));
	for (arc = getArc(getTwinNode(node)); arc != NULL;
	     arc = getNextArc(arc))
		result += connectDomain(getDestination(arc));

	return result;

}

static void destroyConnectedDomain(Node * node, Graph * graph)
{
	Arc *arc;

	if (getNodeStatus(node))
		return;
	setNodeStatus(node, true);

	for (arc = getArc(node); arc != NULL; arc = getNextArc(arc))
		destroyConnectedDomain(getDestination(arc), graph);
	for (arc = getArc(getTwinNode(node)); arc != NULL;
	     arc = getNextArc(arc))
		destroyConnectedDomain(getDestination(arc), graph);

	destroyNode(node, graph);

}

void destroyDisconnectedElements(Graph * graph)
{
	Node *node;
	IDnum index;
	Coordinate domainSize;
	FibHeap *heap = newFibHeap();
	Coordinate *domainSizes =
	    callocOrExit(1 + nodeCount(graph), Coordinate);

	resetNodeStatus(graph);

	puts("Destroying disconnected domains");

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL || getNodeStatus(node))
			continue;
		domainSize = connectDomain(node);
		printf("CONNECT\t%lld\n", (long long) domainSize);
		insertNodeIntoHeap(heap, domainSize, node);
		domainSizes[index] = domainSize;
	}

	resetNodeStatus(graph);

	while (true) {
		node = removeNextNodeFromHeap(heap);
		if (node == NULL || domainSizes[getNodeID(node)] > 1200)
			break;

		destroyConnectedDomain(node, graph);
	}


	destroyHeap(heap);
	free(domainSizes);
	puts("Done");

	concatenateGraph(graph);
}

static Coordinate connectDomainNodeCount(Node * node)
{
	Coordinate result = 1;
	Arc *arc;

	if (getNodeStatus(node))
		return 0;
	setNodeStatus(node, true);

	for (arc = getArc(node); arc != NULL; arc = getNextArc(arc))
		result += connectDomain(getDestination(arc));
	for (arc = getArc(getTwinNode(node)); arc != NULL;
	     arc = getNextArc(arc))
		result += connectDomain(getDestination(arc));

	return result;

}

void measureTangleSizes(Graph * graph, Coordinate maxLength)
{
	Node *node;
	IDnum index;
	Coordinate domainSize;

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL)
			continue;
		if (getNodeLength(node) >= maxLength)
			destroyNode(node, graph);
	}

	renumberNodes(graph);
	resetNodeStatus(graph);

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL || getNodeStatus(node))
			continue;
		domainSize = connectDomainNodeCount(node);
		printf("CONNECT\t%lld\n", (long long) domainSize);
	}

	puts("Done");
}

void destroyEmptyNodes(Graph * graph)
{
	IDnum index;
	Node *node;

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL)
			continue;

		if (getMarker(node) == NULL)
			destroyNode(node, graph);
	}

	concatenateGraph(graph);
}

static Coordinate getReadLength(PassageMarker * marker)
{
	PassageMarker *current;
	Coordinate totalLength = 0;

	for (current = marker; current != NULL;
	     current = getNextInSequence(current))
		totalLength += getPassageMarkerLength(current);

	return totalLength;
}

static void destroyRead(PassageMarker * marker)
{
	PassageMarker *current, *next;

	for (current = marker; current != NULL; current = next) {
		next = getNextInSequence(current);
		destroyPassageMarker(current);
	}
}

void removeShortReads(Graph * graph)
{
	IDnum index;
	Node *node;
	PassageMarker *marker, *next;

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL)
			continue;


		for (marker = getMarker(node); marker != NULL;
		     marker = next) {
			if (getPassageMarkerSequenceID(marker) > 0
			    && isInitial(marker)
			    && getReadLength(marker) < 400) {
				destroyRead(marker);
				next = getMarker(node);
			} else if (getPassageMarkerSequenceID(marker) < 0
				   && isTerminal(marker)
				   && getReadLength(getTwinMarker(marker))
				   < 400) {
				destroyRead(getTwinMarker(marker));
				next = getMarker(node);
			} else
				next = getNextInNode(marker);

		}
	}
}

Coordinate totalGraphLength(Graph * graph)
{
	IDnum index;
	Coordinate totalLength = 0;

	for (index = 1; index <= nodeCount(graph); index++)
		totalLength += getNodeLength(getNodeInGraph(graph, index));

	return totalLength;
}

void destroySinglePoolNodesStrict(Graph * graph)
{
	IDnum index;
	Node *node;

	puts("Destroying single pool nodes");
	resetNodeStatus(graph);

	// Remove empty, single pool nodes, mark other single pool nodes
	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL)
			continue;

		if (getVirtualCoverage(node, 0) != 0
		    && getVirtualCoverage(node, 1) != 0)
			continue;

		destroyNode(node, graph);
	}

	puts("Done");

	concatenateGraph(graph);
}

void contigStats(Node ** contigs, IDnum readCount)
{
	FibHeap *heap = newFibHeap();
	IDnum index;
	Coordinate totalLength = 0;
	Coordinate sumLength = 0;
	Node *node;
	Coordinate halfLength;

	for (index = 0; index <= readCount; index++) {
		if (contigs[index] != NULL) {
			node = contigs[index];
			printf("CONTIG %lld\n", (long long) getNodeLength(node));
			insertNodeIntoHeap(heap, getNodeLength(node),
					   node);
			totalLength += getNodeLength(node);
		}
	}
	halfLength = totalLength / 2;

	node = removeNextNodeFromHeap(heap);
	while (node != NULL) {
		sumLength += getNodeLength(node);
		if (sumLength >= halfLength)
			break;
		node = removeNextNodeFromHeap(heap);
	}

	destroyHeap(heap);
	printf("N50 %lld Total %lld\n", (long long) getNodeLength(node), (long long) totalLength);
}

void exportContigs(Node ** contigs, ReadSet * reads, char *filename,
		   int WORDLENGTH, int pairedReadsCount)
{
	TightString **sequences =
	    mallocOrExit(reads->readCount, TightString *);
	IDnum i;

	for (i = 0; i < pairedReadsCount; i++) {
		if (contigs[i] == NULL)
			sequences[i] = NULL;
		else
			sequences[i] = expandNode(contigs[i], WORDLENGTH);
	}

	exportSequenceArray(filename, sequences, reads->readCount);
}

boolean *removeLowCoverageNodesAndDenounceDubiousReads(Graph * graph,
						       double minCov)
{
	IDnum index;
	Node *node;
	boolean denounceReads = readStartsAreActivated(graph);
	boolean *res = NULL; 
	ShortReadMarker *nodeArray, *shortMarker;
	PassageMarker *marker;
	IDnum maxIndex;
	IDnum readID;
	IDnum index2;
	// Original
	double nodeDensity = 0.0;
	int countNodeTotal = 0, countNodeUnderMinCov = 0, countNodeOnLongRead = 0; 
	int countNodeUMCbutOnLongRead = 0, countNodeSupportedByLongRead = 0;
	int countNodeEscapedByLongRead = 0;
	boolean escapeByLongRead = false;
	// Original

	printf("Removing contigs with coverage < %f...\n", minCov);
		
	if (denounceReads)
		res = callocOrExit(sequenceCount(graph), boolean);

	// Original
	countNodeTotal = nodeCount(graph);
	// Original

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (getNodeLength(node) == 0)
			continue;

		// Original 
		if ((marker = getMarker(node))){
			countNodeOnLongRead++;
		}
		// Original
		
		if ((double)getTotalCoverage(node) / getNodeLength(node) < minCov) {
		        // Original
		        countNodeUnderMinCov++;
			if ((marker = getMarker(node))) {
				nodeDensity = (double)getTotalCoverage(node) / getNodeLength(node);
				// printf("[Under]\tnodeID : %d\tnodeLen : %d\tnodeCov : %f\n",
				//       getNodeID(node), getNodeLength(node), nodeDensity);
				countNodeUMCbutOnLongRead++;

				if (markerCount(node) >= 1) {
					countNodeSupportedByLongRead++;
				}
				if (markerCount(node) >= 1 && !(nodeDensity == 0.0)) {
					escapeByLongRead = true;
					countNodeEscapedByLongRead++;
				} else {
					escapeByLongRead = false;
				}
			}
			// Original

			
			// Original
			//if (! escapeByLongRead) {
			// Original

			if (denounceReads) {
				nodeArray = getNodeReads(node, graph);
				maxIndex = getNodeReadCount(node, graph);
				for (index2 = 0; index2 < maxIndex; index2++) {
					shortMarker =
					    getShortReadMarkerAtIndex(nodeArray,
								      index2);
					readID = getShortReadMarkerID(shortMarker);
					//printf("Dubious %d\n", readID);
					if (readID > 0)
						res[readID - 1] = true;
					else
						res[-readID - 1] = true;
				}

				nodeArray = getNodeReads(getTwinNode(node), graph);
				maxIndex =
				    getNodeReadCount(getTwinNode(node), graph);
				for (index2 = 0; index2 < maxIndex; index2++) {
					shortMarker =
					    getShortReadMarkerAtIndex(nodeArray,
								      index2);
					readID = getShortReadMarkerID(shortMarker);
					//printf("Dubious %d\n", readID);
					if (readID > 0)
						res[readID - 1] = true;
					else
						res[-readID - 1] = true;
				}
			}
			
			// Original
			//}
			// Original

			while ((marker = getMarker(node))) {
			        if (!isInitial(marker)
				    && !isTerminal(marker))
					disconnectNextPassageMarker
					    (getPreviousInSequence(marker),
					     graph);
				destroyPassageMarker(marker);
			}
			destroyNode(node, graph);

		}
		// Original
		// else {
		//	nodeDensity = (double)getTotalCoverage(node) / getNodeLength(node);
		//	printf("[Over]\tnodeID : %d\tnodeLen : %d\tnodeCov : %f\n",
		//	       getNodeID(node), getNodeLength(node), nodeDensity);
		//}
		// Original
	}
	
	/*
	// Original
	printf("No. of Nodes : %d\n", countNodeTotal);
	printf("No. of Nodes Under MinCov : %d\n", countNodeUnderMinCov);
	printf("No. of Nodes On Long Reads : %d\n", countNodeOnLongRead);
	printf("No. of Nodes Under MinCov but On Long Reads : %d\n", countNodeOnLongRead);
	printf("No. of Nodes Under MinCov but Supported by Long Reads : %d\n",
	       countNodeSupportedByLongRead);
	printf("No. of Nodes 0 < Cov < minCov but Escaped (Supported) by Long Reads : %d\n",
	       countNodeEscapedByLongRead);
	// Original
	*/

	concatenateGraph(graph);
	return res;
}

void removeLowCoverageNodes(Graph * graph, double minCov)
{
	IDnum index;
	Node *node;
	PassageMarker *marker;

	if (minCov < 0)
		return;

	printf("Applying a coverage cutoff of %f...\n", minCov);

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (getNodeLength(node) > 0
		    && getTotalCoverage(node) / getNodeLength(node) <
		    minCov) {
			while ((marker = getMarker(node))) {
				if (!isInitial(marker)
				    && !isTerminal(marker))
					disconnectNextPassageMarker
					    (getPreviousInSequence(marker),
					     graph);
				destroyPassageMarker(marker);
			}
			destroyNode(node, graph);
		}
	}

	concatenateGraph(graph);
}

void removeHighCoverageNodes(Graph * graph, double maxCov)
{
	IDnum index;
	Node *node;
	PassageMarker *marker;

	if (maxCov < 0)
		return;

	printf("Applying an upper coverage cutoff of %f...\n", maxCov);

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (getNodeLength(node) > 0
		    && getTotalCoverage(node) / getNodeLength(node) >
		    maxCov) {
			while ((marker = getMarker(node))) {
				if (!isInitial(marker)
				    && !isTerminal(marker))
					disconnectNextPassageMarker
					    (getPreviousInSequence(marker),
					     graph);
				destroyPassageMarker(marker);
			}
			destroyNode(node, graph);
		}
	}

	concatenateGraph(graph);
}

void removeMissingStrain(Graph * graph, Category cat)
{
	IDnum ID;
	Node *node;

	for (ID = 1; ID <= nodeCount(graph); ID++) {
		node = getNodeInGraph(graph, ID);

		if (node == NULL)
			continue;

		if (getVirtualCoverage(node, cat) == 0)
			destroyNode(node, graph);
	}

	concatenateGraph(graph);
}

static void exportAMOSLib(FILE * outfile, Graph * graph, Category cat)
{
	Coordinate distance = getInsertLength(graph, cat * 2);
	double variance = getInsertLength_var(graph, cat * 2);

	if (distance == -1)
		return;

	fprintf(outfile, "{LIB\n");
	fprintf(outfile, "iid:%d\n", (int) (cat + 1));
	fprintf(outfile, "{DST\n");
	fprintf(outfile, "mea:%lld\n", (long long) distance);
	fprintf(outfile, "std:%lld\n", (long long) sqrt(variance));
	fprintf(outfile, "}\n");
	fprintf(outfile, "}\n");
}

static void exportAMOSMarker(FILE * outfile, PassageMarker * marker,
			     Coordinate nodeLength, Coordinate start,
			     Coordinate finish, int wordShift)
{
	Coordinate sequenceStart, sequenceFinish;

	if (getStartOffset(marker) >= finish
	    || getFinishOffset(marker) > nodeLength - start)
		return;

	sequenceStart = getPassageMarkerStart(marker);
	if (start > getStartOffset(marker)) {
		if (getPassageMarkerSequenceID(marker) > 0)
			sequenceStart += start - getStartOffset(marker);
		else
			sequenceStart -= start - getStartOffset(marker);
	}

	sequenceFinish = getPassageMarkerFinish(marker);
	if (nodeLength - finish > getFinishOffset(marker)) {
		if (getPassageMarkerSequenceID(marker) > 0)
			sequenceFinish -=
			    nodeLength - finish - getFinishOffset(marker);
		else
			sequenceFinish +=
			    nodeLength - finish - getFinishOffset(marker);
	}

	if (getPassageMarkerSequenceID(marker) > 0)
		sequenceFinish += wordShift;
	else
		sequenceStart += wordShift;

	fprintf(outfile, "{TLE\n");
	fprintf(outfile, "src:%d\n", getAbsolutePassMarkerSeqID(marker));
	if (getStartOffset(marker) > start)
		fprintf(outfile, "off:%lld\n",
			(long long) (getStartOffset(marker) - start));
	else
		fprintf(outfile, "off:0\n");
	fprintf(outfile, "clr:%lld,%lld\n", (long long) sequenceStart, (long long) sequenceFinish);
	fprintf(outfile, "}\n");
}

static void exportAMOSShortMarker(FILE * outfile, ShortReadMarker * marker,
				  ReadSet * reads, Coordinate start,
				  Coordinate finish)
{
	Coordinate offset =
	    getShortReadMarkerPosition(marker) -
	    getShortReadMarkerOffset(marker);
	TightString *sequence =
	    reads->tSequences[getShortReadMarkerID(marker) - 1];

	if (getShortReadMarkerPosition(marker) == -1)
		return;

	if (offset >= finish || offset + getLength(sequence) < start)
		return;

	fprintf(outfile, "{TLE\n");
	fprintf(outfile, "src:%d\n", getShortReadMarkerID(marker));
	fprintf(outfile, "off:%lld\n", (long long) (offset - start));
	fprintf(outfile, "clr:0,%lld\n", (long long) getLength(sequence));
	fprintf(outfile, "}\n");
}

static void exportAMOSReverseShortMarker(FILE * outfile,
					 ShortReadMarker * marker,
					 Coordinate nodeLength,
					 int wordShift, ReadSet * reads,
					 Coordinate start,
					 Coordinate finish)
{
	TightString *sequence =
	    reads->tSequences[getShortReadMarkerID(marker) - 1];

	Coordinate offset =
	    nodeLength - getShortReadMarkerPosition(marker) +
	    getShortReadMarkerOffset(marker) - getLength(sequence) +
	    wordShift;

	if (getShortReadMarkerPosition(marker) == -1)
		return;

	if (offset >= finish || offset + getLength(sequence) < start)
		return;

	fprintf(outfile, "{TLE\n");
	fprintf(outfile, "src:%d\n", getShortReadMarkerID(marker));
	fprintf(outfile, "off:%lld\n", (long long) (offset - start));
	fprintf(outfile, "clr:%lld,0\n", (long long) getLength(sequence));
	fprintf(outfile, "}\n");
}

static void exportAMOSContig(FILE * outfile, ReadSet * reads, Node * node,
			     Graph * graph, Coordinate contigStart,
			     Coordinate contigFinish, IDnum iid,
			     IDnum internalIndex)
{
	Coordinate start;
	char str[100];
	PassageMarker *marker;
	ShortReadMarker *shortMarkerArray, *shortMarker;
	Coordinate index, maxIndex;
	int wordShift = getWordLength(graph) - 1;
	char *string = expandNodeFragment(node, contigStart, contigFinish,
					  getWordLength(graph));
	Coordinate length = contigFinish - contigStart + wordShift;

	fprintf(outfile, "{CTG\n");
	fprintf(outfile, "iid:%d\n", iid);
	fprintf(outfile, "eid:%d-%d\n", getNodeID(node), internalIndex);

	fprintf(outfile, "seq:\n");
	for (start = 0; start <= length; start += 60) {
		strncpy(str, &(string[start]), 60);
		str[60] = '\0';
		fprintf(outfile, "%s\n", str);
	}
	fprintf(outfile, ".\n");

	fprintf(outfile, "qlt:\n");
	for (start = 0; start <= length; start += 60) {
		strncpy(str, &(string[start]), 60);
		str[60] = '\0';
		fprintf(outfile, "%s\n", str);
	}
	fprintf(outfile, ".\n");

	free(string);

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker))
		exportAMOSMarker(outfile, marker, getNodeLength(node),
				 contigStart, contigFinish, wordShift);

	if (readStartsAreActivated(graph)) {
		shortMarkerArray = getNodeReads(node, graph);
		maxIndex = getNodeReadCount(node, graph);
		for (index = 0; index < maxIndex; index++) {
			shortMarker =
			    getShortReadMarkerAtIndex(shortMarkerArray,
						      index);
			exportAMOSShortMarker(outfile, shortMarker, reads,
					      contigStart, contigFinish);
		}

		shortMarkerArray = getNodeReads(getTwinNode(node), graph);
		maxIndex = getNodeReadCount(getTwinNode(node), graph);
		for (index = 0; index < maxIndex; index++) {
			shortMarker =
			    getShortReadMarkerAtIndex(shortMarkerArray,
						      index);
			exportAMOSReverseShortMarker(outfile, shortMarker,
						     getNodeLength(node),
						     wordShift, reads,
						     contigStart,
						     contigFinish);
		}
	}

	fprintf(outfile, "}\n");
}

static void exportAMOSNode(FILE * outfile, ReadSet * reads, Node * node,
			   Graph * graph)
{
	Coordinate start = 0;
	Coordinate finish;
	GapMarker *gap;
	IDnum smallIndex = 0;
	static IDnum iid = 1;
	IDnum contigIndex = iid;
	int wordShift = getWordLength(graph) - 1;

	for (gap = getGap(node, graph); gap; gap = getNextGap(gap)) {
		finish = getGapStart(gap);
		exportAMOSContig(outfile, reads, node, graph, start,
				 finish, iid++, smallIndex++);
		start = getGapFinish(gap);
	}

	finish = getNodeLength(node);
	exportAMOSContig(outfile, reads, node, graph, start, finish, iid++,
			 smallIndex);

	if (!getGap(node, graph))
		return;

	start = 0;

	fprintf(outfile, "{SCF\n");
	fprintf(outfile, "eid:%d\n", getNodeID(node));
	for (gap = getGap(node, graph); gap; gap = getNextGap(gap)) {
		finish = getGapStart(gap);
		fprintf(outfile, "{TLE\n");
		fprintf(outfile, "off:%lld\n", (long long) start);
		fprintf(outfile, "clr:0,%lld\n",
			(long long) (finish - start + (long long) wordShift));
		fprintf(outfile, "src:%d\n", contigIndex++);
		fprintf(outfile, "}\n");
		start = getGapFinish(gap);
	}
	finish = getNodeLength(node);
	fprintf(outfile, "{TLE\n");
	fprintf(outfile, "off:%lld\n", (long long) start);
	fprintf(outfile, "clr:0,%lld\n", (long long) (finish - start));
	fprintf(outfile, "src:%d\n", contigIndex++);
	fprintf(outfile, "}\n");

	fprintf(outfile, "}\n");
}

static void exportAMOSRead(FILE * outfile, TightString * tString,
			   IDnum index, IDnum frg_index)
{
	Coordinate start, finish;
	char str[100];

	fprintf(outfile, "{RED\n");
	fprintf(outfile, "iid:%d\n", index);
	fprintf(outfile, "eid:%d\n", index);
	if (frg_index > 0)
		fprintf(outfile, "frg:%d\n", frg_index);

	fprintf(outfile, "seq:\n");
	start = 0;
	while (start <= getLength(tString)) {
		finish = start + 60;
		readTightStringFragment(tString, start, finish, str);
		fprintf(outfile, "%s\n", str);
		start = finish;
	}
	fprintf(outfile, ".\n");

	fprintf(outfile, "qlt:\n");
	start = 0;
	while (start <= getLength(tString)) {
		finish = start + 60;
		readTightStringFragment(tString, start, finish, str);
		fprintf(outfile, "%s\n", str);
		start = finish;
	}
	fprintf(outfile, ".\n");

	fprintf(outfile, "}\n");
}

void exportAMOSContigs(char *filename, Graph * graph,
		       Coordinate cutoff_length, ReadSet * reads)
{
	IDnum index;
	Category cat;
	Node *node;
	FILE *outfile;

	printf("Writing into AMOS file %s...\n", filename);
	outfile = fopen(filename, "w");

	if (outfile == NULL)
		exitErrorf(EXIT_FAILURE, true, "Could not write to AMOS file %s",
		       filename);

	for (cat = 0; cat <= CATEGORIES; cat++)
		exportAMOSLib(outfile, graph, cat);

	for (index = 1; index <= reads->readCount; index++) {
		if (reads->categories[index - 1] % 2 != 0 &&
		    getInsertLength(graph,
				    reads->categories[index - 1]) >= 0) {
			fprintf(outfile, "{FRG\n");
			fprintf(outfile, "lib:%d\n",
				(int) ((reads->categories[index - 1] / 2) + 1));
			fprintf(outfile, "rds:%d,%d\n", index,
				index + 1);
			fprintf(outfile, "eid:%d\n", index);
			fprintf(outfile, "iid:%d\n", index);
			fprintf(outfile, "typ:I\n");
			fprintf(outfile, "}\n");
			index++;
		}
	}

	for (index = 1; index <= reads->readCount; index++) {
		if (reads->categories[index - 1] % 2 != 0 &&
		    getInsertLength(graph,
				    reads->categories[index - 1]) >= 0) {
			exportAMOSRead(outfile,
				       reads->tSequences[index - 1], index,
				       index);
			index++;
			exportAMOSRead(outfile,
				       reads->tSequences[index - 1], index,
				       index - 1);
		} else {
			exportAMOSRead(outfile,
				       reads->tSequences[index - 1], index,
				       -1);
		}
	}

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (node == NULL)
			continue;

		if (getNodeLength(node) >= cutoff_length)
			exportAMOSNode(outfile, reads, node, graph);
	}

	fclose(outfile);

}

boolean isNatural(Graph * graph)
{
	Node *node;
	IDnum index;

	for (index = 1; index < nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL)
			continue;

		if (getNodeLength(node) == 0)
			return false;

		if (simpleArcCount(node) > 4)
			return false;

		if (simpleArcCount(getTwinNode(node)) > 4)
			return false;
	}

	return true;
}

static Node *followShortPath(Arc * arc)
{
	Node *node = getDestination(arc);
	if (simpleArcCount(node) != 1)
		return NULL;
	node = getDestination(getArc(node));
	return getTwinNode(node);

}

static void checkNodeForHallidayJunction(Node * node, Graph * graph)
{
	Node *nodeA, *nodeB, *nodeC, *nodeD;
	Arc *arc1, *arc2;

	if (simpleArcCount(node) != 2)
		return;

	arc1 = getArc(node);
	arc2 = getNextArc(arc1);

	nodeA = followShortPath(arc1);
	if (nodeA == NULL || simpleArcCount(nodeA) != 2
	    || !isUniqueBasic(nodeA))
		return;

	nodeB = followShortPath(arc2);
	if (nodeB == NULL || simpleArcCount(nodeB) != 2
	    || !isUniqueBasic(nodeB))
		return;

	if (nodeA == nodeB) {
		return;
	}

	arc1 = getArc(nodeA);
	arc2 = getNextArc(arc1);
	nodeC = followShortPath(arc1);
	if (nodeC == NULL)
		return;
	if (nodeC == node) {
		nodeC = followShortPath(arc2);
		if (nodeC == NULL || nodeC == node
		    || simpleArcCount(nodeC) != 2
		    || !isUniqueBasic(nodeC)) {
			printf("NO %d %d %d %d\n", getNodeID(node),
			       getNodeID(nodeA), getNodeID(nodeB),
			       getNodeID(nodeC));
			return;
		}
	} else {
		if (simpleArcCount(nodeC) != 2 || !isUniqueBasic(nodeC)) {
			puts("2");
			return;
		}
		nodeD = followShortPath(arc2);
		if (nodeD != node) {
			puts("3");
			return;
		}
	}

	puts("A");

	arc1 = getArc(nodeB);
	arc2 = getNextArc(arc1);
	nodeD = followShortPath(arc1);
	if (nodeD != nodeC && nodeD != node)
		return;
	nodeD = followShortPath(arc2);
	if (nodeD != nodeC && nodeD != node)
		return;

	arc1 = getArc(nodeB);
	arc2 = getNextArc(arc1);
	nodeD = followShortPath(arc1);
	if (nodeD != nodeC && nodeD != node)
		return;
	nodeD = followShortPath(arc2);
	if (nodeD != nodeC && nodeD != node)
		return;

	printf("JOHNNY HALLIDAY JUNCTION %d %d %d %d\n",
	       getNodeID(node), getNodeID(nodeC), getNodeID(nodeA),
	       getNodeID(nodeB));
}

void searchForHallidayJunction(Graph * graph)
{
	IDnum index;
	Node *node;

	setBaseCoverage(8);

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (isUniqueBasic(node)) {
			checkNodeForHallidayJunction(node, graph);
			checkNodeForHallidayJunction(getTwinNode(node),
						     graph);
		}
	}
}

Coordinate totalAssemblyLength(Graph * graph)
{
	IDnum index;
	Node *node;
	Coordinate total = 0;

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node)
			total += getNodeLength(node);
	}

	return total;
}

IDnum usedReads(Graph * graph, Coordinate minContigLength) 
{
	IDnum res = 0;
	boolean * used = callocOrExit(sequenceCount(graph) + 1, boolean);
	IDnum nodeID, readID;
	Node * node;
	PassageMarker * marker;
	ShortReadMarker * shortReadArray, * shortReadMarker;
	IDnum shortReadCount, shortReadIndex;

	for(nodeID = 1; nodeID <= nodeCount(graph); nodeID++) {
		node = getNodeInGraph(graph, nodeID);
		if (node == NULL || getNodeLength(node) < minContigLength)
			continue;
		
		// Long reads
		for(marker = getMarker(node); marker != NULL; marker = getNextInNode(marker)) {
			readID = getPassageMarkerSequenceID(marker);
			if (readID < 0)
				readID = -readID;
			used[readID] = true;	
		}	

		// Short reads		
		if (!readStartsAreActivated(graph))
			continue;

		shortReadArray = getNodeReads(node, graph);
		shortReadCount = getNodeReadCount(node, graph);
		for (shortReadIndex = 0; shortReadIndex < shortReadCount; shortReadIndex++) {
			shortReadMarker = getShortReadMarkerAtIndex(shortReadArray, shortReadIndex);
			readID = getShortReadMarkerID(shortReadMarker);
			used[readID] = true;	
		}
		
		shortReadArray = getNodeReads(getTwinNode(node), graph);
		shortReadCount = getNodeReadCount(getTwinNode(node), graph);
		for (shortReadIndex = 0; shortReadIndex < shortReadCount; shortReadIndex++) {
			shortReadMarker = getShortReadMarkerAtIndex(shortReadArray, shortReadIndex);
			readID = getShortReadMarkerID(shortReadMarker);
			used[readID] = true;	
		}
	}

	for (readID = 1; readID <= sequenceCount(graph); readID++) 
		if (used[readID])
			res++;

	free(used);	

	return res;
}

void logFinalStats(Graph * graph, Coordinate minContigKmerLength, char *directory)
{
	char *logFilename =
	    mallocOrExit(strlen(directory) + 100, char);
	char *statsLine = 
	    mallocOrExit(5000, char);
	FILE *logFile;

	strcpy(logFilename, directory);
	strcat(logFilename, "/Log");
	logFile = fopen(logFilename, "a");

	if (logFile == NULL)
		exitErrorf(EXIT_FAILURE, true, "Could not write to %s",
		       logFilename);

	sprintf
	    (statsLine, "Final graph has %ld nodes and n50 of %lld, max %lld, total %lld, using %ld/%ld reads\n",
	     (long) nodeCount(graph), (long long) n50(graph), (long long) maxLength(graph),
	     (long long) totalAssemblyLength(graph), (long) usedReads(graph, minContigKmerLength),
	     (long) sequenceCount(graph));

	fprintf(logFile, "%s", statsLine);
	fprintf(stdout, "%s", statsLine);

	fclose(logFile);
	free(logFilename);
	free(statsLine);
}

void exportUnusedReads(Graph* graph, ReadSet * reads, Coordinate minContigKmerLength, char* directory) {
	char *outFilename =
	    mallocOrExit(strlen(directory) + 100, char);
	FILE * outfile;
	boolean * used = callocOrExit(sequenceCount(graph) + 1, boolean);
	IDnum nodeID, readID;
	Node * node;
	PassageMarker * marker;
	ShortReadMarker * shortReadArray, * shortReadMarker;
	IDnum shortReadCount, shortReadIndex;
	// Original
	IDnum numUnusedReads;
	// Original

	strcpy(outFilename, directory);
	strcat(outFilename, "/UnusedReads.fa");
	outfile = fopen(outFilename, "w");

	printf("Printing unused reads into %s\n", outFilename);
	
	for(nodeID = 1; nodeID <= nodeCount(graph); nodeID++) {
		node = getNodeInGraph(graph, nodeID);
		if (node == NULL || getNodeLength(node) < minContigKmerLength)
			continue;
		
		// Long reads
		for(marker = getMarker(node); marker != NULL; marker = getNextInNode(marker)) {
			readID = getPassageMarkerSequenceID(marker);
			if (readID < 0)
				readID = -readID;
			used[readID] = true;	
		}	

		// Short reads		
		if (!readStartsAreActivated(graph))
			continue;

		shortReadArray = getNodeReads(node, graph);
		shortReadCount = getNodeReadCount(node, graph);
		for (shortReadIndex = 0; shortReadIndex < shortReadCount; shortReadIndex++) {
			shortReadMarker = getShortReadMarkerAtIndex(shortReadArray, shortReadIndex);
			readID = getShortReadMarkerID(shortReadMarker);
			used[readID] = true;	
		}
		
		shortReadArray = getNodeReads(getTwinNode(node), graph);
		shortReadCount = getNodeReadCount(getTwinNode(node), graph);
		for (shortReadIndex = 0; shortReadIndex < shortReadCount; shortReadIndex++) {
			shortReadMarker = getShortReadMarkerAtIndex(shortReadArray, shortReadIndex);
			readID = getShortReadMarkerID(shortReadMarker);
			used[readID] = true;	
		}
	}

	for (readID = 1; readID <= sequenceCount(graph); readID++) 
		if (!used[readID]) {
			exportTightString(outfile, reads->tSequences[readID - 1], readID);	
			numUnusedReads++;
		}

	// Original
	//printf("%d\n", numUnusedReads);
	// Original

	free(outFilename);
	free(used);	
	fclose(outfile);
}

// Original
double getNodeDensity(Node * node)
{
	Coordinate nodeLength, nodeCoverage;
	
	nodeLength = getNodeLength(node);
	nodeCoverage = (getVirtualCoverage(node, 0) 
			+ getVirtualCoverage(node, 1));

	return nodeCoverage /(double) nodeLength;
}

int * makeDummySubgraphMask(Graph * graph)
{
	int lenSubgraphMask = 2 * nodeCount(graph) + 1;
	int *subgraphMask = callocOrExit(lenSubgraphMask, int);
	int index;

	for (index = 0; index < lenSubgraphMask; index++)
		subgraphMask[index] = 1;

	return subgraphMask;
}

int estimated_cov_multi(Graph * graph, int * subgraphMask, double expCovMulti[100])
{
	double histo[LEN_HISTO_X];
	Node *node;
	int index, ecmIndex = 0;
	double binWidth = 0.2;
	int bin = 0;
	double peakCov = 0.0, peakHeight = 0.0;
	double lastPeakHeight = 0.0;
	double SNratio = 10, thresNoiseHeight = 0.0;
	int noiseCount = 0, thresNoiseCount = 5;
	double thresMinPeak = 2.0;

	puts("Starting peak detection...");

	// Initialize expCovMulti[] and histo[]
	for (index = 0; index < 100; index++)
		expCovMulti[index] = -1;
	for (index = 0; index < LEN_HISTO_X; index++)
		histo[index] = 0.0;

	// Make histogram
	for (index = 1; index <= nodeCount(graph); index++) {
		if (subgraphMask[index + nodeCount(graph)] == 1) {
			node = getNodeInGraph(graph, index);
			node = getNodeInGraph(graph, index);
			if (node == NULL || getNodeLength(node) <= 0)
				continue;
			bin = (int) floor(getNodeDensity(node) / binWidth);
			if (bin >= LEN_HISTO_X - 1)
				bin = LEN_HISTO_X - 1;
			histo[bin] += getNodeLength(node);
		}
	}
	
	// Define length threshold of noise
	 // Skip index = 0 to avoid the influence of long reads
	for (index = LEN_HISTO_X - 2; index >= 1; index--) {
		if (histo[index] > peakHeight)
			peakHeight = histo[index];
	}
	thresNoiseHeight = peakHeight / (double) SNratio;
	peakHeight = 0.0;

	// Detect peaks
	for (index = LEN_HISTO_X - 2; index >= 1; index--) {
		if (histo[index] > thresNoiseHeight) {
			if (histo[index] > peakHeight 
			    && histo[index] > lastPeakHeight) {
				peakHeight = histo[index];
				peakCov = (double) (index + 0.5) * binWidth;
				noiseCount = 0;
				continue;
			}
			else {
				noiseCount++;
			}
		}
		else {
			lastPeakHeight = 0.0;
			noiseCount++;
		}

		if (peakHeight > 0.0 && noiseCount >= thresNoiseCount) {
			if (peakCov < thresMinPeak)
				break;

			expCovMulti[ecmIndex++] = peakCov;

			peakCov = 0.0;
			lastPeakHeight = peakHeight;
			peakHeight = 0.0;
			noiseCount = 0;
		}
	}

	// Output detedted peaks
	if (ecmIndex == 0) {
		puts("Error!! Couldn't detect any peaks");
		exit(1);
	}
	for (index = 0; index < ecmIndex; index++)
		printf("Detected Peak Coverage : %f\n", expCovMulti[index]);

	puts("Peak detection finished");

	return ecmIndex;
}

static void eliminateNullNodes(Graph * graph, int * subgraphMask)
{
	Node *node;
	int index;
	int lenSubgraphMask = 2 * nodeCount(graph) + 1;
	
	for (index = 0; index < lenSubgraphMask; index++) {
		node = getNodeInGraph(graph, index - nodeCount(graph));
		if (node == NULL || getNodeID(node) == 0)
			subgraphMask[index] = -2;
	}
}

static boolean checkLongReadExistence(Graph * graph)
{
	int index;

	for (index = 1; index <= nodeCount(graph); index++) {
		if (getMarker(getNodeInGraph(graph, index)) != NULL)
			return true;
	}
       
	return false;
}

static void depthFirstSearchSubgraph(int currentIndex, Graph * graph, int * subgraphMask)
{
	Arc *activeArc = NULL;
	int nextIndex = 0;

	if (subgraphMask[currentIndex + nodeCount(graph)] == 0) {
		// Mark "Visiting"
		subgraphMask[currentIndex + nodeCount(graph)] = 1;
		
		// Find "Unvisited" Node
		for (activeArc = getArc(getNodeInGraph(graph, currentIndex)); 
		     activeArc != NULL; activeArc = getNextArc(activeArc)) {
			nextIndex = getNodeID(getDestination(activeArc));
			if (subgraphMask[nextIndex] == 0) {
				depthFirstSearchSubgraph(nextIndex, graph, subgraphMask);
				depthFirstSearchSubgraph((nextIndex * -1), graph, subgraphMask);
			}
		}
	}
}

void resetUniqueness(Graph * graph)
{
	Node *node;
	int index;

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL)
			continue;
		setUniqueness(node, false);
	}
}

static void setSubgraphMask(int * subgraphMask, int lenSubgraphMask,
			    int before, int after)
{
	int index;

	for (index = 0; index < lenSubgraphMask; index++) {
		if (subgraphMask[index] == before)
			subgraphMask[index] = after;
	}
}

static int getUnvisitedNodeID(int * subgraphMask, int lenSubgraphMask)
{
	int index;

	for (index = 0; index < lenSubgraphMask; index++) {
		if (index == (lenSubgraphMask - 1) / 2)
			continue;
		if (subgraphMask[index] == 0)
			return index - (lenSubgraphMask - 1) / 2;
	}

	// Visited all nodes
	return 0;
}

static void shelveSubgraphMask(int * subgraphMask, int lenSubgraphMask,
			       int exception)
{
	int index;

	for (index = 0; index < lenSubgraphMask; index++) {
		if (subgraphMask[index] == exception)
			subgraphMask[index] = 0;
		else
			subgraphMask[index] += 100;
	}
}

static void unshelveSubgraphMask(int * subgraphMask, int lenSubgraphMask)
{
	int index;

	for (index = 0; index < lenSubgraphMask; index++) {
		if (subgraphMask[index] >= 50)
			subgraphMask[index] -= 100;
	}
}

static int estimated_cov_subgraph(Graph * graph, int * subgraphMask, double expCovPandS[2],
				  double rateChimericSubgraph)
{
	int nodeIndex;
	long int sumLenPrimary = 0, sumLenSecondary = 0, sumLenTotal;
	double perPrimary, perSecondary;
	Node *node;
	double cov = 0.0;
	double primary = expCovPandS[0], secondary = expCovPandS[1];

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		if (subgraphMask[nodeIndex + nodeCount(graph)] == 1) {
			node = getNodeInGraph(graph, nodeIndex);
			if (node == NULL)
				continue;
			cov = getNodeDensity(node);
			if (fabs(cov - primary) <= fabs(cov - secondary))
				sumLenPrimary += getNodeLength(node);
			else
				sumLenSecondary += getNodeLength(node);
		}
	}
       
	sumLenTotal = sumLenPrimary + sumLenSecondary;
	perPrimary = (double) sumLenPrimary / sumLenTotal;
	perSecondary = (double) sumLenSecondary / sumLenTotal;
	if (perSecondary <= rateChimericSubgraph) {
		// Non Chimeric Subgraph, belongs to Primary Species
		return 1;
	} else if (perPrimary <= rateChimericSubgraph) {
		// Non Chimeric Subgraph, belongs to Secondary Species
		return -1;
	} else {
		// Chimeric Subgraph
		return 2;
	}
}

static void forceSeparateChimericSubgraph(Graph * graph, int * subgraphMask,
					  double expCovPandS[2])
{
	int maskIndex;
	int lenSubgraphMask = nodeCount(graph) * 2 + 1;
	Node *node;
	double cov, primary = expCovPandS[0], secondary = expCovPandS[1];

	for (maskIndex = 0; maskIndex < lenSubgraphMask; maskIndex++) {
		if (subgraphMask[maskIndex] == 1) {
			node = getNodeInGraph(graph, maskIndex - nodeCount(graph));
			if (node == NULL)
				continue;
			cov = getNodeDensity(node);
			if (fabs(cov - primary) <= fabs(cov - secondary))
				subgraphMask[maskIndex] = 2;
			else
				subgraphMask[maskIndex] = -1;
		}
	}
}

static void judgeChimericSubgraph(Graph * graph, int * subgraphMask, double expCovPandS[2],
				  double rateChimericSubgraph, boolean discardChimericSubgraph)
{
	int nodeIndex;
	int judgeResult = 0;
	int lenSubgraphMask = nodeCount(graph) * 2 + 1;
	boolean flagVisitedAllNodes = false;
	int numSubgraph = 0, thresNumSubgraph = lenSubgraphMask;

	// Shelve not "Visiting" subgraphs
	shelveSubgraphMask(subgraphMask, lenSubgraphMask, 1);
	
	while (!flagVisitedAllNodes) {
		// Check Infinite Loop
		numSubgraph++;
		if (numSubgraph >= thresNumSubgraph) {
			puts("Resolving Repeat Error!! Infinite Loop");
			free(subgraphMask);
			exit(1);
		}
		
		// Choice starting unvisited node
		nodeIndex = getUnvisitedNodeID(subgraphMask, lenSubgraphMask);
		
		printf("nodeIndex = %d\n", nodeIndex);
		
		// Depth-first search (node & twin)
		depthFirstSearchSubgraph(nodeIndex, graph, subgraphMask);
		depthFirstSearchSubgraph((nodeIndex * -1), graph, subgraphMask);
		
		// Estimate exp_cov in the Subgraph
		judgeResult = estimated_cov_subgraph(graph, subgraphMask, expCovPandS,
						     rateChimericSubgraph);
		
		if (judgeResult == 1) {
			printf("NonChimeric Subgraph, belongs to Primary Species\n");
			setSubgraphMask(subgraphMask, lenSubgraphMask, 1, 2);
		}
		else if (judgeResult == -1) {
			printf("NonChimeric Subgraph, belongs to Secondary Species\n");
			setSubgraphMask(subgraphMask, lenSubgraphMask, 1, -1);
		}
		else {
			printf("Chimeric Subgraph!\n");
			if (discardChimericSubgraph)
				setSubgraphMask(subgraphMask, lenSubgraphMask, 1, -2);
			else
				forceSeparateChimericSubgraph(graph, subgraphMask, expCovPandS);
		}
	
		// Judge whether all nodes in Subgraph have visited or not
		if (getUnvisitedNodeID(subgraphMask, lenSubgraphMask) == 0)
			flagVisitedAllNodes = true;
	}

	// Unshelve not "Visiting" subgraphs
	unshelveSubgraphMask(subgraphMask, lenSubgraphMask);

	printf("Separating Chimeric Subgraphs Finished!\n");
}

static boolean checkPrimaryExpCovExistence(double expCovMulti[2], double expCovPandS[2])
{
	int index;
	double expCov, primary = expCovPandS[0], secondary = expCovPandS[1];

	if (secondary == -1)
		secondary = primary / (double) 2;

	for (index = 0; index < 2; index++) {
		expCov = expCovMulti[index];
		if (expCov == -1)
			break;
		if (fabs(expCov - primary) <= fabs(expCov - secondary)) {
			printf("Primary exp_cov (%f <-> %f), Exist\n",
			       primary, secondary);
			return true;
		}
	}
	
	printf("Primary exp_cov (%f <-> %f), NOT exist\n", primary, secondary);
	return false;
}

static boolean judgeSkip(Graph * graph, int * subgraphMask)
{
	int nodeIndex;
	int lenSubgraphMask = nodeCount(graph) * 2 + 1;
	int countSkip = 0;
	double insertLen = getInsertLength(graph, 1);
	double skipCandidateNodeLen = 0.0;
	Node *node = NULL;
		
	for (nodeIndex = 1; nodeIndex < nodeCount(graph); nodeIndex++) {
		if (subgraphMask[nodeIndex + nodeCount(graph)] == 1) {
			countSkip++;
			node = getNodeInGraph(graph, nodeIndex);
			skipCandidateNodeLen = getNodeLength(node);
		}
	}

	if (countSkip <= 1 && skipCandidateNodeLen < insertLen) {
		setSubgraphMask(subgraphMask, lenSubgraphMask, 1, -2);
		printf("Skipped\n");
		return true;
	}
	else
		return false;
}

static void printActiveNodes(int * subgraphMask, int lenSubgraphMask)
{
	int index;

	printf("Active Nodes : ");
	for (index = 0; index < lenSubgraphMask; index++) {
		if (subgraphMask[index] == 2)
			printf("%d ", index - (lenSubgraphMask - 1) / 2);
	}
       
	printf("\n");
}

void resolveRepeatOfAllSubgraphs(Graph * graph, ReadSet * reads, double expCovMulti[100],
				 boolean * dubious, boolean force_jumps, int argPebbleRounds, 
				 double rateChimericSubgraph, boolean discardChimericSubgraph,
				 double repeatNodeCovSD)
{
	int nodeIndex = 1, ecmIndex = 0;
	int numSubgraph = 0;
	int thresNumSubgraph = nodeCount(graph) * 2;
	int lenSubgraphMask = 2 * nodeCount(graph) + 1;
	int *subgraphMask = callocOrExit(lenSubgraphMask, int);
	double expCovSubgraph = 0.0;
	double expCovPandS[2];
	int numPeaks = 0;
	int countInterRepeatLoop = 0, thresInterRepeatLoop = 20;
	int pebbleRounds = argPebbleRounds;
	boolean flagLongRead = false, flagVisitedAllNodes = false;
	
	puts("\nResolving Repeats for each subgraph\n");

	// Eliminate NULL nodes
	eliminateNullNodes(graph, subgraphMask);

	// Check whether long reads are in the input sequences
	flagLongRead = checkLongReadExistence(graph);

	// Print Expected Coverages
	for (ecmIndex = 0; expCovMulti[ecmIndex] > 0.001; ecmIndex++)
		printf("Expected Coverage %d : %f\n", ecmIndex+1, expCovMulti[ecmIndex]);
	ecmIndex = 0;

	/*
	// Detect peaks from whole Graph
	setSubgraphMask(subgraphMask, lenSubgraphMask, 0, 1);
	estimated_cov_multi(graph, subgraphMask, expCovMulti);
	setSubgraphMask(subgraphMask, lenSubgraphMask, 1, 0);
	*/

	while (!flagVisitedAllNodes) {
		// Set expCovPandS
		if (expCovMulti[ecmIndex] != -1) {
			expCovPandS[0] = expCovMulti[ecmIndex++];
			expCovPandS[1] = expCovMulti[ecmIndex];
		}
		printf("\nPrimary exp_cov : %f\n", expCovPandS[0]);
			
		// Resolve repeats for each Subgraph
		while (true) {
			// Check Infinite Loop
			numSubgraph++;
			if (numSubgraph >= thresNumSubgraph) {
				puts("Resolving Repeat Error!! Infinite Loop");
				free(subgraphMask);
				exit(1);
			}
			
			// Choice starting unvisited node
			nodeIndex = getUnvisitedNodeID(subgraphMask, lenSubgraphMask);

			printf("nodeIndex = %d\n", nodeIndex);

			// Depth-first search (node & twin)
			depthFirstSearchSubgraph(nodeIndex, graph, subgraphMask);
			depthFirstSearchSubgraph((nodeIndex * -1), graph, subgraphMask);

			// Estimate the number of peaks
			numPeaks = estimated_cov_subgraph(graph, subgraphMask, expCovPandS,
							  rateChimericSubgraph);

			// Judge whether the Subgraph is chimeric or not
			if (numPeaks >= 2) {
				puts("Multiple Peaks Detected!");
				// Identify and Separate InterRepeats
				while (identifyAndSeparateInterRepeats(graph, expCovPandS, 
								       repeatNodeCovSD)) {
					// Check Infinite Loop
					if (countInterRepeatLoop++ >= thresInterRepeatLoop) {
						puts("Force-quitted to Identify InterRepeats");
						eliminateNullNodes(graph, subgraphMask);
						break;
						//puts("Identifying InterRepeat Error! Infinite Loop");
						//free(subgraphMask);
						//exit(1);
					}
					// Eliminate NULL nodes
					eliminateNullNodes(graph, subgraphMask);
				}
				// Judge whether each Subgraph is chimeric or not
				judgeChimericSubgraph(graph, subgraphMask, expCovPandS,
						      rateChimericSubgraph, discardChimericSubgraph);
			}
			else if (numPeaks == 1)
				setSubgraphMask(subgraphMask, lenSubgraphMask, 1, 2);
						
			// Judge whether all nodes in Subgraphs have visited or not
			if (getUnvisitedNodeID(subgraphMask, lenSubgraphMask) != 0) {
				printf("Unvisited Node : %d\n", 
				       getUnvisitedNodeID(subgraphMask, lenSubgraphMask));
				continue;
			}
			else {
				printf("\nGo to Assembly!\n");
				//printActiveNodes(subgraphMask, lenSubgraphMask);
				expCovSubgraph = expCovPandS[0];
				printf("exp_cov = %f\n", expCovSubgraph);
			}

			// -------------------- Assemble in the Subgraph --------------------
			// Judge unique or repeat
			identifyUniqueNodesSubgraph(graph, subgraphMask, 
						    isUniqueSolexaSubgraph, expCovSubgraph); 
			// Rock Band in the Subgraph
			if (flagLongRead)
				readCoherentSubgraph(graph, expCovSubgraph, reads, subgraphMask);
			// Pebble in the Subgraph
			for (pebbleRounds = argPebbleRounds; pebbleRounds > 0; pebbleRounds--)
				exploitShortReadPairs(graph, reads, dubious, force_jumps);
			// Print "Finished"
			printf("Subgraph Assembly Finished!\n\n");
			// ------------------------------------------------------------------

			// Eliminate NULL Nodes
			eliminateNullNodes(graph, subgraphMask);
			
			// Reset uniqueness
			resetUniqueness(graph);
			
			// Set "2" -> "-2", "-1" -> "0"
			setSubgraphMask(subgraphMask, lenSubgraphMask, 2, -2);
			setSubgraphMask(subgraphMask, lenSubgraphMask, -1, 0);
			
			// Judge whether all nodes in Graph have visited or not
			if (getUnvisitedNodeID(subgraphMask, lenSubgraphMask) == 0)
				flagVisitedAllNodes = true;
			break;
		}
	}
	
	// Resolved Successfully
	puts("Resolved Successfully!\n");
	free(subgraphMask);	
}
// Original
