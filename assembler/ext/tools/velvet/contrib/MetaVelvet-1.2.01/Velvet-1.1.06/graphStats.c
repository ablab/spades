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
#include <ctype.h>

#include "globals.h"
#include "graph.h"
#include "graphStats.h"
#include "readSet.h"
#include "tightString.h"
#include "passageMarker.h"
#include "concatenatedGraph.h"
#include "readCoherentGraph.h"
#include "fibHeap.h"
#include "utility.h"
#include "recycleBin.h"
#include "passageMarker.h"
// Original
#include "shortReadPairs.h"
// Original


static PassageMarkerList *copyMarkers(Node * node)
{
	PassageMarkerList *list = NULL;
	PassageMarkerList *_new;
	PassageMarkerI currentMarker;

	for (currentMarker = getMarker(node); currentMarker != NULL_IDX;
	     currentMarker = getNextInNode(currentMarker)) {
		_new = newPassageMarkerList(currentMarker, list);
		list = _new;
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

// Counts the number of markers for one node
int nodeGenomicMultiplicity(Node * node, IDnum firstStrain)
{
	int counter = 0;
	PassageMarkerI marker;

	if (node == NULL)
		return 0;

	for (marker = getMarker(node); marker != NULL_IDX;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) < firstStrain)
			counter++;

	return counter;
}

boolean isOnlyGenome(Node * node, IDnum firstStrain)
{
	PassageMarkerI marker;

	for (marker = getMarker(node); marker != NULL_IDX;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) >= firstStrain)
			return false;

	return true;
}

boolean isOnlyStrain(Node * node, IDnum firstStrain)
{
	PassageMarkerI marker;

	for (marker = getMarker(node); marker != NULL_IDX;
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

	if (getMarker(node) == NULL_IDX)
		return false;

	if (getAbsolutePassMarkerSeqID(getMarker(node)) >= firstStrain)
		return false;

	if (getNextInNode(getMarker(node)) != NULL_IDX)
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

	velvetLog("SNP\t%lld\t%ld\n", (long long) position, (long) sequence);

	return true;
}

void removeStrainMarkers(Node * node, IDnum firstStrain)
{
	PassageMarkerI marker;
	PassageMarkerI tmp = NULL_IDX;

	marker = getMarker(node);
	while (marker != NULL_IDX) {
		tmp = getNextInNode(marker);

		if (getAbsolutePassMarkerSeqID(marker) >= firstStrain)
			destroyPassageMarker(marker);
		marker = tmp;
	}

}

Coordinate commonLength(Node * node, IDnum firstStrain)
{
	PassageMarkerI marker = getMarker(node);
	int orig = 0;
	int strain = 0;

	while (marker != NULL_IDX) {
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

boolean isMixed(Node * node, IDnum firstStrain)
{
	return !isOnlyStrain(node, firstStrain)
	    && !isOnlyGenome(node, firstStrain);
}

int countLocalBreakpoints(PassageMarkerI marker, IDnum firstStrain)
{
	PassageMarkerI localMarker;
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
	for (localMarker = getMarker(localNode); localMarker != NULL_IDX;
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

IDnum genomeMarkerCount(Node * node, IDnum firstStrain)
{
	PassageMarkerI marker;
	IDnum counter = 0;

	for (marker = getMarker(node); marker != NULL_IDX;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) < firstStrain)
			counter++;

	return counter;
}

Coordinate readCoverage(Node * node)
{
	PassageMarkerI marker;
	Coordinate sum = 0;

	for (marker = getMarker(node); marker != NULL_IDX;
	     marker = getNextInNode(marker)) {
		if (getTwinMarker(marker) == NULL_IDX) {
			velvetLog("Node %li screwed up\n", (long) getNodeID(node));
			velvetLog("Sequence %li\n",
			       (long) getPassageMarkerSequenceID(marker));
			abort();
		}
		sum += getPassageMarkerLength(marker);
	}

	return sum;
}

Coordinate refReadCoverage(Node * node, IDnum firstStrain)
{
	PassageMarkerI marker;
	Coordinate sum = 0;

	for (marker = getMarker(node); marker != NULL_IDX;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) < firstStrain)
			sum += getPassageMarkerLength(marker);

	return sum;
}

Coordinate newReadCoverage(Node * node, IDnum firstStrain)
{
	PassageMarkerI marker;
	Coordinate sum = 0;

	for (marker = getMarker(node); marker != NULL_IDX;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) >= firstStrain) {
			sum += getPassageMarkerLength(marker);
			if (getPassageMarkerLength(marker) < 0)
				velvetLog("Bizarre marker %li at node %li\n",
				       (long) getPassageMarkerSequenceID(marker),
				       (long) getNodeID(node));
		}

	return sum;
}

static void printShortCounts(FILE * outfile, Node * node, Graph * graph, ReadSet * reads)
{
	IDnum counts[CATEGORIES];
	Category cat;
	IDnum shortReadIndex;
	IDnum readID;
	IDnum shortReadCount;
	ShortReadMarker *array;
	ShortReadMarker *marker;

	if (!readStartsAreActivated(graph)) {
		for (cat = 0; cat < CATEGORIES; cat++)
			velvetFprintf(outfile, "\tN/A");
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
		velvetFprintf(outfile, "\t%li", (long) counts[cat]);
}

void displayGeneralStatistics(Graph * graph, char *filename, ReadSet * reads)
{
	IDnum nodeIndex;
	Node *node;
	Category cat;
	FILE *outfile;

	outfile = fopen(filename, "w");
	if (outfile == NULL) {
		velvetLog("Couldn't open file %s, sorry\n", filename);
		return;
	} else
		velvetLog("Writing into stats file %s...\n", filename);

	velvetFprintf(outfile, "ID\tlgth\tout\tin\tlong_cov");

#ifndef SINGLE_COV_CAT
	for (cat = 0; cat < CATEGORIES; cat++) {
		velvetFprintf(outfile, "\tshort%i_cov", (int) (cat + 1));
		velvetFprintf(outfile, "\tshort%i_Ocov", (int) (cat + 1));
	}
#else
	velvetFprintf(outfile, "\tshort_cov");
#endif

	velvetFprintf(outfile, "\tlong_nb");
	for (cat = 0; cat < CATEGORIES; cat++) {
		velvetFprintf(outfile, "\tshort%i_nb", (int) (cat + 1));
	}

	velvetFprintf(outfile, "\n");

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex);
		if (node == NULL)
			continue;
		velvetFprintf
		    (outfile, "%ld\t%lld\t%i\t%i",
		     (long) getNodeID(node), (long long) getNodeLength(node), arcCount(node),
		     arcCount(getTwinNode(node)));

		if (getNodeLength(node) > 0) {
			velvetFprintf(outfile, "\t%f",
				readCoverage(node) /
				(double) getNodeLength(node));
#ifndef SINGLE_COV_CAT
			for (cat = 0; cat < CATEGORIES; cat++) {
				velvetFprintf(outfile, "\t%f",
					getVirtualCoverage(node, cat) /
					(double) getNodeLength(node));
				velvetFprintf(outfile, "\t%f",
					getOriginalVirtualCoverage(node, cat) /
					(double) getNodeLength(node));
			}
#else
			velvetFprintf(outfile, "\t%f",
					getVirtualCoverage(node) /
					(double) getNodeLength(node));
#endif
		} else {
			velvetFprintf(outfile, "\tInf");
#ifndef SINGLE_COV_CAT
			for (cat = 0; cat < CATEGORIES; cat++)
				velvetFprintf(outfile, "\tInf\tInf");
#else
			velvetFprintf(outfile, "\tInf");
#endif
		}

		velvetFprintf(outfile, "\t%li", (long) markerCount(node));
		printShortCounts(outfile, node, graph, reads); 

		velvetFprintf(outfile, "\n");
	}

	fclose(outfile);
}

void displayLocalBreakpoint(PassageMarkerI strainMarker,
			    IDnum firstStrain,
			    PassageMarkerI genomeMarker,
			    Node ** genomeDestination,
			    Node ** strainDestination, IDnum * counter,
			    IDnum nodeCount)
{
	boolean isTranslocation;
	PassageMarkerI marker;
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
//              velvetLog("Parallel paths\n");
		return;
	}

	destination2A = getNode(getNextInSequence(genomeMarker));

	if (destination2A == NULL)
		return;

	velvetLog("Lengths %lld %lld\n", (long long) getNodeLength(destinationA),
	       (long long) getNodeLength(destination2A));

	// Hop to another genomic node
//      if (getNodeLength(destinationA) > 24) {
	//velvetLog("wrong length %d %d\n", getNodeLength(destination) , getNodeID(destination));
//              return;
//      }

	destination =
	    getNode(getNextInSequence(getNextInSequence(strainMarker)));

	if (destination == NULL)
		return;

	// Eliminate those that point to uniquely strain sequences 
	if (nodeGenomicMultiplicity(destination, firstStrain) != 1) {
//              velvetLog("Multiple genome reads\n");
		return;
	}
	// Hop to another genomic node
//      if (getNodeLength(destination2A) != 24) {
	//velvetLog("wrong length 2\n");
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

//      velvetLog("Assigning %p and %p to %d\n", destination, destination2, localID);
	velvetLog("lengths %lld\t%lld\n", (long long) getNodeLength(destinationA),
	       (long long) getNodeLength(destination2A));

	// Detect translocation
	isTranslocation = true;
	for (marker = getMarker(destination); marker != NULL_IDX;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) ==
		    getAbsolutePassMarkerSeqID(genomeMarker)) {
			isTranslocation = false;
			break;
		}

	if (isTranslocation) {
		velvetLog("BREAK TRANS\t%ld\t%lld\t%lld\t%lld\n",
		       (long) getAbsolutePassMarkerSeqID(genomeMarker),
		       (long long) getPassageMarkerStart(genomeMarker),
		       (long long) getNodeLength(destinationA),
		       (long long) getNodeLength(destination2A));
		counter[2]++;
		return;
	}
	// Detect breakpoint
	velvetLog("BREAK INTRA\t%ld\t%lld\t%lld\t%lld\n",
	       (long) getAbsolutePassMarkerSeqID(genomeMarker),
	       (long long) getPassageMarkerStart(genomeMarker),
	       (long long) getNodeLength(destinationA), (long long) getNodeLength(destination2A));
	counter[1]++;

	// Check for inversion
	if (getPassageMarkerSequenceID(marker) !=
	    -getPassageMarkerSequenceID(genomeMarker))
		return;

//      velvetLog("potential!!\n");

	node1 = getTwinNode(destination);

	if (getNodeStatus(node1)) {
		node2 =
		    getTwinNode(genomeDestination
				[getNodeID(node1) + nodeCount]);
		if (getNodeStatus(node2))
			if (strainDestination[getNodeID(node2) + nodeCount]
			    == destination2) {
//                                      velvetLog("Safe\n");
				counter[1] -= 4;
				counter[0]++;
			} else;
//                              velvetLog("stopped 3\n");
		else;
//                      velvetLog("stopped 2\n");
	} else;
//              velvetLog("stopped 1\n");
}

PassageMarkerI genomeMarker(Node * node, IDnum firstStrain)
{
	PassageMarkerI marker;

	if (genomeMarkerCount(node, firstStrain) != 1)
		return NULL_IDX;

	for (marker = getMarker(node); marker != NULL_IDX;
	     marker = getNextInNode(marker))
		if (getAbsolutePassMarkerSeqID(marker) < firstStrain)
			return marker;

	return NULL_IDX;
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
	velvetFprintf(outfile, "> ARC from NODE %li", (long) getNodeID(getOrigin(arc)));
	velvetFprintf(outfile, "%s\n", str);
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
		velvetLog("Could not open %s, sorry\n", filename);
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

typedef struct mask_st Mask;

struct mask_st {
	Coordinate start;
	Coordinate finish;
	Mask* next;
}; 

static RecycleBin * maskMemory = NULL;

static Mask *allocateMask()
{
	if (maskMemory == NULL)
		maskMemory = newRecycleBin(sizeof(Mask), 10000);

	return (Mask *) allocatePointer(maskMemory);
}

static void deallocateMask(Mask * mask)
{
	deallocatePointer(maskMemory, mask);
}


static Mask * newMask(Coordinate position)
{
	Mask * mask = allocateMask();
	mask->start = position;
	mask->finish = position;
	mask->next = NULL;
	return mask;
}

static Mask * lowCoverageRegions(Coordinate * starts, Coordinate * stops, size_t length, IDnum cutoff, Coordinate nodeLength) {
	size_t indexStart = 0;
	size_t indexStop = 0;
	int currentValue = 0;
	Mask * regions = NULL;
	Mask * lastRegion = NULL;
	boolean openMask = false;

	while (indexStart < length && indexStop < length) {
		if (starts[indexStart] == stops[indexStop]) {
			indexStart++;
			indexStop++;
		} else if (starts[indexStart] < stops[indexStop]) {
			if (currentValue == cutoff - 1 && lastRegion) {
				lastRegion->finish = starts[indexStart] - 1;
				openMask = false;
			}
			currentValue++;
			indexStart++;
		} else {
			if (currentValue == cutoff && stops[indexStop] != nodeLength) {
				if (regions) {
					lastRegion->next = newMask(stops[indexStop]);
					lastRegion = lastRegion->next;
				} else { 
					regions = newMask(stops[indexStop]);
					lastRegion = regions;
				}
				openMask = true;
			}
			currentValue--;
			indexStop++;
		}
	}

	while (indexStart < length) {
		if (currentValue == cutoff - 1 && lastRegion) {
			lastRegion->finish = starts[indexStart] - 1;
			openMask = false;
		} else if (currentValue >= cutoff) {
			break;
		}
		currentValue++;
		indexStart++;
	}

	while (indexStop < length) {
		if (currentValue == cutoff + 1) {
			if (regions) {
				lastRegion->next = newMask(stops[indexStop]);
				lastRegion = lastRegion->next;
			} else { 
				regions = newMask(stops[indexStop]);
				lastRegion = regions;
			}
			openMask = true;
		} else if (currentValue < cutoff)
			break;
		currentValue--;
		indexStop++;
	}

	if (openMask)
		lastRegion->finish = nodeLength;

	free(starts);
	free(stops);
	return regions;
}

static int compareCoords(const void * A, const void * B) {
	Coordinate * a_p = (Coordinate *) A;
	Coordinate * b_p = (Coordinate *) B;
	Coordinate a = * a_p;
	Coordinate b = * b_p;
	if (a < b)
		return -1;
	else if (a > b)
		return 1;
	else
		return 0;
}

static void sortCoords(Coordinate * array, IDnum length) {
	qsort(array, (size_t) length, sizeof(Coordinate), compareCoords);
} 

static void getShortReadCoords(Coordinate * starts, Coordinate * stops, Node * node, Graph * graph, ShortLength * readLengths) {
	ShortReadMarker * markers = getNodeReads(node, graph);
	ShortReadMarker* marker;
	IDnum index;
	for (index = 0; index < getNodeReadCount(node, graph); index++) {
		marker = getShortReadMarkerAtIndex(markers, index);
		starts[index] = getShortReadMarkerPosition(marker);
		stops[index] = starts[index] - getShortReadMarkerOffset(marker) + readLengths[getShortReadMarkerID(marker) - 1] - 1;
	}
}

static void getShortReadTwinCoords(Coordinate * starts, Coordinate * stops, Node * node, Graph * graph, ShortLength * readLengths, IDnum offset) {
	ShortReadMarker * markers = getNodeReads(getTwinNode(node), graph);
	ShortReadMarker* marker;
	IDnum index;
	for (index = 0; index < getNodeReadCount(getTwinNode(node), graph); index++) {
		marker = getShortReadMarkerAtIndex(markers, index);
		stops[index + offset] = getNodeLength(node) - 1 - getShortReadMarkerPosition(marker);
		starts[index + offset] = stops[index + offset] + getShortReadMarkerOffset(marker) - readLengths[getShortReadMarkerID(marker) - 1] + 1;
	}
}

static void getLongReadCoords(Coordinate * starts, Coordinate * stops, Node * node, Graph * graph, ReadSet * reads, IDnum offset) {
	PassageMarkerI marker;
	IDnum index = offset;

	for (marker = getMarker(node); marker; marker = getNextInNode(marker)) {
		if (reads->categories[getAbsolutePassMarkerSeqID(marker) - 1] != REFERENCE) {
			starts[index] = getStartOffset(marker);
			stops[index++] = getNodeLength(node) - 1 - getFinishOffset(marker);
		} else {
			starts[index] = -5;
			stops[index++] = -10;
		}
	}
}

static Mask * findLowCoverageRegions(Node * node, Graph * graph, IDnum cutoff, ReadSet * reads, ShortLength * readLengths) {
	// Fill arrays
	IDnum nodeReads = getNodeReadCount(node, graph);
	IDnum twinReads =  getNodeReadCount(getTwinNode(node), graph);
	IDnum longReads = markerCount(node);
	IDnum length = nodeReads + twinReads + longReads;
	Coordinate * starts = callocOrExit(length, Coordinate);
	Coordinate * stops = callocOrExit(length, Coordinate);
	getShortReadCoords(starts, stops, node, graph, readLengths);
	getShortReadTwinCoords(starts, stops, node, graph, readLengths, nodeReads);
	getLongReadCoords(starts, stops, node, graph, reads, nodeReads + twinReads);

	// Sort arrays
	sortCoords(starts, length);
	sortCoords(stops, length);

	// Go through array
	return lowCoverageRegions(starts, stops, length, cutoff, getNodeLength(node));
}

static void exportLongNodeSequence(FILE * outfile, Node * node, Graph * graph, ReadSet * reads, ShortLength * readLengths, IDnum cutoff) {
	TightString *tString;
	Coordinate position;
	char nucleotide;
	int WORDLENGTH = getWordLength(graph);
	GapMarker *gap;
	IDnum nodeIndex = getNodeID(node);
	Mask * mask = NULL;
	Mask * next;

	if (readStartsAreActivated(graph) && cutoff > 0)
		mask = findLowCoverageRegions(node, graph, cutoff, reads, readLengths);

	tString = expandNode(node, WORDLENGTH);
	velvetFprintf(outfile, ">NODE_%ld_length_%lld_cov_%f\n",
		(long) nodeIndex, (long long) getNodeLength(node),
		(getTotalCoverage(node) + readCoverage(node)) /
		(float) getNodeLength(node));

	gap = getGap(node, graph);
	for (position = 0; position < WORDLENGTH; position++) {
		if (position % 60 == 0 && position > 0)
			velvetFprintf(outfile, "\n"); 
		nucleotide = getNucleotideChar(position, tString);
		velvetFprintf(outfile, "%c", nucleotide);
	}

	gap = getGap(node, graph);
	for (; position < getLength(tString); position++) {
		if (position % 60 == 0)
			velvetFprintf(outfile, "\n");

		while (gap
		    && position - WORDLENGTH + 1 >=
		    getGapFinish(gap))
			gap = getNextGap(gap);

		while (mask
		    && position - WORDLENGTH + 1 >=
		    mask->finish) {
			next = mask->next;
			deallocateMask(mask);
			mask = next;	
		}

		if (gap
		    && position - WORDLENGTH + 1 >=
		    getGapStart(gap)) {
			velvetFprintf(outfile, "N");
		} else if (mask &&
			position - WORDLENGTH + 1 >=
			mask->start) {
			nucleotide =
			    getNucleotideChar(position, tString);
			velvetFprintf(outfile, "%c", tolower(nucleotide));

		} else {
			nucleotide =
			    getNucleotideChar(position, tString);
			velvetFprintf(outfile, "%c", nucleotide);
		}
	}
	velvetFprintf(outfile, "\n");
	destroyTightString (tString);
}

void exportLongNodeSequences(char *filename, Graph * graph,
			     Coordinate minLength, ReadSet * reads, ShortLength * readLengths, IDnum minCov)
{
	FILE *outfile = fopen(filename, "w");
	IDnum nodeIndex;
	Node *node;
	//double sensitivity, specificity;

	if (outfile == NULL) {
		velvetLog("Could not write into %s, sorry\n", filename);
		return;
	} else {
		velvetLog("Writing contigs into %s...\n", filename);
	}

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex);

		if (node == NULL || getNodeLength(node) < minLength)
			continue;

		exportLongNodeSequence(outfile, node, graph, reads, readLengths, minCov);
	}

	if (maskMemory)
		destroyRecycleBin(maskMemory);
	maskMemory = NULL;

	fclose(outfile);
}

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
		velvetLog("EMPTY GRAPH\n");
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

	//velvetLog("Measuring median coverage depth...\n");

	if (nodeCount(graph) == 0) {
		velvetLog("EMPTY GRAPH\n");
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
		  //velvetLog("Median coverage depth = %f\n", getTotalCoverage(node) / (double) getNodeLength(node));
		  //velvetFprintf(logFile, "Median coverage depth = %f\n", getTotalCoverage(node) / (double) getNodeLength(node));
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

static boolean terminalReferenceMarker(Node * node, ReadSet * reads) {
	PassageMarkerI marker;

	for (marker = getMarker(node); marker; marker = getNextInNode(marker))
		if (reads->categories[getAbsolutePassMarkerSeqID(marker) - 1] == REFERENCE
		    && (!getNextInSequence(marker)
		        || !getPreviousInSequence(marker)))
			return true;

	return false;
}

static boolean hasReferenceMarker(Node * node, ReadSet * reads) {
	PassageMarkerI marker;
	
	for (marker = getMarker(node); marker != NULL_IDX; marker = getNextInNode(marker))
		if (reads->categories[getAbsolutePassMarkerSeqID(marker) - 1] == REFERENCE)
			return true;

	return false;
}

inline static void
destroyNodePassageMarkers(Graph *graph,
			     Node* node)
{
	PassageMarkerI marker;

	while ((marker = getMarker(node)) != NULL_IDX) {
		if (!isInitial(marker) && !isTerminal(marker))
			deleteNextPassageMarker(getPreviousInSequence(marker), graph);
		destroyPassageMarker(marker);
	}
}

inline static void
removeNodeAndDenounceDubiousReads(Graph *graph,
				  Node  *node,
				  boolean denounceReads,
				  boolean *res,
				  Coordinate minLength,
				  FILE *outfile)
{
	if (denounceReads) {
		ShortReadMarker *nodeArray;
		ShortReadMarker *shortMarker;
		IDnum maxIndex;
		IDnum index;
		IDnum readID;

		nodeArray = getNodeReads(node, graph);
		maxIndex = getNodeReadCount(node, graph);
		for (index = 0; index < maxIndex; index++) {
			shortMarker = getShortReadMarkerAtIndex(nodeArray, index);
			readID = getShortReadMarkerID(shortMarker);
			if (readID > 0)
				res[readID - 1] = true;
			else
				res[-readID - 1] = true;
		}

		nodeArray = getNodeReads(getTwinNode(node), graph);
		maxIndex = getNodeReadCount(getTwinNode(node), graph);
		for (index = 0; index < maxIndex; index++) {
			shortMarker = getShortReadMarkerAtIndex(nodeArray, index);
			readID = getShortReadMarkerID(shortMarker);
			if (readID > 0)
				res[readID - 1] = true;
			else
			  res[-readID - 1] = true;
		}
	}

	destroyNodePassageMarkers(graph, node);

	if (outfile != NULL && getNodeLength(node) > minLength)
		exportLongNodeSequence(outfile, node, graph, NULL, NULL, -1);

	destroyNode(node, graph);
}

boolean *removeLowCoverageNodesAndDenounceDubiousReads(Graph * graph,
						       double minCov,
						       ReadSet * reads,
						       boolean _export,
						       Coordinate minLength,
						       char *filename)
{
	IDnum index;
	Node *node;
	boolean denounceReads = readStartsAreActivated(graph);
	boolean *res = NULL; 
	FILE * outfile = NULL;

	velvetLog("Removing contigs with coverage < %f...\n", minCov);
		
	if (denounceReads)
		res = callocOrExit(sequenceCount(graph), boolean);
		
	if (_export) {
		outfile = fopen(filename, "w");

		if (outfile == NULL) {
			velvetLog("Could not write into %s, sorry\n", filename);
			return res;
		} else {
			velvetLog("Writing contigs into %s...\n", filename);
		}
	}


	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (getNodeLength(node) == 0)
			continue;

		if (getTotalCoverage(node) / getNodeLength(node) < minCov 
		    && !hasReferenceMarker(node, reads))
			removeNodeAndDenounceDubiousReads(graph,
							  node,
							  denounceReads,
							  res,
							  minLength,
							  outfile);
	}

	concatenateGraph(graph);

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (getNodeLength(node) == 0)
			continue;

		if (getTotalCoverage(node) / getNodeLength(node) < minCov 
		    && !terminalReferenceMarker(node, reads))
			removeNodeAndDenounceDubiousReads(graph,
							  node,
							  denounceReads,
							  res,
							  minLength,
							  outfile);
	}

	if (_export)
		fclose(outfile);

	concatenateGraph(graph);
	return res;
}

static Coordinate getLongCoverage(Node * node) {
	PassageMarkerI marker;
	Coordinate total = 0;

	for (marker = getMarker(node); marker; marker = getNextInNode(marker))
		total += getPassageMarkerLength(marker);
	
	return total;
}

void removeLowCoverageReferenceNodes(Graph * graph, double minCov, double minLongCov, ReadSet * reads)
{
	IDnum index;
	Node *node;

	velvetLog("Removing reference contigs with coverage < %f...\n", minCov);

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (getNodeLength(node) == 0)
			continue;

		if ((getTotalCoverage(node) / getNodeLength(node) < minCov 
		    || getLongCoverage(node) / getNodeLength(node) < minLongCov)
		    && hasReferenceMarker(node, reads)) {
			destroyNodePassageMarkers(graph, node);

			destroyNode(node, graph);
		}
	}

	concatenateGraph(graph);
}

void removeLowLongCoverageNodesAndDenounceDubiousReads(Graph * graph,
						       double minCov,
						       ReadSet * reads,
						       boolean * res,
						       boolean _export,
						       Coordinate minLength,
						       char *filename)
{
	IDnum index;
	Node *node;
	boolean denounceReads = readStartsAreActivated(graph);
	FILE * outfile = NULL;

	if (minCov < 0)
		return;

	velvetLog("Removing contigs with coverage < %f...\n", minCov);
		
	if (_export) {
		outfile = fopen(filename, "a");

		if (outfile == NULL) {
			velvetLog("Could not write into %s, sorry\n", filename);
			return;
		} else {
			velvetLog("Writing contigs into %s...\n", filename);
		}
	}

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (getNodeLength(node) == 0)
			continue;

		if (getLongCoverage(node) / getNodeLength(node) < minCov 
		    && !hasReferenceMarker(node, reads))
			removeNodeAndDenounceDubiousReads(graph,
							  node,
							  denounceReads,
							  res,
							  minLength,
							  outfile);
	}

	concatenateGraph(graph);

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (getNodeLength(node) == 0)
			continue;

		if (getLongCoverage(node) / getNodeLength(node) < minCov 
		    && !terminalReferenceMarker(node, reads))
			removeNodeAndDenounceDubiousReads(graph,
							  node,
							  denounceReads,
							  res,
							  minLength,
							  outfile);
	}

	if (_export)
		fclose(outfile);

	concatenateGraph(graph);
}

void removeHighCoverageNodes(Graph * graph, double maxCov, boolean _export, Coordinate minLength, char *filename)
{
	IDnum index;
	Node *node;
	FILE * outfile = NULL;

	if (maxCov < 0)
		return;

	velvetLog("Applying an upper coverage cutoff of %f...\n", maxCov);
		
	if (_export) {
		outfile = fopen(filename, "w");

		if (outfile == NULL) {
			velvetLog("Could not write into %s, sorry\n", filename);
			return;
		} else {
			velvetLog("Writing contigs into %s...\n", filename);
		}
	}

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (getNodeLength(node) > 0
		    && getTotalCoverage(node) / getNodeLength(node) > maxCov) {
			destroyNodePassageMarkers(graph, node);

			if (_export && getNodeLength(node) > minLength) 
				exportLongNodeSequence(outfile, node, graph, NULL, NULL, -1);

			destroyNode(node, graph);
		}
	}

	if (_export)
		fclose(outfile);

	concatenateGraph(graph);
}

static void exportAMOSLib(FILE * outfile, Graph * graph, Category cat)
{
	Coordinate distance = getInsertLength(graph, cat * 2);
	double variance = getInsertLength_var(graph, cat * 2);

	if (distance == -1)
		return;

	velvetFprintf(outfile, "{LIB\n");
	velvetFprintf(outfile, "iid:%d\n", (int) (cat + 1));
	velvetFprintf(outfile, "{DST\n");
	velvetFprintf(outfile, "mea:%lld\n", (long long) distance);
	velvetFprintf(outfile, "std:%lld\n", (long long) sqrt(variance));
	velvetFprintf(outfile, "}\n");
	velvetFprintf(outfile, "}\n");
}

static void exportAMOSMarker(FILE * outfile, PassageMarkerI marker,
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

	velvetFprintf(outfile, "{TLE\n");
	velvetFprintf(outfile, "src:%li\n", (long) getAbsolutePassMarkerSeqID(marker));
	if (getStartOffset(marker) > start)
		velvetFprintf(outfile, "off:%lld\n",
			(long long) (getStartOffset(marker) - start));
	else
		velvetFprintf(outfile, "off:0\n");
	velvetFprintf(outfile, "clr:%lld,%lld\n", (long long) sequenceStart, (long long) sequenceFinish);
	velvetFprintf(outfile, "}\n");
}

static void exportAMOSShortMarker(FILE * outfile, ShortReadMarker * marker,
				  ReadSet * reads, Coordinate start,
				  Coordinate finish)
{
	Coordinate offset =
	    getShortReadMarkerPosition(marker) -
	    getShortReadMarkerOffset(marker);
	TightString *sequence =
	    getTightStringInArray (reads->tSequences, getShortReadMarkerID(marker) - 1);

	if (getShortReadMarkerPosition(marker) == -1)
		return;

	if (offset >= finish || offset + getLength(sequence) < start)
		return;

	velvetFprintf(outfile, "{TLE\n");
	velvetFprintf(outfile, "src:%li\n", (long) getShortReadMarkerID(marker));
	velvetFprintf(outfile, "off:%lld\n", (long long) (offset - start));
	velvetFprintf(outfile, "clr:0,%lld\n", (long long) getLength(sequence));
	velvetFprintf(outfile, "}\n");
}

static void exportAMOSReverseShortMarker(FILE * outfile,
					 ShortReadMarker * marker,
					 Coordinate nodeLength,
					 int wordShift, ReadSet * reads,
					 Coordinate start,
					 Coordinate finish)
{
	TightString *sequence =
	    getTightStringInArray (reads->tSequences, getShortReadMarkerID(marker) - 1);

	Coordinate offset =
	    nodeLength - getShortReadMarkerPosition(marker) +
	    getShortReadMarkerOffset(marker) - getLength(sequence) +
	    wordShift;

	if (getShortReadMarkerPosition(marker) == -1)
		return;

	if (offset >= finish || offset + getLength(sequence) < start)
		return;

	velvetFprintf(outfile, "{TLE\n");
	velvetFprintf(outfile, "src:%li\n", (long) getShortReadMarkerID(marker));
	velvetFprintf(outfile, "off:%lld\n", (long long) (offset - start));
	velvetFprintf(outfile, "clr:%lld,0\n", (long long) getLength(sequence));
	velvetFprintf(outfile, "}\n");
}

static void exportAMOSContig(FILE * outfile, ReadSet * reads, Node * node,
			     Graph * graph, Coordinate contigStart,
			     Coordinate contigFinish, IDnum iid,
			     IDnum internalIndex)
{
	Coordinate start;
	char str[100];
	PassageMarkerI marker;
	ShortReadMarker *shortMarkerArray, *shortMarker;
	Coordinate index, maxIndex;
	int wordShift = getWordLength(graph) - 1;
	char *string = expandNodeFragment(node, contigStart, contigFinish,
					  getWordLength(graph));
	Coordinate length = contigFinish - contigStart + wordShift;

	velvetFprintf(outfile, "{CTG\n");
	velvetFprintf(outfile, "iid:%li\n", (long) iid);
	velvetFprintf(outfile, "eid:%li-%li\n", (long) getNodeID(node), (long) internalIndex);

	velvetFprintf(outfile, "seq:\n");
	for (start = 0; start <= length; start += 60) {
		strncpy(str, &(string[start]), 60);
		str[60] = '\0';
		velvetFprintf(outfile, "%s\n", str);
	}
	velvetFprintf(outfile, ".\n");

	velvetFprintf(outfile, "qlt:\n");
	for (start = 0; start <= length; start += 60) {
		strncpy(str, &(string[start]), 60);
		str[60] = '\0';
		velvetFprintf(outfile, "%s\n", str);
	}
	velvetFprintf(outfile, ".\n");

	free(string);

	for (marker = getMarker(node); marker != NULL_IDX;
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

	velvetFprintf(outfile, "}\n");
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

	velvetFprintf(outfile, "{SCF\n");
	velvetFprintf(outfile, "eid:%li\n", (long) getNodeID(node));
	for (gap = getGap(node, graph); gap; gap = getNextGap(gap)) {
		finish = getGapStart(gap);
		velvetFprintf(outfile, "{TLE\n");
		velvetFprintf(outfile, "off:%lld\n", (long long) start);
		velvetFprintf(outfile, "clr:0,%lld\n",
			(long long) (finish - start + (long long) wordShift));
		velvetFprintf(outfile, "src:%li\n", (long) contigIndex++);
		velvetFprintf(outfile, "}\n");
		start = getGapFinish(gap);
	}
	finish = getNodeLength(node);
	velvetFprintf(outfile, "{TLE\n");
	velvetFprintf(outfile, "off:%lld\n", (long long) start);
	velvetFprintf(outfile, "clr:0,%lld\n", (long long) (finish - start));
	velvetFprintf(outfile, "src:%li\n", (long) contigIndex++);
	velvetFprintf(outfile, "}\n");

	velvetFprintf(outfile, "}\n");
}

static void exportAMOSRead(FILE * outfile, TightString * tString,
			   IDnum index, IDnum frg_index)
{
	Coordinate start, finish;
	char str[100];

	velvetFprintf(outfile, "{RED\n");
	velvetFprintf(outfile, "iid:%li\n", (long) index);
	velvetFprintf(outfile, "eid:%li\n", (long) index);
	if (frg_index > 0)
		velvetFprintf(outfile, "frg:%li\n", (long) frg_index);

	velvetFprintf(outfile, "seq:\n");
	start = 0;
	while (start <= getLength(tString)) {
		finish = start + 60;
		readTightStringFragment(tString, start, finish, str);
		velvetFprintf(outfile, "%s\n", str);
		start = finish;
	}
	velvetFprintf(outfile, ".\n");

	velvetFprintf(outfile, "qlt:\n");
	start = 0;
	while (start <= getLength(tString)) {
		finish = start + 60;
		readTightStringFragment(tString, start, finish, str);
		velvetFprintf(outfile, "%s\n", str);
		start = finish;
	}
	velvetFprintf(outfile, ".\n");

	velvetFprintf(outfile, "}\n");
}

void exportAMOSContigs(char *filename, Graph * graph,
		       Coordinate cutoff_length, ReadSet * reads)
{
	IDnum index;
	Category cat;
	Node *node;
	FILE *outfile;

	velvetLog("Writing into AMOS file %s...\n", filename);
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
			velvetFprintf(outfile, "{FRG\n");
			velvetFprintf(outfile, "lib:%d\n",
				(int) ((reads->categories[index - 1] / 2) + 1));
			velvetFprintf(outfile, "rds:%li,%li\n", (long) index,
				(long) index + 1);
			velvetFprintf(outfile, "eid:%li\n", (long) index);
			velvetFprintf(outfile, "iid:%li\n", (long) index);
			velvetFprintf(outfile, "typ:I\n");
			velvetFprintf(outfile, "}\n");
			index++;
		}
	}

	for (index = 1; index <= reads->readCount; index++) {
		if (reads->categories[index - 1] % 2 != 0 &&
		    getInsertLength(graph,
				    reads->categories[index - 1]) >= 0) {
			exportAMOSRead(outfile,
				       getTightStringInArray(reads->tSequences, index - 1), index,
				       index);
			index++;
			exportAMOSRead(outfile,
				       getTightStringInArray(reads->tSequences, index - 1), index,
				       index - 1);
		} else {
			exportAMOSRead(outfile,
				       getTightStringInArray(reads->tSequences, index - 1), index,
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
	PassageMarkerI marker;
	ShortReadMarker * shortReadArray, * shortReadMarker;
	IDnum shortReadCount, shortReadIndex;

	for(nodeID = 1; nodeID <= nodeCount(graph); nodeID++) {
		node = getNodeInGraph(graph, nodeID);
		if (node == NULL || getNodeLength(node) < minContigLength)
			continue;
		
		// Long reads
		for(marker = getMarker(node); marker != NULL_IDX; marker = getNextInNode(marker)) {
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

	velvetFprintf(logFile, "%s", statsLine);
	velvetFprintf(stdout, "%s", statsLine);

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
	PassageMarkerI marker;
	ShortReadMarker * shortReadArray, * shortReadMarker;
	IDnum shortReadCount, shortReadIndex;

	strcpy(outFilename, directory);
	strcat(outFilename, "/UnusedReads.fa");
	outfile = fopen(outFilename, "w");

	velvetLog("Printing unused reads into %s\n", outFilename);

	for(nodeID = 1; nodeID <= nodeCount(graph); nodeID++) {
		node = getNodeInGraph(graph, nodeID);
		if (node == NULL || getNodeLength(node) < minContigKmerLength)
			continue;
		
		// Long reads
		for(marker = getMarker(node); marker != NULL_IDX; marker = getNextInNode(marker)) {
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
		if (!used[readID])
			exportTightString(outfile, getTightStringInArray(reads->tSequences, readID - 1), readID);

	free(outFilename);
	free(used);	
	fclose(outfile);
}

static IDnum getReferenceCount(ReadSet * reads) {
	IDnum index;

	for (index = 0; index < reads->readCount; index++) 
		if (reads->categories[index] != REFERENCE)
			break;

	return index;
}

//////////////////////////////////////////////////////////////////////////
// Reference identifiers
//////////////////////////////////////////////////////////////////////////

typedef struct referenceCoord_st ReferenceCoord;

struct referenceCoord_st {
	IDnum start;
	IDnum finish;
	char * name;
	boolean positive_strand;
} ATTRIBUTE_PACKED;

static ReferenceCoord * collectReferenceCoords(char * sequencesFilename, IDnum referenceCount) {
	FILE * file = fopen(sequencesFilename, "r");
	char line[MAXLINE];
	char name[5000];
	Coordinate start, finish;
	long long longlongvar;
	IDnum refIndex = 0;
	ReferenceCoord * refCoords = callocOrExit(referenceCount, ReferenceCoord);
	int i;

	while (fgets(line, MAXLINE, file)) {
		if (line[0] == '>') {
			if (strchr(line, ':')) {
				sscanf(strtok(line, ":-\r\n"), ">%s", name);
				sscanf(strtok(NULL, ":-\r\n"), "%lli", &longlongvar);
				start = longlongvar;
				sscanf(strtok(NULL, ":-\r\n"), "%lli", &longlongvar);
				finish = longlongvar;
				refCoords[refIndex].name = callocOrExit(strlen(name) + 1, char);  
				if (start <= finish) {
					strcpy(refCoords[refIndex].name, name);
					refCoords[refIndex].start = start;
					refCoords[refIndex].finish = finish;
					refCoords[refIndex].positive_strand = true;
				} else {
					strcpy(refCoords[refIndex].name, name);
					refCoords[refIndex].start = finish;
					refCoords[refIndex].finish = start;
					refCoords[refIndex].positive_strand = false;
				}
			} else {
				for (i = strlen(line) - 1;
				     i >= 0 && (line[i] == '\n' || line[i] == '\r'); i--) {
					line[i] = '\0';
				}

				strcpy(name, line + 1);
				refCoords[refIndex].name = callocOrExit(strlen(name) + 1, char);  
				strcpy(refCoords[refIndex].name, name);
				refCoords[refIndex].start = 1;
				refCoords[refIndex].finish = -1;
				refCoords[refIndex].positive_strand = true;
			}
			if (++refIndex == referenceCount)
				break;	
		}
	}
	
	fclose(file);
	return refCoords;
}

typedef struct refMap_st {
	IDnum start;
	IDnum finish;
	IDnum refID;
	IDnum refStart;
	IDnum refFinish;
} ATTRIBUTE_PACKED ReferenceMapping; 

static int compareReferenceMappings(const void * A, const void * B) {
	ReferenceMapping * refMapA = (ReferenceMapping *) A;
	ReferenceMapping * refMapB = (ReferenceMapping *) B;
	
	if (refMapA->start < refMapB->start)
		return -1;
	else if (refMapA->start == refMapB->start)
		return 0;
	else 
		return 1;
}

static void initializeReferenceMapping(ReferenceMapping * refMap, PassageMarkerI marker, Node * node) {
	refMap->start = getStartOffset(marker);
	refMap->finish = getNodeLength(node) - getFinishOffset(marker); 
	refMap->refID = getPassageMarkerSequenceID(marker);
	refMap->refStart = getPassageMarkerStart(marker);
	refMap->refFinish = getPassageMarkerFinish(marker);
}

static void velvetFprintfReferenceMapping(FILE * file, ReferenceMapping * mapping, ReferenceCoord * refCoords, int wordLength) {
	ReferenceCoord * refCoord;
	Coordinate start, finish;

	if (mapping->refID > 0) 
		refCoord = &refCoords[mapping->refID - 1];
	else
		refCoord = &refCoords[-mapping->refID - 1];

	if (mapping->refID > 0) {
		if (refCoord->positive_strand) {
			start = refCoord->start + mapping->refStart;
			finish = refCoord->start + mapping->refFinish + wordLength - 2;
		} else {
			start = refCoord->finish - mapping->refStart + wordLength - 1;
			finish = refCoord->finish - mapping->refFinish + 1;
		}
	} else {
		if (refCoord->positive_strand) {
			start = refCoord->start + mapping->refStart + wordLength - 1;
			finish = refCoord->start + mapping->refFinish + 1;
		} else {
			start = refCoord->finish - mapping->refStart; 
			finish = refCoord->finish - mapping->refFinish + wordLength;  
		}
	}
		
	velvetFprintf(file, "%lli\t%lli\t%s\t%lli\t%lli\n",
		(long long) mapping->start + 1, (long long) mapping->finish + wordLength - 1, 
		refCoord->name, (long long) start, (long long) finish);
}

static void exportLongNodeMapping(FILE * outfile, Node * node, ReadSet * reads, ReferenceCoord * refCoords, int wordLength) {
	PassageMarkerI marker;
	ReferenceMapping * referenceMappings;
	IDnum index;
	IDnum referenceCount = 0;

	// Count reference sequences
	for (marker = getMarker(node); marker != NULL_IDX; marker = getNextInNode(marker))
		if (reads->categories[getAbsolutePassMarkerSeqID(marker) - 1] == REFERENCE)
			referenceCount++;

	// Header
	velvetFprintf(outfile, ">contig_%li\n", (long) getNodeID(node));

	// Create table
	referenceMappings = callocOrExit(referenceCount, ReferenceMapping);	

	// Initialize table
	referenceCount = 0;
	for (marker = getMarker(node); marker != NULL_IDX; marker = getNextInNode(marker))
		if (reads->categories[getAbsolutePassMarkerSeqID(marker) - 1] == REFERENCE)
			initializeReferenceMapping(&referenceMappings[referenceCount++], marker, node);

	// Sort table
	qsort(referenceMappings, referenceCount, sizeof(ReferenceMapping), compareReferenceMappings);

	// Print table
	for (index = 0; index < referenceCount; index++)
		velvetFprintfReferenceMapping(outfile, &referenceMappings[index], refCoords, wordLength);

	// Clean table
	free(referenceMappings);
}

void exportLongNodeMappings(char *filename, Graph * graph, ReadSet * reads,
			     Coordinate minLength, char * sequencesFilename)
{
	FILE * outfile;
	IDnum nodeIndex, refIndex;
	Node *node;
	ReferenceCoord * refCoords;
	IDnum referenceCount = getReferenceCount(reads); 

	if (referenceCount == 0)	
		return;

	refCoords = collectReferenceCoords(sequencesFilename, referenceCount);

	outfile = fopen(filename, "w");
	if (outfile == NULL) {
		velvetLog("Could not write into %s, sorry\n", filename);
		return;
	} else {
		velvetLog("Writing contigs into %s...\n", filename);
	}

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex);

		if (node == NULL || getNodeLength(node) < minLength)
			continue;
		
		exportLongNodeMapping(outfile, node, reads, refCoords, getWordLength(graph));
	}

	for (refIndex = 0; refIndex < referenceCount; refIndex++)
		free(refCoords[refIndex].name);
	free(refCoords);
	fclose(outfile);
}

static void removeLowArcsFromNode(Node* node, Graph * graph, double cutoff) {
	Arc * arc, * next;
	if (node == NULL)
		return;

	for (arc = getArc(node); arc; arc = next) {
		next = getNextArc(arc);
		if (getMultiplicity(arc) <= cutoff)
			destroyArc(arc, graph);
	}
}

void removeLowArcs(Graph * graph, double cutoff) {
	velvetLog("Removing single arcs\n");
	IDnum index;
	for (index = -nodeCount(graph); index <= nodeCount(graph); index++)
		removeLowArcsFromNode(getNodeInGraph(graph, index), graph, cutoff);

	concatenateGraph(graph);
}
