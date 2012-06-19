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

#include "globals.h"
#include "graph.h"
#include "tightString.h"
#include "dfibHeap.h"
#include "recycleBin.h"
#include "passageMarker.h"
#include "shortReadPairs.h"
#include "locallyCorrectedGraph.h"
#include "utility.h"

static const Time INDEL = 0;
static const Time SIM[4][4] = {
	{1, 0, 0, 0},
	{0, 1, 0, 0},
	{0, 0, 1, 0},
	{0, 0, 0, 1}
};

//Global variables used throughout this procedure(internal use only !)
static int MAXREADLENGTH = 100;
static int MAXNODELENGTH = 200;
static double MAXDIVERGENCE = 0.2;
static int MAXGAPS = 3;

static Time *times;
static Node **previous;

static DFibHeapNode **dheapNodes;
static DFibHeap *dheap;

static TightString *fastSequence;
static TightString *slowSequence;

static int WORDLENGTH;
static int SELF_LOOP_CUTOFF;
static Graph *graph;
static Node *start;

static PassageMarker *fastPath;
static PassageMarker *slowPath;

static double **Fmatrix;
//End of global variables;

static void setNodeTime(Node * node, Time time)
{
	times[getNodeID(node) + nodeCount(graph)] = time;
}

static Time getNodeTime(Node * node)
{
	return times[getNodeID(node) + nodeCount(graph)];
}

static Node *getNodePrevious(Node * node)
{
	return previous[getNodeID(node) + nodeCount(graph)];
}

static boolean isPreviousToNode(Node * previous, Node * target)
{
	Node *currentNode = target;
	Node *previousNode = NULL;
	Time targetTime = getNodeTime(target);

	//printf("Testing if %li is previous to %li\n", getNodeID(previous), getNodeID(target));

	while (true) {
		//printf("CCC %li %f\n", getNodeID(currentNode), getNodeTime(currentNode));

		if (currentNode == previous)
			return true;

		if (currentNode == previousNode)
			return false;

		if (getNodeID(currentNode) > nodeCount(graph)
		    || getNodeID(currentNode) < -nodeCount(graph)) {
			printf("Node ID??? %d %d\n",
			       getNodeID(currentNode),
			       getNodeID(previousNode));
		}

		if (getNodeTime(currentNode) != targetTime)
			return false;

		previousNode = currentNode;
		currentNode = getNodePrevious(currentNode);
	}
}

static boolean
extractSequence(PassageMarker * path, TightString * sequence)
{
	PassageMarker *marker;
	Coordinate seqLength = 0;
	Coordinate writeIndex = 0;

	//printf("Extracting sequence %li ... ", pathLength);

	//Measure length
	for (marker = getNextInSequence(path); !isTerminal(marker);
	     marker = getNextInSequence(marker))
		seqLength += getNodeLength(getNode(marker));

	if (seqLength > MAXREADLENGTH)
		return false;
	else
		setTightStringLength(sequence, seqLength);

	//Copy sequences
	for (marker = getNextInSequence(path); !isTerminal(marker);
	     marker = getNextInSequence(marker)) {
		appendNodeSequence(getNode(marker), sequence, writeIndex);
		writeIndex += getNodeLength(getNode(marker));
	}

	return true;
}

static Time max(Time A, Time B, Time C)
{
	if (A >= B && A >= C)
		return A;
	else if (B >= C)
		return B;
	else
		return C;
}

static boolean
compareSequences(TightString * sequence1, TightString * sequence2)
{
	Coordinate i, j;
	Coordinate length1 = getLength(sequence1);
	Coordinate length2 = getLength(sequence2);
	Coordinate maxLength;
	Time Choice1, Choice2, Choice3;
	Time maxScore;

	if (length1 == 0 || length2 == 0)
		return false;

	maxLength = (length1 > length2 ? length1 : length2);

	if (length1 < WORDLENGTH || length2 < WORDLENGTH)
		if (maxLength - length1 > MAXGAPS
		    || maxLength - length2 > MAXGAPS)
			return false;

	for (i = 0; i <= length1; i++)
		Fmatrix[i][0] = 0;
	for (j = 0; j <= length2; j++)
		Fmatrix[0][j] = 0;

	for (i = 1; i <= length1; i++) {
		for (j = 1; j <= length2; j++) {
			Choice1 =
			    Fmatrix[i - 1][j - 1] +
			    SIM[(int) getNucleotide(i - 1, sequence1)]
			    [(int) getNucleotide(j - 1, sequence2)];
			Choice2 = Fmatrix[i - 1][j] + INDEL;
			Choice3 = Fmatrix[i][j - 1] + INDEL;
			Fmatrix[i][j] = max(Choice1, Choice2, Choice3);
		}
	}

	maxScore = Fmatrix[length1][length2];

	if (maxScore < maxLength - MAXGAPS)
		return false;

	if ((1 - maxScore / maxLength) > MAXDIVERGENCE)
		return false;

	return true;
}

static void destroyPaths()
{
	PassageMarker *marker;

	while (slowPath != NULL) {
		marker = slowPath;
		getNodeTime(getNode(marker));
		getNodeTime(getTwinNode(getNode(marker)));

		slowPath = getNextInSequence(marker);
		destroyPassageMarker(marker);
	}

	while (fastPath != NULL) {
		marker = fastPath;
		getNodeTime(getNode(marker));
		getNodeTime(getTwinNode(getNode(marker)));
		fastPath = getNextInSequence(marker);
		destroyPassageMarker(marker);
	}
}

static void cleanUpRedundancy_local()
{
	PassageMarker *current;

	for (current = getNextInSequence(slowPath); !isTerminal(current);
	     current = getNextInSequence(current))
		handicapNode(getNode(current));

	destroyPaths();
}

static void comparePaths_local(Node * destination, Node * origin)
{
	IDnum slowLength, fastLength;
	Node *fastNode, *slowNode;
	IDnum i;
	PassageMarker *marker;

	//Measure lengths
	slowLength = fastLength = 0;
	fastNode = destination;
	slowNode = origin;

	//puts("Looking into separate paths");

	while (fastNode != slowNode) {
		//printf("Fast node %li Slow node %li\n", getNodeID(fastNode), getNodeID(slowNode));

		if (getNodeTime(fastNode) > getNodeTime(slowNode)) {
			fastLength++;
			fastNode = getNodePrevious(fastNode);
		} else if (getNodeTime(fastNode) < getNodeTime(slowNode)) {
			slowLength++;
			slowNode = getNodePrevious(slowNode);
		} else if (isPreviousToNode(slowNode, fastNode)) {
			while (fastNode != slowNode) {
				fastLength++;
				fastNode = getNodePrevious(fastNode);
			}
		} else if (isPreviousToNode(fastNode, slowNode)) {
			while (slowNode != fastNode) {
				slowLength++;
				slowNode = getNodePrevious(slowNode);
			}
		} else {
			fastLength++;
			fastNode = getNodePrevious(fastNode);
			slowLength++;
			slowNode = getNodePrevious(slowNode);
		}

		if (slowLength > MAXNODELENGTH
		    || fastLength > MAXNODELENGTH) {
			//printf("Paths too fragmented %li %li\n", slowLength, fastLength);
			return;
		}
	}

	if (fastLength == 0)
		return;

	//Backtracking to record actual paths
	fastPath = addUncertainPassageMarker(1, destination);
	setPassageMarkerStatus(fastPath, true);

	for (i = 0; i < fastLength; i++) {
		marker =
		    addUncertainPassageMarker(1,
					      getNodePrevious(getNode
							      (fastPath)));
		setPassageMarkerStatus(marker, true);
		connectPassageMarkers(marker, fastPath, graph);
		fastPath = marker;
	}

	slowPath = addUncertainPassageMarker(2, destination);
	setPassageMarkerStatus(slowPath, true);

	marker = addUncertainPassageMarker(2, origin);
	setPassageMarkerStatus(marker, true);
	connectPassageMarkers(marker, slowPath, graph);
	slowPath = marker;

	for (i = 0; i < slowLength; i++) {
		marker =
		    addUncertainPassageMarker(2,
					      getNodePrevious(getNode
							      (slowPath)));
		setPassageMarkerStatus(marker, true);
		connectPassageMarkers(marker, slowPath, graph);
		slowPath = marker;
	}

	//Extract sequences
	if (!extractSequence(fastPath, fastSequence)
	    || !extractSequence(slowPath, slowSequence)) {
		//puts("Paths too long");
		destroyPaths();
		return;
	}
	//Compare sequences
	if (compareSequences(fastSequence, slowSequence)) {
		//puts("Correcting discrepancy");
		cleanUpRedundancy_local();
		return;
	}
	//puts("\tFinished comparing paths, changes made");
	destroyPaths();
}

static void tourBusArc_local(Node * origin, Arc * arc, Time originTime)
{
	Node *destination = getDestination(arc);
	Time arcTime, totalTime, destinationTime;
	IDnum nodeIndex = getNodeID(destination) + nodeCount(graph);
	Node *oldPrevious = previous[nodeIndex];

	//printf("Trying arc from %li -> %li\n", getNodeID(origin), getNodeID(destination)); 

	if (oldPrevious == origin)
		return;

	arcTime =
	    ((Time) getNodeLength(origin)) / ((Time) getMultiplicity(arc));
	totalTime = originTime + arcTime;

	destinationTime = times[nodeIndex];

	if (destinationTime == -1) {
		//puts("New destination");
		setNodeTime(destination, totalTime);
		dheapNodes[nodeIndex] =
		    insertNodeIntoDHeap(dheap, totalTime, destination);
		previous[nodeIndex] = origin;
		return;
	} else if (destinationTime > totalTime) {
		//printf("Previously visited from slower node %li\n", getNodeID(getNodePrevious(destination))); 
		if (dheapNodes[nodeIndex] == NULL) {
			return;
		}

		setNodeTime(destination, totalTime);
		replaceKeyInDHeap(dheap, dheapNodes[nodeIndex], totalTime);
		previous[nodeIndex] = origin;

		comparePaths_local(destination, oldPrevious);
		return;
	} else {
		//printf("Previously visited by faster node %li\n", getNodeID(getNodePrevious(destination))); 
		comparePaths_local(destination, origin);
	}
}

static void tourBusNode_local(Node * node)
{
	Arc *arc;
	Node *destination;
	Time nodeTime = getNodeTime(node);

	//printf("Node %li %f %i %p\n", getNodeID(node),
	//       times[getNodeID(node) + nodeCount(graph)], simpleArcCount(node),
	//       node);

	for (arc = getArc(node); arc != NULL; arc = getNextArc(arc)) {
		destination = getDestination(arc);

		// Node doesn't belong to the marked node area 
		if (getNodeStatus(getDestination(arc)) != 1)
			continue;

		tourBusArc_local(node, arc, nodeTime);

		if (getNodeStatus(node) != 1)
			break;
	}
}

static boolean isLocalDeadEnd(Node * node)
{
	Arc *arc;

	for (arc = getArc(node); arc != NULL; arc = getNextArc(arc))
		if (getNodeStatus(getDestination(arc)) == 1)
			return false;

	return true;
}

static boolean isLocalTwinDeadEnd(Node * node)
{
	Arc *arc;

	for (arc = getArc(getTwinNode(node)); arc != NULL;
	     arc = getNextArc(arc))
		if (getNodeStatus(getTwinNode(getDestination(arc))) == 1)
			return false;

	return true;
}

static void clipTipsVeryHardLocally()
{
	NodeList *nodeList, *next;
	Node *current, *twin;
	boolean modified = true;

	//puts("Clipping short tips off graph HARD");

	while (modified) {
		modified = false;

		for (nodeList = getMarkedNodeList(); nodeList != NULL;
		     nodeList = next) {
			next = nodeList->next;
			current = nodeList->node;

			if (current == NULL || getNodeStatus(current) != 1)
				continue;

			if (getUniqueness(current))
				continue;

			//printf("Checking node HARD %li %i\n", getNodeID(current), simpleArcCount(current));

			twin = getTwinNode(current);

			if (isLocalDeadEnd(current)
			    || isLocalTwinDeadEnd(current)) {
				//printf("Found tip at node %li\n", getNodeID(current));
				handicapNode(current);
				modified = true;
			}
		}
	}
}

static void tourBus_local(Node * startingPoint)
{
	Node *currentNode = startingPoint;
	IDnum nodeID = getNodeID(startingPoint) + nodeCount(graph);

	//printf("Tour bus from node %li...\n", getNodeID(startingPoint));

	times[nodeID] = 0;
	previous[nodeID] = currentNode;

	while (currentNode != NULL) {
		dheapNodes[getNodeID(currentNode) + nodeCount(graph)] =
		    NULL;
		tourBusNode_local(currentNode);
		currentNode = removeNextNodeFromDHeap(dheap);
	}
}

void prepareGraphForLocalCorrections(Graph * argGraph)
{
	IDnum nodes = nodeCount(argGraph);
	IDnum index;

	//Setting global params
	graph = argGraph;
	WORDLENGTH = getWordLength(graph);;
	SELF_LOOP_CUTOFF = WORDLENGTH;
	// Done with global params

	// Original
	/*
	printf("Preparing to correct graph with cutoff %f\n",
	       MAXDIVERGENCE);
	*/
	// Original

	// Allocating memory
	times = mallocOrExit(2 * nodes + 1, Time);
	previous = mallocOrExit(2 * nodes + 1, Node *);

	dheapNodes = mallocOrExit(2 * nodes + 1, DFibHeapNode *);

	dheap = newDFibHeap();

	fastSequence = newTightString(MAXREADLENGTH);
	slowSequence = newTightString(MAXREADLENGTH);

	for (index = 0; index < (2 * nodeCount(graph) + 1); index++) {
		times[index] = -1;
		dheapNodes[index] = NULL;
		previous[index] = NULL;
	}

	Fmatrix = callocOrExit(MAXREADLENGTH + 1, double *);
	for (index = 0; index < MAXREADLENGTH + 1; index++)
		Fmatrix[index] = callocOrExit(MAXREADLENGTH + 1, double);
	//Done with memory 
}

void correctGraphLocally(Node * argStart)
{
	IDnum index, nodeIndex;
	NodeList *nodeList;

	start = argStart;
	//printf("Correcting graph from node %li\n", getNodeID(start));

	clipTipsVeryHardLocally();

	index = 0;
	for (nodeList = getMarkedNodeList(); nodeList != NULL;
	     nodeList = nodeList->next) {
		nodeIndex = getNodeID(nodeList->node) + nodeCount(graph);
		times[nodeIndex] = -1;
		dheapNodes[nodeIndex] = NULL;
		previous[nodeIndex] = NULL;
	}

	tourBus_local(start);
}

void deactivateLocalCorrectionSettings()
{
	// Original
	/*
	puts("Deactivating local correction settings");
	*/
	// Original

	IDnum index;

	for (index = 0; index <= MAXREADLENGTH; index++) {
		free(Fmatrix[index]);
	}
	free(Fmatrix);

	free(times);
	free(previous);
	free(dheapNodes);
	destroyDHeap(dheap);

	destroyTightString(fastSequence);
	destroyTightString(slowSequence);
}

void setLocalMaxReadLength(int value)
{
	MAXREADLENGTH = value;
	MAXNODELENGTH = 2 * value;
}

void setLocalMaxGaps(int value)
{
	MAXGAPS = value;
}

void setLocalMaxDivergence(double value)
{
	MAXDIVERGENCE = value;
}
