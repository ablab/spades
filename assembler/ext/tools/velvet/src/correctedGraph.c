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
#include <time.h>

#include "globals.h"
#include "graph.h"
#include "tightString.h"
#include "dfibHeap.h"
#include "fibHeap.h"
#include "recycleBin.h"
#include "passageMarker.h"
#include "concatenatedGraph.h"
#include "graphStats.h"
#include "utility.h"

#define TICKET_BLOCK_SIZE 10000

static const Time INDEL = 0;
static const Time SIM[4][4] = {
	{1, 0, 0, 0},
	{0, 1, 0, 0},
	{0, 0, 1, 0},
	{0, 0, 0, 1}
};

typedef struct tkt_st Ticket;

struct tkt_st {
	Ticket *next;
	IDnum id_a;
} ATTRIBUTE_PACKED;

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

static Node *activeNode;
static Node *startingNode;
static int WORDLENGTH;
static Graph *graph;
static IDnum dbgCounter;

static PassageMarkerI fastPath;
static PassageMarkerI slowPath;

static IDnum *eligibleStartingPoints;

static double **Fmatrix;
static Coordinate *slowToFastMapping;
static Coordinate *fastToSlowMapping;

static RecycleBin *ticketMemory;

static Ticket **todoLists;
static Ticket **todo;
static Ticket *done;
static boolean *progressStatus;

static ShortLength *sequenceLengths;
static Category *sequenceCategories;

//End of global variables;

static void setNodeTime(Node * node, Time time)
{
	times[getNodeID(node) + nodeCount(graph)] = time;
}

static Time getNodeTime(Node * node)
{
	return times[getNodeID(node) + nodeCount(graph)];
}

static void setNodePrevious(Node * previousNode, Node * node)
{
	previous[getNodeID(node) + nodeCount(graph)] = previousNode;
}

static Node *getNodePrevious(Node * node)
{
	return previous[getNodeID(node) + nodeCount(graph)];
}

static void setNodeDHeapNode(Node * node, DFibHeapNode * dheapNode)
{
	dheapNodes[getNodeID(node) + nodeCount(graph)] = dheapNode;
}

static DFibHeapNode *getNodeDHeapNode(Node * node)
{
	return dheapNodes[getNodeID(node) + nodeCount(graph)];
}

static Ticket *newTicket()
{
	if (ticketMemory == NULL)
		ticketMemory =
		    newRecycleBin(sizeof(Ticket), TICKET_BLOCK_SIZE);

	return allocatePointer(ticketMemory);
}

static boolean isPreviousToNode(Node * previous, Node * target)
{
	Node *currentNode = target;
	Node *previousNode = NULL;
	Time targetTime = getNodeTime(target);

	//velvetLog("Testing if %ld is previous to %ld\n", getNodeID(previous), getNodeID(target));

	while (true) {
		if (currentNode == previous)
			return true;

		if (currentNode == previousNode)
			return false;

		if (getNodeTime(currentNode) != targetTime)
			return false;

		previousNode = currentNode;
		currentNode = getNodePrevious(currentNode);
	}
}

static void concatenateCommonTodoLists(Node * nodeA, Node * nodeB)
{
	Ticket **listA = &todoLists[getNodeID(nodeA) + nodeCount(graph)];
	Ticket **listB = &todoLists[getNodeID(nodeB) + nodeCount(graph)];
	Ticket *head = NULL;
	Ticket *tail = NULL;
	Ticket *tmp;
	IDnum idA, idB;
	IDnum targetID = getNodeID(nodeA);
	IDnum indexA, indexB;
	IDnum nodes = nodeCount(graph);

	//velvetLog("Merging todo list %ld into %ld\n", getNodeID(nodeB),
	//       getNodeID(nodeA));

	if (*listB == NULL)
		return;

	if (*listA == NULL) {
		*listA = *listB;
		*listB = NULL;
		return;
	}

	while (*listA != NULL && *listB != NULL) {
		idA = (*listA)->id_a;
		idB = (*listB)->id_a;
		indexA = idA + nodes;
		indexB = idB + nodes;

		if (previous[indexA] == nodeA) {
			tmp = *listA;
			*listA = (*listA)->next;
			deallocatePointer(ticketMemory, tmp);
			continue;
		}

		if (idB == targetID || previous[indexB] == nodeA) {
			tmp = *listB;
			*listB = (*listB)->next;
			deallocatePointer(ticketMemory, tmp);
			continue;
		}

		if (idA > idB) {
			tmp = *listB;
			*listB = (*listB)->next;
		} else if (idA < idB) {
			tmp = *listA;
			*listA = (*listA)->next;
		} else {
			tmp = *listB;
			*listB = (*listB)->next;
			deallocatePointer(ticketMemory, tmp);

			tmp = *listA;
			*listA = (*listA)->next;
		}

		if (tail == NULL) {
			tail = tmp;
			head = tail;
		} else {
			tail->next = tmp;
			tail = tail->next;
		}
	}

	while (*listA != NULL) {
		idA = (*listA)->id_a;
		indexA = idA + nodes;

		if (previous[indexA] == nodeA) {
			tmp = *listA;
			*listA = (*listA)->next;
			deallocatePointer(ticketMemory, tmp);
		} else if (tail != NULL) {
			tail->next = *listA;
			*listA = (*listA)->next;
			tail = tail->next;
		} else {
			head = *listA;
			*listA = (*listA)->next;
			tail = head;
		}
	}

	while (*listB != NULL) {
		idB = (*listB)->id_a;
		indexB = idB + nodes;

		if (idB == targetID || previous[indexB] == nodeA) {
			tmp = *listB;
			*listB = (*listB)->next;
			deallocatePointer(ticketMemory, tmp);
		} else if (tail != NULL) {
			tail->next = *listB;
			*listB = (*listB)->next;
			tail = tail->next;
		} else {
			head = *listB;
			*listB = (*listB)->next;
			tail = head;

		}
	}

	if (tail != NULL)
		tail->next = NULL;

	*listA = head;
	*listB = NULL;
}

static void concatenateTodoListIntoActive(Node * nodeB)
{
	Ticket **listB = &todoLists[getNodeID(nodeB) + nodeCount(graph)];
	Ticket *head = NULL;
	Ticket *tail = NULL;
	Ticket *tmp;
	IDnum nodes = nodeCount(graph);
	IDnum idA, idB;
	IDnum activeID = getNodeID(activeNode);
	IDnum indexB, indexA;

	//velvetLog("Merging todo list %ld into active node %ld\n",
	//       getNodeID(nodeB), getNodeID(activeNode));

	if (*listB == NULL)
		return;

	if (*todo == NULL) {
		*todo = *listB;
		*listB = NULL;
		return;
	}

	while (*todo != NULL && *listB != NULL) {
		idA = (*todo)->id_a;
		idB = (*listB)->id_a;
		indexA = idA + nodes;
		indexB = idB + nodes;

		if (previous[indexA] == activeNode
		    || progressStatus[indexA]) {
			tmp = *todo;
			*todo = (*todo)->next;
			deallocatePointer(ticketMemory, tmp);
			continue;
		}

		if (idB == activeID || previous[indexB] == activeNode
		    || progressStatus[indexB]) {
			tmp = *listB;
			*listB = (*listB)->next;
			deallocatePointer(ticketMemory, tmp);
			continue;
		}

		if (idA > idB) {
			tmp = *listB;
			*listB = (*listB)->next;
		} else if (idA < idB) {
			tmp = *todo;
			*todo = (*todo)->next;
		} else {
			tmp = *listB;
			*listB = (*listB)->next;
			deallocatePointer(ticketMemory, tmp);

			tmp = *todo;
			*todo = (*todo)->next;
		}

		if (tail == NULL) {
			tail = tmp;
			head = tail;
		} else {
			tail->next = tmp;
			tail = tmp;
		}
	}

	while (*todo != NULL) {
		idA = (*todo)->id_a;
		indexA = idA + nodes;

		if (previous[indexA] == activeNode
		    || progressStatus[indexA]) {
			tmp = *todo;
			*todo = (*todo)->next;
			deallocatePointer(ticketMemory, tmp);
		} else if (tail != NULL) {
			tail->next = *todo;
			*todo = (*todo)->next;
			tail = tail->next;
		} else {
			head = *todo;
			*todo = (*todo)->next;
			tail = head;
		}
	}

	while (*listB != NULL) {
		idB = (*listB)->id_a;
		indexB = idB + nodes;

		if (idB == activeID || previous[indexB] == activeNode
		    || progressStatus[indexB]) {
			tmp = *listB;
			*listB = (*listB)->next;
			deallocatePointer(ticketMemory, tmp);
		} else if (tail != NULL) {
			tail->next = *listB;
			*listB = (*listB)->next;
			tail = tail->next;
		} else {
			head = *listB;
			*listB = (*listB)->next;
			tail = head;

		}
	}

	if (tail != NULL)
		tail->next = NULL;
	*todo = head;
	*listB = NULL;
}

static void concatenateTodoLists(Node * nodeA, Node * nodeB)
{
	if (nodeA == activeNode)
		concatenateTodoListIntoActive(nodeB);
	else
		concatenateCommonTodoLists(nodeA, nodeB);
}

static IDnum nextTodoTicket()
{
	Ticket *tkt;
	IDnum index;

	while (*todo != NULL) {
		tkt = *todo;
		*todo = tkt->next;

		index = tkt->id_a + nodeCount(graph);

		if (previous[index] == activeNode) {
			deallocatePointer(ticketMemory, tkt);
			continue;
		}

		progressStatus[index] = true;

		tkt->next = done;
		done = tkt;

		return tkt->id_a;
	}

	return 0;
}

static void freeDoneTickets()
{
	Ticket *tkt;
	IDnum nodes = nodeCount(graph);

	while (done != NULL) {
		tkt = done;
		done = tkt->next;
		progressStatus[tkt->id_a + nodes] = false;
		deallocatePointer(ticketMemory, tkt);
	}
}

static void updateNodeStatus(Node * node)
{
	FibHeap *heap = newFibHeap();
	Arc *arc;
	Node *currentNode = node;
	Node *destination;

	setNodeStatus(currentNode, true);

	while (currentNode != NULL) {
		for (arc = getArc(currentNode); arc != NULL;
		     arc = getNextArc(arc)) {
			destination = getDestination(arc);
			if (getNodeStatus(destination) > 1) {
				setNodeStatus(destination, true);
				insertNodeIntoHeap(heap,
						   getNodeID(destination),
						   destination);
			}
		}

		currentNode = removeNextNodeFromHeap(heap);
	}

	destroyHeap(heap);
}

static void determineEligibleStartingPoints()
{
	IDnum nodeIndex;
	IDnum maxmult;
	Node *node;
	Arc *arc;
	IDnum counter = 0;
	FibHeap *heap = newFibHeap();

	velvetLog("Determining eligible starting points\n");

	for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex);
		if (node == NULL)
			continue;

		maxmult = 0;
		for (arc = getArc(node); arc != NULL;
		     arc = getNextArc(arc))
			if (getMultiplicity(arc) > maxmult)
				maxmult = getMultiplicity(arc);

		insertNodeIntoHeap(heap, -maxmult, node);

		// Same for twin node
		node = getNodeInGraph(graph, -nodeIndex);
		maxmult = 0;
		for (arc = getArc(node); arc != NULL;
		     arc = getNextArc(arc))
			if (getMultiplicity(arc) > maxmult)
				maxmult = getMultiplicity(arc);

		insertNodeIntoHeap(heap, -maxmult, node);
	}

	while ((node = removeNextNodeFromHeap(heap)) != NULL)
		eligibleStartingPoints[counter++] = getNodeID(node);

	destroyHeap(heap);
	velvetLog("Done listing starting nodes\n");
}

static Node *nextStartingPoint()
{
	static IDnum counter = 0;
	Node *result = NULL;

	while (result == NULL || getNodeStatus(result) > 0) {
		if (counter >= nodeCount(graph) * 2)
			return NULL;

		result =
		    getNodeInGraph(graph,
				   eligibleStartingPoints[counter++]);
	}

	return result;
}

static boolean
extractSequence(PassageMarkerI path, TightString * sequence)
{
	PassageMarkerI marker;
	Coordinate seqLength = 0;
	Coordinate writeIndex = 0;

	//velvetLog("Extracting sequence %ld ... ", pathLength);

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

	if (length1 < WORDLENGTH || length2 < WORDLENGTH) {
		if (maxLength - length1 > MAXGAPS
		    || maxLength - length2 > MAXGAPS)
			return false;
		if (WORDLENGTH - length1 > MAXGAPS
		    || WORDLENGTH - length2 > MAXGAPS)
			return false;
	}

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

static void mapSlowOntoFast()
{
	Coordinate slowIndex = getLength(slowSequence);
	Coordinate fastIndex = getLength(fastSequence);
	int fastn, slown;

	if (slowIndex == 0) {
		slowToFastMapping[0] = fastIndex;

		while (fastIndex >= 0)
			fastToSlowMapping[fastIndex--] = 0;

		return;
	}

	if (fastIndex == 0) {
		while (slowIndex >= 0)
			slowToFastMapping[slowIndex--] = 0;

		fastToSlowMapping[0] = slowIndex;

		return;
	}

	while (slowIndex > 0 && fastIndex > 0) {
		fastn = (int) getNucleotide(fastIndex - 1, fastSequence);
		slown = (int) getNucleotide(slowIndex - 1, slowSequence);

		if (Fmatrix[fastIndex][slowIndex] ==
		    Fmatrix[fastIndex - 1][slowIndex - 1] +
		    SIM[fastn][slown]) {
			fastToSlowMapping[--fastIndex] = --slowIndex;
			slowToFastMapping[slowIndex] = fastIndex;
		} else if (Fmatrix[fastIndex][slowIndex] ==
			   Fmatrix[fastIndex - 1][slowIndex] + INDEL)
			fastToSlowMapping[--fastIndex] = slowIndex - 1;

		else if (Fmatrix[fastIndex][slowIndex] ==
			 Fmatrix[fastIndex][slowIndex - 1] + INDEL)
			slowToFastMapping[--slowIndex] = fastIndex - 1;

		else {
			velvetLog("Error\n");
			fflush(stdout);
			abort();
		}
	}

	while (slowIndex > 0)
		slowToFastMapping[--slowIndex] = -1;
	while (fastIndex > 0)
		fastToSlowMapping[--fastIndex] = -1;

	slowToFastMapping[getLength(slowSequence)] =
	    getLength(fastSequence);
	fastToSlowMapping[getLength(fastSequence)] =
	    getLength(slowSequence);
}

static void remapNodeArcsOntoTarget(Node * source, Node * target)
{
	Arc *arc;

	if (source == activeNode) {
		activeNode = target;
		todo =
		    &todoLists[getNodeID(activeNode) + nodeCount(graph)];
	}
	concatenateTodoLists(target, source);

	arc = getArc(source);
	while (arc != NULL) {
		createAnalogousArc(target, getDestination(arc), arc, graph);
		destroyArc(arc, graph);
		arc = getArc(source);
	}
}

static void remapNodeArcsOntoNeighbour(Node * source, Node * target)
{
	remapNodeArcsOntoTarget(source, target);
	remapNodeArcsOntoTarget(getTwinNode(source), getTwinNode(target));
}

static void remapNodeMarkersOntoNeighbour(Node * source,
					  PassageMarkerI sourceMarker,
					  Node * target,
					  PassageMarkerI targetMarker)
{
	PassageMarkerI marker;
	Coordinate offset;
	IDnum sourceLength, index;
	ShortReadMarker *sourceArray, *shortMarker;
	Coordinate position;

	Coordinate targetStart = getPassageMarkerStart(targetMarker);
	Coordinate targetFinish = getPassageMarkerFinish(targetMarker);
	Coordinate sourceStart = getPassageMarkerStart(sourceMarker);
	Coordinate sourceFinish = getPassageMarkerFinish(sourceMarker);

	Coordinate alignedTargetLength = targetFinish - targetStart;
	Coordinate alignedSourceLength = sourceFinish - sourceStart;

	Coordinate realTargetLength = getNodeLength(target);
	Coordinate realSourceLength = getNodeLength(source);

	while (getMarker(source) != NULL_IDX) {
		marker = getMarker(source);
		extractPassageMarker(marker);
		transposePassageMarker(marker, target);

		if (realSourceLength != 0 && alignedTargetLength != 0) {
			if (isInitial(marker)) {
				offset = getStartOffset(marker);
				offset *= alignedSourceLength;
				offset /= realSourceLength;
				offset += sourceStart;
				offset = slowToFastMapping[offset];
				offset -= targetStart;
				offset *= realTargetLength;
				offset /= alignedTargetLength;

				if (offset < 0)
					offset = 0;
				if (offset > realTargetLength)
					offset = realTargetLength;
			} else
				offset = 0;

			setStartOffset(marker, offset);

			if (isTerminal(marker)) {
				offset = getFinishOffset(marker);
				offset *= alignedSourceLength;
				offset /= realSourceLength;
				offset = sourceFinish - offset;
				offset = slowToFastMapping[offset];
				offset = targetFinish - offset;
				offset *= realTargetLength;
				offset /= alignedTargetLength;

				if (offset < 0)
					offset = 0;
				if (offset > realTargetLength)
					offset = realTargetLength;
			} else
				offset = 0;

			setFinishOffset(marker, offset);
		} else {
			setStartOffset(marker, 0);
			setFinishOffset(marker, 0);
		}
	}

	// Short read markers 
	if (readStartsAreActivated(graph)) {
		// Update Coordinates
		sourceArray = getNodeReads(source, graph);
		sourceLength = getNodeReadCount(source, graph);

		for (index = 0; index < sourceLength; index++) {
			shortMarker =
			    getShortReadMarkerAtIndex(sourceArray, index);
			position = getShortReadMarkerPosition(shortMarker);

			if (position > -1) {
				if (realSourceLength > 0
				    && alignedTargetLength > 0) {
					position *= alignedSourceLength;
					position /= realSourceLength;
					position += sourceStart;
					position =
					    slowToFastMapping[position];
					position -= targetStart;
					position *= realTargetLength;
					position /= alignedTargetLength;

					if (position < 0)
						position = 0;
					if (position > realTargetLength)
						position =
						    realTargetLength;
				} else
					position = 0;
			}

			setShortReadMarkerPosition(shortMarker, position);
		}
		mergeNodeReads(target, source, graph);

		// Same but for symmetrical reads
		sourceArray = getNodeReads(getTwinNode(source), graph);
		sourceLength =
		    getNodeReadCount(getTwinNode(source), graph);

		for (index = 0; index < sourceLength; index++) {
			shortMarker =
			    getShortReadMarkerAtIndex(sourceArray, index);
			position = getShortReadMarkerPosition(shortMarker);

			if (position > -1) {
				if (realSourceLength > 0
				    && alignedTargetLength > 0) {
					position =
					    realSourceLength - position;
					position *= alignedSourceLength;
					position /= realSourceLength;
					position += sourceStart;
					position =
					    slowToFastMapping[position];
					position -= targetStart;
					position *= realTargetLength;
					position /= alignedTargetLength;
					position =
					    realTargetLength - position;

					if (position < 0)
						position = 0;
					if (position > realTargetLength)
						position =
						    realTargetLength;
				} else
					position = 0;
			}

			setShortReadMarkerPosition(shortMarker, position);
		}
		mergeNodeReads(getTwinNode(target), getTwinNode(source),
			       graph);
	}
	// Virtual reads
#ifndef SINGLE_COV_CAT
	Category cat;
	for (cat = 0; cat < CATEGORIES; cat++)
		incrementVirtualCoverage(target, cat,
					 getVirtualCoverage(source, cat));
#else
	incrementVirtualCoverage(target, getVirtualCoverage(source));
#endif
}

static void remapBackOfNodeArcsOntoNeighbour(Node * source, Node * target)
{
	Arc *arc;

	remapNodeArcsOntoTarget(getTwinNode(source), getTwinNode(target));
	for (arc = getArc(source); arc != NULL; arc = getNextArc(arc))
		createAnalogousArc(target, source, arc, graph);

}

static Coordinate
remapBackOfNodeMarkersOntoNeighbour(Node * source,
				    PassageMarkerI sourceMarker,
				    Node * target,
				    PassageMarkerI targetMarker,
				    boolean slowToFast)
{
	PassageMarkerI marker, newMarker, previousMarker, nextMarker;
	Coordinate halfwayPoint, halfwayPointOffset, breakpoint,
	    newStartOffset, newFinishOffset;
	Coordinate *targetToSourceMapping, *sourceToTargetMapping;
	ShortReadMarker *selectedShortReads, *shortRead;
	IDnum selectedShortReadCount, shortReadIndex;
	Coordinate position;

	Coordinate targetStart = getPassageMarkerStart(targetMarker);
	Coordinate targetFinish = getPassageMarkerFinish(targetMarker);
	Coordinate sourceStart = getPassageMarkerStart(sourceMarker);
	Coordinate sourceFinish = getPassageMarkerFinish(sourceMarker);

	Coordinate alignedTargetLength = targetFinish - targetStart;
	Coordinate alignedSourceLength = sourceFinish - sourceStart;

	Coordinate realTargetLength = getNodeLength(target);
	Coordinate realSourceLength = getNodeLength(source);

	if (slowToFast) {
		sourceToTargetMapping = slowToFastMapping;
		targetToSourceMapping = fastToSlowMapping;
	} else {
		sourceToTargetMapping = fastToSlowMapping;
		targetToSourceMapping = slowToFastMapping;
	}

	// Calculating source node breakpoint:
	if (alignedSourceLength > 0 && targetFinish > 0) {
		halfwayPoint =
		    targetToSourceMapping[targetFinish - 1] - sourceStart +
		    1;
		halfwayPoint *= realSourceLength;
		halfwayPoint /= alignedSourceLength;
	} else
		halfwayPoint = 0;

	if (halfwayPoint < 0)
		halfwayPoint = 0;
	if (halfwayPoint > realSourceLength)
		halfwayPoint = realSourceLength;
	halfwayPointOffset = realSourceLength - halfwayPoint;

	// Complete markers
	for (marker = getMarker(source); marker != NULL_IDX;
	     marker = nextMarker) {
		nextMarker = getNextInNode(marker);

		// To avoid making loops...
		if (getNode(getPreviousInSequence(marker)) == target)
			continue;

		// Markers which are downstream of the breakpoint
		if (getStartOffset(marker) > halfwayPoint) {
			newStartOffset =
			    getStartOffset(marker) - halfwayPoint;
			setStartOffset(marker, newStartOffset);
			continue;
		}
		// Markers which are upstream of the breakpoint
		if (getFinishOffset(marker) > halfwayPointOffset) {
			if (slowToFast) {
				if (realSourceLength > 0
				    && alignedTargetLength > 0) {
					newFinishOffset =
					    getFinishOffset(marker) -
					    halfwayPointOffset;
					newFinishOffset *=
					    alignedSourceLength;
					newFinishOffset /=
					    realSourceLength;
					newFinishOffset *=
					    realTargetLength;
					newFinishOffset /=
					    alignedTargetLength;
					if (newFinishOffset < 0)
						newFinishOffset = 0;
					else if (newFinishOffset >
						 realTargetLength)
						newFinishOffset =
						    realTargetLength;
				} else
					newFinishOffset = 0;
			} else {
				newFinishOffset =
				    getFinishOffset(marker) -
				    halfwayPointOffset;
			}
			setFinishOffset(marker, newFinishOffset);
			extractPassageMarker(marker);
			transposePassageMarker(marker, target);
			continue;
		}
		// Markers on both sides of the divide
		newMarker =
		    addPassageMarker(getPassageMarkerSequenceID(marker),
				     getPassageMarkerStart(marker),
				     target);

		setPassageMarkerStart(newMarker,
				      getPassageMarkerStart(marker));
		setPassageMarkerStatus(newMarker,
				       getPassageMarkerStatus(marker));

		if (realSourceLength - getStartOffset(marker) -
		    getFinishOffset(marker) > 0) {
			breakpoint = halfwayPoint - getStartOffset(marker);
			breakpoint *= getPassageMarkerLength(marker);
			breakpoint /= realSourceLength -
			    getStartOffset(marker) -
			    getFinishOffset(marker);
			if (breakpoint > getPassageMarkerLength(marker))
				breakpoint = getPassageMarkerLength(marker);
			breakpoint *= passageMarkerDirection(marker);
			breakpoint += getPassageMarkerStart(marker);
		} else {
			breakpoint = getPassageMarkerStart(marker);
		}

		setPassageMarkerFinish(newMarker, breakpoint);
		setPassageMarkerStart(marker, breakpoint);

		if (slowToFast) {
			if (realSourceLength != 0
			    && alignedTargetLength != 0) {
				newStartOffset = getStartOffset(marker);
				newStartOffset *= alignedSourceLength;
				newStartOffset /= realSourceLength;
				newStartOffset *= realTargetLength;
				newStartOffset /= alignedTargetLength;
				if (newStartOffset < 0)
					newStartOffset = 0;
				else if (newStartOffset > realTargetLength)
					newStartOffset = realTargetLength;
			} else {
				newStartOffset = 0;
			}
		} else {
			newStartOffset = getStartOffset(marker);
		}

		setStartOffset(newMarker, newStartOffset);
		setFinishOffset(newMarker, 0);
		setStartOffset(marker, 0);

		previousMarker = getPreviousInSequence(marker);
		setNextInSequence(previousMarker, newMarker);
		setPreviousInSequence(previousMarker, newMarker);

		setPreviousInSequence(newMarker, marker);
		setNextInSequence(newMarker, marker);
	}

	// Read starts
	if (readStartsAreActivated(graph)) {
		selectedShortReads =
		    extractBackOfNodeReads(source, halfwayPoint, graph,
					   &selectedShortReadCount,
					   sourceMarker, sequenceLengths);
		if (slowToFast) {
			if (realSourceLength > 0
			    && alignedTargetLength > 0) {
				for (shortReadIndex = 0;
				     shortReadIndex <
				     selectedShortReadCount;
				     shortReadIndex++) {
					shortRead =
					    getShortReadMarkerAtIndex
					    (selectedShortReads,
					     shortReadIndex);
					position =
					    getShortReadMarkerPosition
					    (shortRead);
					if (position > -1) {
						position *=
						    alignedSourceLength;
						position /=
						    realSourceLength;
						position += sourceStart;
						position =
						    sourceToTargetMapping
						    [position];
						position -= targetStart;
						position *=
						    realTargetLength;
						position /=
						    alignedTargetLength;
						if (position < 0)
							position = 0;
						if (position >
						    realTargetLength)
							position =
							    realTargetLength;
					}
					setShortReadMarkerPosition
					    (shortRead, position);
				}
			} else {
				for (shortReadIndex = 0;
				     shortReadIndex <
				     selectedShortReadCount;
				     shortReadIndex++) {
					shortRead =
					    getShortReadMarkerAtIndex
					    (selectedShortReads,
					     shortReadIndex);
					position =
					    getShortReadMarkerPosition
					    (shortRead);
					if (position > -1)
						setShortReadMarkerPosition
						    (shortRead, 0);
				}

			}
		}
		injectShortReads(selectedShortReads,
				 selectedShortReadCount, target, graph);

		selectedShortReads =
		    extractFrontOfNodeReads(getTwinNode(source),
					    halfwayPoint, graph,
					    &selectedShortReadCount,
					    sourceMarker, sequenceLengths);
		if (slowToFast) {
			if (realSourceLength > 0
			    && alignedTargetLength > 0) {
				for (shortReadIndex = 0;
				     shortReadIndex <
				     selectedShortReadCount;
				     shortReadIndex++) {
					shortRead =
					    getShortReadMarkerAtIndex
					    (selectedShortReads,
					     shortReadIndex);
					position =
					    getShortReadMarkerPosition
					    (shortRead);
					if (position > -1) {
						position =
						    getShortReadMarkerPosition
						    (shortRead);
						position =
						    realSourceLength -
						    position;
						position *=
						    alignedSourceLength;
						position /=
						    realSourceLength;
						position += sourceStart;
						position =
						    sourceToTargetMapping
						    [position];
						position -= targetStart;
						position *=
						    realTargetLength;
						position /=
						    alignedTargetLength;
						position =
						    realTargetLength -
						    position;
						if (position < 0)
							position = 0;
						if (position >
						    realTargetLength)
							position =
							    realTargetLength;
					}
					setShortReadMarkerPosition
					    (shortRead, position);
				}
			} else {
				for (shortReadIndex = 0;
				     shortReadIndex <
				     selectedShortReadCount;
				     shortReadIndex++) {
					shortRead =
					    getShortReadMarkerAtIndex
					    (selectedShortReads,
					     shortReadIndex);
					position =
					    getShortReadMarkerPosition
					    (shortRead);
					if (position > -1)
						setShortReadMarkerPosition
						    (shortRead, 0);
				}

			}
		}
		injectShortReads(selectedShortReads,
				 selectedShortReadCount,
				 getTwinNode(target), graph);
	}
	// Virtual coverage
	if (alignedSourceLength != 0) {
		Coordinate coverage;
#ifndef SINGLE_COV_CAT
		Category cat;
		for (cat = 0; cat < CATEGORIES; cat++) {
			coverage = getVirtualCoverage(source, cat);
			coverage *= halfwayPoint;
			coverage /= alignedSourceLength;
			incrementVirtualCoverage(target, cat, coverage);
			incrementVirtualCoverage(source, cat, -coverage);

			coverage = getOriginalVirtualCoverage(source, cat);
			coverage *= halfwayPoint;
			coverage /= alignedSourceLength;
			incrementOriginalVirtualCoverage(source, cat,
							 -coverage);
		}
#else
		coverage = getVirtualCoverage(source);
		coverage *= halfwayPoint;
		coverage /= alignedSourceLength;
		incrementVirtualCoverage(target, coverage);
		incrementVirtualCoverage(source, -coverage);
#endif
	}

	return halfwayPointOffset;
}

static void remapNodeInwardReferencesOntoNode(Node * source, Node * target)
{
	Arc *arc;
	Node *destination;

	for (arc = getArc(source); arc != NULL; arc = getNextArc(arc)) {
		destination = getDestination(arc);
		if (destination == target || destination == source)
			continue;
		if (getNodePrevious(destination) == source)
			setNodePrevious(target, destination);
	}
}

static void remapNodeTimesOntoTargetNode(Node * source, Node * target)
{
	Time nodeTime = getNodeTime(source);
	Node *previous = getNodePrevious(source);
	Time targetTime = getNodeTime(target);

	if (nodeTime == -1)
		return;

	if (previous == source) {
		setNodeTime(target, nodeTime);
		setNodePrevious(target, target);
	} else if (targetTime == -1
		   || targetTime > nodeTime
		   || (targetTime == nodeTime
		       && !isPreviousToNode(target, previous))) {
		setNodeTime(target, nodeTime);
		if (previous != getTwinNode(source))
			setNodePrevious(previous, target);
		else
			setNodePrevious(getTwinNode(target), target);
	}

	remapNodeInwardReferencesOntoNode(source, target);

	setNodePrevious(NULL, source);
}

static void foldSymmetricalNode(Node * node)
{
	Node *twinNode = getTwinNode(node);
	Node *tmp, *destination;
	Arc *arc;
	PassageMarkerI oldMarker = getMarker(node);
	PassageMarkerI currentMarker, newMarker, previousMarker;
	Coordinate halfwayPoint;
	IDnum totalMult;

	// Reduce time complexity of damn thing
	if (simpleArcCount(node) < simpleArcCount(twinNode)) {
		tmp = twinNode;
		twinNode = node;
		node = tmp;
	}
	// Destroy link to old markers 
	setMarker(node, NULL_IDX);

	// Reinsert markers properly
	while (oldMarker != NULL_IDX) {
		currentMarker = oldMarker;
		oldMarker = getNextInNode(currentMarker);
		previousMarker = getPreviousInSequence(currentMarker);

		if (getNode(previousMarker) != twinNode) {
			newMarker =
			    addUncertainPassageMarker
			    (getPassageMarkerSequenceID(currentMarker),
			     twinNode);
			setPassageMarkerStatus(newMarker,
					       getPassageMarkerStatus
					       (currentMarker));

			setPassageMarkerStart(newMarker,
					      getPassageMarkerStart
					      (currentMarker));

			// For security issues:
			if (currentMarker == slowPath)
				slowPath = newMarker;
			else if (currentMarker == fastPath)
				fastPath = newMarker;

			halfwayPoint =
			    (getPassageMarkerStart(currentMarker) +
			     getPassageMarkerFinish(currentMarker))
			    / 2;
			setPassageMarkerFinish(newMarker, halfwayPoint);

			setPassageMarkerStart(currentMarker, halfwayPoint);

			setStartOffset(newMarker,
				       getStartOffset(currentMarker));
			setFinishOffset(newMarker, 0);
			setStartOffset(currentMarker, 0);

			setPreviousInSequence(previousMarker, newMarker);
			setNextInSequence(previousMarker, newMarker);

			setPreviousInSequence(newMarker, currentMarker);
			setNextInSequence(newMarker, currentMarker);
		}

		transposePassageMarker(currentMarker, node);
	}

	// Read start info
	foldSymmetricalNodeReads(node, graph);

	// Coverage => already balanced!

	// References
	if (getNodeTime(node) == -1 && getNodeTime(twinNode) == -1);
	else if (getNodeTime(node) == -1) {
		setNodeTime(node, getNodeTime(twinNode));
	} else if (getNodeTime(twinNode) == -1) {
		setNodeTime(twinNode, getNodeTime(node));
		setNodePrevious(getNodePrevious(node), twinNode);
	} else if (getNodePrevious(node) == node) {
		setNodeTime(twinNode, getNodeTime(node));
		setNodePrevious(twinNode, twinNode);
	} else if (getNodeTime(node) < getNodeTime(twinNode)) {
		setNodeTime(twinNode, getNodeTime(node));
		setNodePrevious(getNodePrevious(node), twinNode);
	} else if (getNodeTime(node) == getNodeTime(twinNode)
		   && isPreviousToNode(node, twinNode)) {
		setNodePrevious(getNodePrevious(node), twinNode);
	} else {
		setNodeTime(node, getNodeTime(twinNode));
	}

	setNodePrevious(twinNode, node);
	remapNodeInwardReferencesOntoNode(twinNode, node);

	// Active node
	if (twinNode == activeNode) {
		activeNode = node;
		todo =
		    &todoLists[getNodeID(activeNode) + nodeCount(graph)];
	}
	concatenateTodoLists(node, twinNode);

	// Remap arcs properly
	arc = getArc(twinNode);
	totalMult = 0;
	while (arc != NULL) {
		destination = getDestination(arc);
		if (destination != node)
			createAnalogousArc(node, destination, arc, graph);
		totalMult += getMultiplicity(arc);
		destroyArc(arc, graph);
		arc = getArc(twinNode);
	}

	arc = createArc(twinNode, node, graph);
	setMultiplicity(arc, totalMult);

	// Uniqueness
	setUniqueness(node, false);

	// Starting node
	if (startingNode == node)
		startingNode = twinNode;
}

static void remapNodeTimesOntoNeighbour(Node * source, Node * target)
{
	remapNodeTimesOntoTargetNode(source, target);
	remapNodeTimesOntoTargetNode(getTwinNode(source),
				     getTwinNode(target));
}

static void remapNodeTimesOntoForwardMiddlePath(Node * source,
						PassageMarkerI path)
{
	PassageMarkerI marker;
	Node *target;
	Time nodeTime = getNodeTime(source);
	Node *previousNode = getNodePrevious(source);
	Time targetTime;

	//velvetLog("Remapping times from %ld to %ld\n", getNodeID(previousNode), getNodeID(source));

	for (marker = path; getNode(marker) != source;
	     marker = getNextInSequence(marker)) {
		target = getNode(marker);
		targetTime = getNodeTime(target);

		//velvetLog("Through %ld\n", getNodeID(target));

		if (targetTime == -1
		    || targetTime > nodeTime
		    || (targetTime == nodeTime
			&& !isPreviousToNode(target, previousNode))) {
			setNodeTime(target, nodeTime);
			setNodePrevious(previousNode, target);
		}

		previousNode = target;
	}

	setNodePrevious(previousNode, source);

}

static void remapNodeTimesOntoTwinMiddlePath(Node * source,
					     PassageMarkerI path)
{
	PassageMarkerI marker;
	Node *target;
	Node *previousNode = getTwinNode(source);
	Time targetTime;
	PassageMarkerI limit = getTwinMarker(getPreviousInSequence(path));
	Time nodeTime = getNodeTime(getNode(limit));

	//velvetLog("Remapping times from twins %ld to %ld\n", getNodeID(previousNode), getNodeID(getNode(limit)));

	// Revving up
	marker = path;
	while (getNode(marker) != source)
		marker = getNextInSequence(marker);
	marker = getTwinMarker(marker);

	// Going down the path
	while (marker != limit) {
		marker = getNextInSequence(marker);
		target = getNode(marker);
		targetTime = getNodeTime(target);

		//velvetLog("Through %ld\n", getNodeID(target));

		if (targetTime == -1
		    || targetTime > nodeTime
		    || (targetTime == nodeTime
			&& !isPreviousToNode(target, previousNode))) {
			setNodeTime(target, nodeTime);
			getNodeTime(target);
			setNodePrevious(previousNode, target);
		}

		previousNode = target;
	}
}

static void
remapNodeFibHeapReferencesOntoNode(Node * source, Node * target)
{
	DFibHeapNode *sourceDHeapNode = getNodeDHeapNode(source);
	DFibHeapNode *targetDHeapNode = getNodeDHeapNode(target);

	if (sourceDHeapNode == NULL)
		return;

	if (targetDHeapNode == NULL) {
		setNodeDHeapNode(target, sourceDHeapNode);
		replaceValueInDHeap(sourceDHeapNode, target);
	} else if (getKey(targetDHeapNode) > getKey(sourceDHeapNode)) {
		setNodeDHeapNode(target, sourceDHeapNode);
		replaceValueInDHeap(sourceDHeapNode, target);
		destroyNodeInDHeap(targetDHeapNode, dheap);
	} else
		destroyNodeInDHeap(sourceDHeapNode, dheap);

	setNodeDHeapNode(source, NULL);
}

static void remapNodeOntoNeighbour(Node * source,
				   PassageMarkerI sourceMarker,
				   Node * target,
				   PassageMarkerI targetMarker)
{
	//velvetLog("Remapping node %ld onto middle path %ld\n", getNodeID(source), getNodeID(target));
	remapNodeMarkersOntoNeighbour(source, sourceMarker, target,
				      targetMarker);

	remapNodeTimesOntoNeighbour(source, target);
	remapNodeArcsOntoNeighbour(source, target);

	remapNodeFibHeapReferencesOntoNode(getTwinNode(source),
					   getTwinNode(target));
	remapNodeFibHeapReferencesOntoNode(source, target);

	if (startingNode == source)
		startingNode = target;
	if (startingNode == getTwinNode(source))
		startingNode = getTwinNode(target);

	destroyNode(source, graph);
}

static void remapBackOfNodeDescriptorOntoNeighbour(Node * source,
						   PassageMarkerI sourceMarker,
						   Node * target,
						   PassageMarkerI targetMarker,
						   boolean slowToFast,
						   Coordinate offset)
{
	//velvetLog("Splitting node descriptor %ld // %ld\n", getNodeLength(source), offset);

	if (slowToFast)
		splitNodeDescriptor(source, NULL, offset);
	else
		splitNodeDescriptor(source, target, offset);
}

static void remapBackOfNodeTimesOntoNeighbour(Node * source, Node * target)
{
	Time targetTime = getNodeTime(target);
	Time nodeTime = getNodeTime(source);
	Node *twinTarget = getTwinNode(target);
	Node *twinSource = getTwinNode(source);
	Node *previous;

	if (nodeTime != -1) {
		previous = getNodePrevious(source);

		if (previous == source) {
			setNodeTime(target, nodeTime);
			setNodePrevious(target, target);
		} else if (targetTime == -1
			   || targetTime > nodeTime
			   || (targetTime == nodeTime
			       && !isPreviousToNode(target, previous))) {
			setNodeTime(target, nodeTime);
			if (previous != getTwinNode(source))
				setNodePrevious(previous, target);
			else
				setNodePrevious(getTwinNode(target),
						target);
		}

		setNodePrevious(target, source);
	}

	targetTime = getNodeTime(twinTarget);
	nodeTime = getNodeTime(twinSource);

	if (nodeTime != -1) {
		if (targetTime == -1
		    || targetTime > nodeTime
		    || (targetTime == nodeTime
			&& !isPreviousToNode(twinTarget, twinSource))) {
			setNodeTime(twinTarget, nodeTime);
			setNodePrevious(twinSource, twinTarget);
		}
	}

	remapNodeInwardReferencesOntoNode(twinSource, twinTarget);
}

static void
remapBackOfNodeOntoNeighbour(Node * source, PassageMarkerI sourceMarker,
			     Node * target, PassageMarkerI targetMarker,
			     boolean slowToFast)
{
	Coordinate offset;
	//velvetLog("Remapping node %ld onto middle path\n", getNodeID(node));

	offset =
	    remapBackOfNodeMarkersOntoNeighbour(source, sourceMarker,
						target, targetMarker,
						slowToFast);
	remapBackOfNodeDescriptorOntoNeighbour(source, sourceMarker,
					       target, targetMarker,
					       slowToFast, offset);
	remapBackOfNodeTimesOntoNeighbour(source, target);
	remapBackOfNodeArcsOntoNeighbour(source, target);

	remapNodeFibHeapReferencesOntoNode(getTwinNode(source),
					   getTwinNode(target));

	if (getTwinNode(source) == startingNode)
		startingNode = getTwinNode(target);
}

static void remapEmptyPathArcsOntoMiddlePathSimple(PassageMarkerI emptyPath,
						   PassageMarkerI targetPath)
{
	PassageMarkerI pathMarker;
	Node *start = getNode(getPreviousInSequence(emptyPath));
	Node *finish = getNode(emptyPath);
	Node *previousNode = start;
	Node *currentNode;
	Arc *originalArc = getArcBetweenNodes(start, finish, graph);

	for (pathMarker = targetPath; getNode(pathMarker) != finish;
	     pathMarker = getNextInSequence(pathMarker)) {
		currentNode = getNode(pathMarker);
		createAnalogousArc(previousNode, currentNode, originalArc, graph);
		previousNode = currentNode;
	}

	createAnalogousArc(previousNode, finish, originalArc, graph);

	destroyArc(originalArc, graph);
}

static void remapEmptyPathMarkersOntoMiddlePathSimple(PassageMarkerI emptyPath,
						      PassageMarkerI targetPath)
{
	PassageMarkerI marker, newMarker, previousMarker, pathMarker;
	Node *start = getNode(getPreviousInSequence(emptyPath));
	Node *finish = getNode(emptyPath);
	PassageMarkerI oldMarker = getMarker(finish);
	Coordinate markerStart;
	IDnum intersectionLength, twinIntersectionLength;
	ShortReadMarker *intersectionReads =
	    commonNodeReads(start, finish, graph, &intersectionLength);
	ShortReadMarker *twinIntersectionReads =
	    commonNodeReads(getTwinNode(start), getTwinNode(finish), graph,
			    &twinIntersectionLength);

	//velvetLog("SIMPLE %ld\t%ld\t%i\t%i\n", markerCount(finish),
	//       getNodeID(finish), arcCount(finish),
	//       arcCount(getTwinNode(finish)));

	// Destroy link to old nodes
	setMarker(finish, NULL_IDX);

	while (oldMarker != NULL_IDX) {
		marker = oldMarker;
		oldMarker = getNextInNode(marker);
		newMarker = getPreviousInSequence(marker);

		if (getNode(newMarker) != start) {
			transposePassageMarker(marker, finish);
			continue;
		}

		markerStart = getPassageMarkerStart(marker);
		for (pathMarker = targetPath;
		     getNode(pathMarker) != finish;
		     pathMarker = getNextInSequence(pathMarker)) {
			previousMarker = newMarker;

			newMarker =
			    addUncertainPassageMarker
			    (getPassageMarkerSequenceID(marker),
			     getNode(pathMarker));
			setPassageMarkerStatus(newMarker,
					       getPassageMarkerStatus
					       (marker));
			setPassageMarkerStart(newMarker, markerStart);
			setPassageMarkerFinish(newMarker, markerStart);
			setNextInSequence(previousMarker, newMarker);
			setPreviousInSequence(previousMarker, newMarker);

			setStartOffset(newMarker, 0);
			setFinishOffset(newMarker, 0);

		}

		setNextInSequence(newMarker, marker);
		setPreviousInSequence(newMarker, marker);
		transposePassageMarker(marker, finish);
	}

	if (readStartsAreActivated(graph)) {
		for (pathMarker = targetPath;
		     getNode(pathMarker) != finish;
		     pathMarker = getNextInSequence(pathMarker)) {
			// Read starts
			spreadReadIDs(intersectionReads,
				      intersectionLength,
				      getNode(pathMarker), graph);
			spreadReadIDs(twinIntersectionReads,
				      twinIntersectionLength,
				      getTwinNode(getNode(pathMarker)),
				      graph);
		}
	}

	free(intersectionReads);
	free(twinIntersectionReads);
}

static boolean markerFollowsPath(PassageMarkerI marker,
				 PassageMarkerI start,
				 PassageMarkerI finish, Node * stopNode)
{
	PassageMarkerI current, path;

	path = start;
	current = marker;

	while (true) {
		if (current == NULL_IDX || path == finish || path == NULL_IDX)
			return true;

		if (getNode(current) != getNode(path))
			return false;

		current = getNextInSequence(current);
		path = getNextInSequence(path);
	}
}

static PassageMarkerList *getAnchors(PassageMarkerI marker, Node * nodeA,
				     Node * nodeB)
{
	PassageMarkerI current, next;
	Node *twinA = getTwinNode(nodeA);
	Node *twinB = getTwinNode(nodeB);
	PassageMarkerList *result = NULL;

	current = marker;
	while (current != NULL_IDX) {
		next = getNextInSequence(current);
		if (getNode(current) == nodeA && getNode(next) == nodeB) {
			result = newPassageMarkerList(next, result);
		}
		if (getNode(current) == twinB && getNode(next) == twinA) {
			result =
			    newPassageMarkerList(getTwinMarker(current),
						 result);
		}
		current = next;
	}

	return result;
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

static void remapEmptyPathMarkersOntoMiddlePathDevious(PassageMarkerI emptyPath,
						       PassageMarkerI targetPath)
{
	PassageMarkerI marker, newMarker, previousMarker, pathMarker;
	Node *start = getNode(getPreviousInSequence(emptyPath));
	Node *finish = getNode(emptyPath);
	PassageMarkerList *anchors = getAnchors(targetPath, start, finish);
	PassageMarkerList *currentAnchor;
	boolean untouchable = false;
	Coordinate markerStart;

	velvetLog("DEVIOUS %li\t%li\t%li\t%li\n", (long) markerCount(finish),
	       (long) getNodeID(finish), (long) arcCount(finish),
	       (long) arcCount(getTwinNode(finish)));

	for (marker = getMarker(finish); marker != NULL_IDX;
	     marker = getNextInNode(marker)) {
		newMarker = getPreviousInSequence(marker);

		if (getNode(newMarker) != start)
			continue;


		for (currentAnchor = anchors; currentAnchor != NULL;
		     currentAnchor = currentAnchor->next)
			if (markerFollowsPath
			    (marker, currentAnchor->marker, targetPath,
			     finish)) {
				untouchable = true;
				break;
			}

		if (untouchable)
			continue;

		markerStart = getPassageMarkerStart(marker);
		for (pathMarker = targetPath;
		     getNode(pathMarker) != finish;
		     pathMarker = getNextInSequence(pathMarker)) {
			previousMarker = newMarker;
			newMarker =
			    addUncertainPassageMarker
			    (getPassageMarkerSequenceID(marker),
			     getNode(pathMarker));
			setPassageMarkerStatus(newMarker,
					       getPassageMarkerStatus
					       (marker));
			setPassageMarkerStart(newMarker, markerStart);
			setPassageMarkerFinish(newMarker, markerStart);
			setNextInSequence(previousMarker, newMarker);
			setPreviousInSequence(previousMarker, newMarker);

			setStartOffset(newMarker, 0);
			setFinishOffset(newMarker, 0);
		}

		setNextInSequence(newMarker, marker);
		setPreviousInSequence(newMarker, marker);
	}

	destroyPassageMarkerList(&anchors);
}

static boolean markerLeadsToArc(PassageMarkerI marker, Node * nodeA,
				Node * nodeB)
{
	PassageMarkerI current, next;
	Node *twinA = getTwinNode(nodeA);
	Node *twinB = getTwinNode(nodeB);

	current = marker;
	while (current != NULL_IDX) {
		next = getNextInSequence(current);
		if (getNode(current) == nodeA && getNode(next) == nodeB)
			return true;
		if (getNode(current) == twinB && getNode(next) == twinA)
			return true;
		current = next;
	}

	return false;
}

static void
remapEmptyPathOntoMiddlePath(PassageMarkerI emptyPath,
			     PassageMarkerI targetPath)
{
	Node *start = getNode(getPreviousInSequence(emptyPath));
	Node *finish = getNode(emptyPath);

	// Remapping markers
	if (!markerLeadsToArc(targetPath, start, finish)) {
		remapEmptyPathArcsOntoMiddlePathSimple(emptyPath,
						       targetPath);
		remapEmptyPathMarkersOntoMiddlePathSimple(emptyPath,
							  targetPath);
	} else {
		remapEmptyPathMarkersOntoMiddlePathDevious(emptyPath,
							   targetPath);
	}

	//Remap times and previous(if necessary)
	if (getNodePrevious(finish) == start)
		remapNodeTimesOntoForwardMiddlePath(finish, targetPath);

	if (getNodePrevious(getTwinNode(start)) == getTwinNode(finish))
		remapNodeTimesOntoTwinMiddlePath(finish, targetPath);
}

static void reduceSlowNodes(PassageMarkerI slowMarker, Node * finish)
{
	PassageMarkerI marker;

	for (marker = slowMarker; getNode(marker) != finish;
	     marker = getNextInSequence(marker)) {
		reduceNode(getNode(marker));
	}
}

static void destroyPaths()
{
	PassageMarkerI marker;

	while (slowPath != NULL_IDX) {
		marker = slowPath;
		slowPath = getNextInSequence(marker);
		destroyPassageMarker(marker);
	}

	while (fastPath != NULL_IDX) {
		marker = fastPath;
		fastPath = getNextInSequence(marker);
		destroyPassageMarker(marker);
	}
}

static Coordinate mapDistancesOntoPaths()
{
	PassageMarkerI marker;
	Coordinate totalDistance = 0;

	marker = slowPath;
	while (!isTerminal(marker)) {
		marker = getNextInSequence(marker);
		setPassageMarkerStart(marker, totalDistance);
		totalDistance += getNodeLength(getNode(marker));
		setPassageMarkerFinish(marker, totalDistance);
	}

	totalDistance = 0;
	marker = fastPath;
	while (!isTerminal(marker)) {
		marker = getNextInSequence(marker);
		setPassageMarkerStart(marker, totalDistance);
		totalDistance += getNodeLength(getNode(marker));
		setPassageMarkerFinish(marker, totalDistance);
	}

	return totalDistance;
}

static boolean markerLeadsToNode(PassageMarkerI marker, Node * node)
{
	PassageMarkerI currentMarker;

	for (currentMarker = marker; currentMarker != NULL_IDX;
	     currentMarker = getNextInSequence(currentMarker))
		if (getNode(currentMarker) == node)
			return true;

	return false;
}

static void transferNodeData(Node * source, Node * target)
{
	Arc *arc;
	Node *twinSource = getTwinNode(source);
	Node *twinTarget = getTwinNode(target);
	Node *destination;

	// Time & Outward references
	if (getNodePrevious(source) == source) {
		setNodeTime(target, getNodeTime(source));
		setNodePrevious(target, target);
	}

	if (getNodeTime(twinSource) == -1);
	else if (getNodePrevious(twinSource) == twinSource) {
		setNodeTime(twinTarget, getNodeTime(twinSource));
		setNodePrevious(twinTarget, twinTarget);
	} else if (getNodeTime(twinTarget) == -1
		   || getNodeTime(twinSource) < getNodeTime(twinTarget)
		   || (getNodeTime(twinSource) == getNodeTime(twinTarget)
		       && !isPreviousToNode(twinTarget, twinSource))) {
		setNodeTime(twinTarget, getNodeTime(twinSource));
		setNodePrevious(getNodePrevious(twinSource), twinTarget);
	}

	if (getNodePrevious(twinTarget) == source)
		setNodePrevious(target, twinTarget);

	// Inward references:
	for (arc = getArc(source); arc != NULL; arc = getNextArc(arc)) {
		destination = getDestination(arc);
		if (getNodePrevious(destination) == source)
			setNodePrevious(target, destination);
	}

	// Fib Heap refs
	remapNodeFibHeapReferencesOntoNode(source, target);
	remapNodeFibHeapReferencesOntoNode(twinSource, twinTarget);

	// Starting point
	if (startingNode == source)
		startingNode = target;
	else if (startingNode == twinSource)
		startingNode = twinTarget;

	if (getNode(slowPath) == twinSource)
		slowPath = getNextInSequence(slowPath);
	if (getNode(fastPath) == twinSource)
		fastPath = getNextInSequence(fastPath);

	// Next node 
	if (source == activeNode) {
		activeNode = target;
		todo =
		    &todoLists[getNodeID(activeNode) + nodeCount(graph)];
	}
	concatenateTodoLists(target, source);

	if (twinSource == activeNode) {
		activeNode = twinTarget;
		todo =
		    &todoLists[getNodeID(activeNode) + nodeCount(graph)];
	}
}

// Replaces two consecutive nodes into a single equivalent node
// The extra memory is freed
static void concatenateNodesAndVaccinate(Node * nodeA, Node * nodeB,
					 Graph * graph)
{
	PassageMarkerI marker, tmpMarker;
	Node *twinA = getTwinNode(nodeA);
	Node *twinB = getTwinNode(nodeB);
	Arc *arc;

	//velvetLog("Concatenating nodes %ld and %ld\n", getNodeID(nodeA), getNodeID(nodeB));
	// Arc management:
	// Freeing useless arcs
	while (getArc(nodeA) != NULL)
		destroyArc(getArc(nodeA), graph);

	// Correct arcs
	for (arc = getArc(nodeB); arc != NULL; arc = getNextArc(arc)) {
		if (getDestination(arc) != twinB)
			createAnalogousArc(nodeA, getDestination(arc), arc, graph);
		else
			createAnalogousArc(nodeA, twinA, arc, graph);
	}

	// Passage marker management in node A:
	for (marker = getMarker(nodeA); marker != NULL_IDX;
	     marker = getNextInNode(marker))
		if (isTerminal(marker))
			incrementFinishOffset(marker,
					      getNodeLength(nodeB));

	// Swapping new born passageMarkers from B to A
	for (marker = getMarker(nodeB); marker != NULL_IDX; marker = tmpMarker) {
		tmpMarker = getNextInNode(marker);

		if (isInitial(marker)) {
			extractPassageMarker(marker);
			transposePassageMarker(marker, nodeA);
			incrementStartOffset(marker, getNodeLength(nodeA));
		} else
			disconnectNextPassageMarker(getPreviousInSequence
						    (marker), graph);
	}

	// Read starts
	concatenateReadStarts(nodeA, nodeB, graph);

	// Descriptor management 
	appendDescriptors(nodeA, nodeB);

	// Update uniqueness:
	setUniqueness(nodeA, getUniqueness(nodeA) || getUniqueness(nodeB));

#ifndef SINGLE_COV_CAT
	Category cat;
	for (cat = 0; cat < CATEGORIES; cat++) {
		// Update virtual coverage
		incrementVirtualCoverage(nodeA, cat,
					 getVirtualCoverage(nodeB, cat));
		// Update virtual coverage
		incrementOriginalVirtualCoverage(nodeA, cat,
						 getOriginalVirtualCoverage(nodeB, cat));
	}
#else
	incrementVirtualCoverage(nodeA, getVirtualCoverage(nodeB));
#endif

	// Freeing gobbled node
	destroyNode(nodeB, graph);
}

static void simplifyNode(Node * node)
{
	Node *twin = getTwinNode(node);
	Node *destination, *twinDestination;

	if (!hasSingleArc(node))
		return;

	destination = getDestination(getArc(node));
	twinDestination = getTwinNode(destination);

	while (hasSingleArc(node)
	       && hasSingleArc(twinDestination)
	       && destination != twin && destination != node) {
		transferNodeData(destination, node);
		concatenateNodesAndVaccinate(node, destination, graph);

		if (!hasSingleArc(node))
			return;
		destination = getDestination(getArc(node));
		twinDestination = getTwinNode(destination);
	}

}

static void concatenatePathNodes(PassageMarkerI pathStart)
{
	PassageMarkerI pathMarker;

	//velvetLog("Removing null loops\n");
	for (pathMarker = pathStart; pathMarker != NULL_IDX;
	     pathMarker = getNextInSequence(pathMarker)) {
		simplifyNode(getNode(pathMarker));
	}
}

#define SLOW_TO_FAST true
#define FAST_TO_SLOW false

static void cleanUpRedundancy()
{
	PassageMarkerI slowMarker = getNextInSequence(slowPath);
	PassageMarkerI fastMarker = getNextInSequence(fastPath);
	Coordinate slowLength, fastLength;
	Coordinate fastConstraint = 0;
	Coordinate slowConstraint = 0;
	Coordinate finalLength;
	Node *slowNode, *fastNode;

	//velvetLog("Correcting new redundancy\n");
	mapSlowOntoFast();
	finalLength = mapDistancesOntoPaths();

	while (slowMarker != NULL_IDX && fastMarker != NULL_IDX) {
		if (isTerminal(slowMarker))
			slowLength = finalLength;
		else {
			slowLength =
			    slowToFastMapping[getPassageMarkerFinish
					      (slowMarker) - 1];
			if (slowLength < slowConstraint)
				slowLength = slowConstraint;
		}

		fastLength = getPassageMarkerFinish(fastMarker) - 1;
		if (fastLength < fastConstraint)
			fastLength = fastConstraint;

		slowNode = getNode(slowMarker);
		fastNode = getNode(fastMarker);

		if (slowNode == fastNode) {
			if (fastLength > slowLength)
				slowConstraint = fastLength;
			else if (fastLength < slowLength)
				fastConstraint = slowLength;

			slowMarker = getNextInSequence(slowMarker);
			fastMarker = getNextInSequence(fastMarker);
		} else if (slowNode == getTwinNode(fastNode)) {
			if (fastLength > slowLength)
				slowConstraint = fastLength;
			else if (fastLength < slowLength)
				fastConstraint = slowLength;

			slowMarker = getNextInSequence(slowMarker);
			fastMarker = getNextInSequence(fastMarker);
			foldSymmetricalNode(fastNode);
		} else if (markerLeadsToNode(slowMarker, fastNode)) {
			reduceSlowNodes(slowMarker, fastNode);
			remapEmptyPathOntoMiddlePath(fastMarker,
						     slowMarker);
			while (getNode(slowMarker) != fastNode)
				slowMarker = getNextInSequence(slowMarker);
		} else if (markerLeadsToNode(fastMarker, slowNode)) {
			remapEmptyPathOntoMiddlePath(slowMarker,
						     fastMarker);
			while (getNode(fastMarker) != slowNode)
				fastMarker = getNextInSequence(fastMarker);
		} else if (slowLength == fastLength) {
			remapNodeOntoNeighbour(slowNode, slowMarker,
					       fastNode, fastMarker);
			slowMarker = getNextInSequence(slowMarker);
			fastMarker = getNextInSequence(fastMarker);
		} else if (slowLength < fastLength) {
			remapBackOfNodeOntoNeighbour(fastNode, fastMarker,
						     slowNode, slowMarker,
						     FAST_TO_SLOW);
			slowMarker = getNextInSequence(slowMarker);
		} else {
			remapBackOfNodeOntoNeighbour(slowNode, slowMarker,
						     fastNode, fastMarker,
						     SLOW_TO_FAST);
			fastMarker = getNextInSequence(fastMarker);
		}

		fflush(stdout);
	}

	//velvetLog("Done with path\n");

	while (!isInitial(slowPath))
		slowPath = getPreviousInSequence(slowPath);
	while (!isInitial(fastPath))
		fastPath = getPreviousInSequence(fastPath);

	//velvetLog("Concatenation\n");

	// Freeing up memory  
	if (slowMarker != NULL_IDX)
		concatenatePathNodes(slowPath);
	else
		concatenatePathNodes(fastPath);

	//velvetLog("Vaccinatting\n");

	destroyPaths();

	// Cleaning up silly structures
	//vaccinatePath(&returnValue);

	//velvetLog("Clean up done\n");
	//fflush(stdout);
}

static boolean pathContainsReference(PassageMarkerI path) {
	PassageMarkerI marker, marker2;

	for (marker = getNextInSequence(path); !isTerminal(marker);
	     marker = getNextInSequence(marker))
		for (marker2 = getMarker(getNode(marker)); marker2 != NULL_IDX; marker2 = getNextInNode(marker2))
			if (marker2 != marker && sequenceCategories[getAbsolutePassMarkerSeqID(marker2) - 1] == REFERENCE)
				return true;

	return false;

}

static void comparePaths(Node * destination, Node * origin)
{
	IDnum slowLength, fastLength;
	Node *fastNode, *slowNode;
	IDnum i;
	PassageMarkerI marker;

	//Measure lengths
	slowLength = fastLength = 0;
	fastNode = destination;
	slowNode = origin;

	while (fastNode != slowNode) {
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
		    || fastLength > MAXNODELENGTH)
			return;
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

	// Avoid merging parallel Reference sequences
	if (pathContainsReference(fastPath) && pathContainsReference(slowPath)) {
		destroyPaths();
		return;
	}
	//Extract sequences
	if (!extractSequence(fastPath, fastSequence)
	    || !extractSequence(slowPath, slowSequence)) {
		destroyPaths();
		return;
	}
	//Compare sequences
	if (compareSequences(fastSequence, slowSequence)) {
		cleanUpRedundancy();
		return;
	}
	//velvetLog("\tFinished comparing paths, changes made\n");
	destroyPaths();
}

static void tourBusArc(Node * origin, Arc * arc, Time originTime)
{
	Node *destination = getDestination(arc);
	Time arcTime, totalTime, destinationTime;
	IDnum nodeIndex = getNodeID(destination) + nodeCount(graph);
	Node *oldPrevious = previous[nodeIndex];

	if (oldPrevious == origin || getNodeStatus(destination) == 1)
		return;

	arcTime =
	    ((Time) getNodeLength(origin)) / ((Time) getMultiplicity(arc));
	totalTime = originTime + arcTime;

	destinationTime = times[nodeIndex];

	if (destinationTime == -1) {
		setNodeTime(destination, totalTime);
		dheapNodes[nodeIndex] =
		    insertNodeIntoDHeap(dheap, totalTime, destination);
		previous[nodeIndex] = origin;
		return;
	} else if (destinationTime > totalTime) {
		if (dheapNodes[nodeIndex] == NULL) {
			//velvetLog("Already expanded though\n");
			return;
		}

		setNodeTime(destination, totalTime);
		replaceKeyInDHeap(dheap, dheapNodes[nodeIndex], totalTime);
		previous[nodeIndex] = origin;

		comparePaths(destination, oldPrevious);
		return;
	} else {
		if (destinationTime == getNodeTime(origin)
		    && isPreviousToNode(destination, origin)) {
			return;
		}

		comparePaths(destination, origin);
	}
}

static void initializeTodoLists()
{
	IDnum index;
	Node *node;
	Arc *arc;
	Ticket *tkt;
	IDnum nodes = nodeCount(graph);
	Ticket **currentList;
	Ticket *currentTicket, *tmp;
	Node *destination;

	velvetLog("Initializing todo lists\n");

	for (index = -nodes; index <= nodes; index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL)
			continue;

		currentList = &todoLists[index + nodes];
		*currentList = NULL;

		for (arc = getArc(node); arc != NULL;
		     arc = getNextArc(arc)) {
			destination = getDestination(arc);

			if (destination == node)
				continue;

			tkt = newTicket();
			tkt->id_a = getNodeID(destination);

			currentTicket = *currentList;
			if (currentTicket == NULL
			    || currentTicket->id_a > tkt->id_a) {
				tkt->next = currentTicket;
				*currentList = tkt;
				continue;
			}

			while (currentTicket->next != NULL
			       && currentTicket->next->id_a < tkt->id_a)
				currentTicket = currentTicket->next;

			tmp = currentTicket->next;
			currentTicket->next = tkt;
			tkt->next = tmp;
		}
	}

	velvetLog("Done with initilization\n");
}

static void tourBusNode(Node * node)
{
	Arc *arc;
	Node *destination;
	Time nodeTime = getNodeTime(node);
	IDnum id;

	dbgCounter++;
	if (dbgCounter % 10000 == 0) {
		velvetLog("%li / %li nodes visited\n", (long) dbgCounter, (long) nodeCount(graph));
		fflush(stdout);
	}

	setSingleNodeStatus(node, 2);
	activeNode = node;
	todo = &todoLists[getNodeID(activeNode) + nodeCount(graph)];
	done = NULL;

	while ((id = nextTodoTicket()) != 0) {
		destination = getNodeInGraph(graph, id);

		// Node doesn't exist anymore
		if (destination == NULL)
			continue;

		arc = getArcBetweenNodes(activeNode, destination, graph);

		// Arc does not exist for some reason (?)
		if (arc == NULL)
			continue;

		tourBusArc(activeNode, arc, nodeTime);
	}

	freeDoneTickets();
}

static Coordinate getTipLength(Node * node)
{
	Node *current = getTwinNode(node);
	Coordinate length = 0;

	if (simpleArcCount(current) > 1)
		return getNodeLength(node);

	while (current != NULL && simpleArcCount(getTwinNode(current)) < 2
	       && simpleArcCount(current) < 2) {
		length += getNodeLength(current);
		current = getDestination(getArc(current));
	}

	return length;
}

void clipTipsHard(Graph * graph, boolean conserveLong)
{
	IDnum index;
	Node *current, *twin;
	boolean modified = true;
	int Wordlength = getWordLength(graph);
	PassageMarkerI marker;

	velvetLog("Clipping short tips off graph, drastic\n");

	while (modified) {
		modified = false;
		for (index = 1; index <= nodeCount(graph); index++) {
			current = getNodeInGraph(graph, index);

			if (current == NULL)
				continue;
	
			if (conserveLong && getMarker(current))
				continue;

			twin = getTwinNode(current);

			if (getArc(current) == NULL
			    && getTipLength(current) < 2 * Wordlength) {
				while ((marker = getMarker(current))) {
					if (!isInitial(marker)
					    && !isTerminal(marker))
						deleteNextPassageMarker
						    (getPreviousInSequence
						     (marker), graph);
					destroyPassageMarker(marker);
				}
				destroyNode(current, graph);
				modified = true;
			} else if (getArc(twin) == NULL
				   && getTipLength(twin) <
				   2 * Wordlength) {
				while ((marker = getMarker(current))) {
					if (!isInitial(marker)
					    && !isTerminal(marker))
						deleteNextPassageMarker
						    (getPreviousInSequence
						     (marker), graph);
					destroyPassageMarker(marker);
				}
				destroyNode(twin, graph);
				modified = true;
			}
		}
	}

	concatenateGraph(graph);
	velvetLog("%li nodes left\n", (long) nodeCount(graph));
}

static void tourBus(Node * startingPoint)
{
	Node *currentNode = startingPoint;
	IDnum nodeID = getNodeID(startingPoint) + nodeCount(graph);

	//velvetLog("Tour bus from node %ld...\n", (long) getNodeID(startingPoint));

	times[nodeID] = 0;
	previous[nodeID] = currentNode;

	while (currentNode != NULL) {
		dheapNodes[getNodeID(currentNode) + nodeCount(graph)] =
		    NULL;
		tourBusNode(currentNode);
		currentNode = removeNextNodeFromDHeap(dheap);
	}
}

void correctGraph(Graph * argGraph, ShortLength * argSequenceLengths, Category * argSequenceCategories, boolean conserveLong)
{
	IDnum nodes;
	IDnum index;
	double *FmatrixMem;

	//Setting global params
	graph = argGraph;
	WORDLENGTH = getWordLength(graph);
	sequenceLengths = argSequenceLengths;
	sequenceCategories = argSequenceCategories;
	dbgCounter = 0;
	// Done with global params

	velvetLog("Correcting graph with cutoff %f\n", MAXDIVERGENCE);

	nodes = nodeCount(graph);

	// Allocating memory
	times = mallocOrExit(2 * nodes + 1, Time);
	previous = mallocOrExit(2 * nodes + 1, Node *);
	dheapNodes = mallocOrExit(2 * nodes + 1, DFibHeapNode *);

	for (index = 0; index < (2 * nodeCount(graph) + 1); index++) {
		times[index] = -1;
		previous[index] = NULL;
		dheapNodes[index] = NULL;
	}

	dheap = newDFibHeap();

	fastSequence = newTightString(MAXREADLENGTH);
	slowSequence = newTightString(MAXREADLENGTH);
	fastToSlowMapping = callocOrExit(MAXREADLENGTH + 1, Coordinate);
	slowToFastMapping = callocOrExit(MAXREADLENGTH + 1, Coordinate);
	Fmatrix = mallocOrExit(MAXREADLENGTH + 1, double *);
	FmatrixMem = callocOrExit((MAXREADLENGTH + 1) * (MAXREADLENGTH + 1), double);
	for (index = 0; index < MAXREADLENGTH + 1; index++)
		Fmatrix[index] = FmatrixMem + index * (MAXREADLENGTH + 1);

	eligibleStartingPoints = mallocOrExit(2 * nodes + 1, IDnum);
	progressStatus = callocOrExit(2 * nodes + 1, boolean);
	todoLists = callocOrExit(2 * nodes + 1, Ticket *);
	//Done with memory 

	resetNodeStatus(graph);
	determineEligibleStartingPoints();
	initializeTodoLists();
	activateArcLookupTable(graph);

	while ((startingNode = nextStartingPoint()) != NULL) {
		//velvetLog("Going through the cycle...\n");
		tourBus(startingNode);
		updateNodeStatus(startingNode);
	}

	deactivateArcLookupTable(graph);
	concatenateGraph(graph);

	clipTipsHard(graph, conserveLong);

	//Deallocating globals
	free(times);
	free(previous);
	free(dheapNodes);
	destroyDHeap(dheap);

	destroyTightString(fastSequence);
	destroyTightString(slowSequence);
	free(fastToSlowMapping);
	free(slowToFastMapping);
	free(Fmatrix);
	free(FmatrixMem);

	free(eligibleStartingPoints);
	free(progressStatus);
	free(todoLists);

	if (ticketMemory != NULL)
		destroyRecycleBin(ticketMemory);

	free(sequenceLengths);
	//Done deallocating
}

void setMaxReadLength(int value)
{
	if (value < 0) {
		velvetLog("Negative branch length %i!\n", value);
		velvetLog("Exiting...\n");
#ifdef DEBUG 
		abort();
#endif 
		exit(1);
	}
	MAXREADLENGTH = value;
	MAXNODELENGTH = 2 * value;
}

void setMaxGaps(int value)
{
	if (value < 0) {
		velvetLog("Negative max gap count %i!\n", value);
		velvetLog("Exiting...\n");
#ifdef DEBUG 
		abort();
#endif 
		exit(1);
	}
	MAXGAPS = value;
}

void setMaxDivergence(double value)
{
	if (value < 0 || value > 1) {
		velvetLog("Divergence rate %lf out of bounds [0,1]!\n",
		       value);
		velvetLog("Exiting...\n");
#ifdef DEBUG 
		abort();
#endif 
		exit(1);
	}
	MAXDIVERGENCE = value;
}
