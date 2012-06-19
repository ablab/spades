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
#include "recycleBin.h"
#include "passageMarker.h"
#include "graphStats.h"
#include "concatenatedGraph.h"
#include "readSet.h"
#include "utility.h"

#define LONG_NODE_CUTOFF 50
#define LN2 0.693147
#define PROBABILITY_CUTOFF 5
#define MAX_READ_COUNT 100
#define MAX_READ_LENGTH 2000

static Graph *graph = NULL;
static PassageMarkerI path = NULL_IDX;
static RecycleBin *listMemory = NULL;
static double expected_coverage = 1;
static TightString *sequences = NULL;
static int MULTIPLICITY_CUTOFF = 2;

static IDnum multCounter = 0;
static IDnum dbgCounter = 0;
static IDnum nullCounter = 0;

typedef struct rb_connection_st RBConnection;

struct rb_connection_st {
	Node *node;
	PassageMarkerI marker;
	RBConnection *next;
	IDnum multiplicity;
}  ATTRIBUTE_PACKED;

static RecycleBin *nodeListMemory = NULL;

#define BLOCKSIZE 1000

static RBConnection *allocateRBConnection()
{
	if (nodeListMemory == NULL)
		nodeListMemory =
		    newRecycleBin(sizeof(RBConnection), BLOCKSIZE);

	return (RBConnection*)allocatePointer(nodeListMemory);
}

static void deallocateRBConnection(RBConnection * nodeList)
{
	deallocatePointer(nodeListMemory, nodeList);
}

boolean isUniqueSolexa(Node * node)
{

	Coordinate nodeLength = getNodeLength(node);
	Coordinate nodeCoverage;
	double nodeDensity, probability;

	nodeCoverage = getTotalCoverage(node);

	if (nodeLength > LONG_NODE_CUTOFF) {
		nodeDensity = nodeCoverage / (double) nodeLength;

		probability =
		    LN2 / 2 +
		    nodeLength / (2 * expected_coverage) *
		    (expected_coverage * expected_coverage -
		     nodeDensity * nodeDensity / 2);
		return probability > PROBABILITY_CUTOFF;
	}
	return false;
}

static void identifyUniqueNodes(boolean(*isUniqueFunction) (Node *))
{
	IDnum index;
	Node *node;
	IDnum counter = 0;

	velvetLog("Identifying unique nodes\n");

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (node == NULL)
			continue;

		setUniqueness(node, isUniqueFunction(node));

		if (getUniqueness(node))
			counter++;
	}

	velvetLog("Done, %li unique nodes counted\n", (long) counter);
}

boolean uniqueNodesConnect(Node * startingNode)
{
	Node *destination = NULL;
	PassageMarkerI startMarker, currentMarker;
	RBConnection *newList;
	RBConnection *list = NULL;
	boolean multipleHits = false;

	if (arcCount(startingNode) == 0)
		return false;

	if (getMarker(startingNode) == NULL_IDX)
		return false;

	dbgCounter++;

	// Checking for multiple destinations
	for (startMarker = getMarker(startingNode); startMarker != NULL_IDX;
	     startMarker = getNextInNode(startMarker)) {
		if (getFinishOffset(startMarker) >
		    2 * getWordLength(graph))
			continue;

		for (currentMarker = getNextInSequence(startMarker);
		     currentMarker != NULL_IDX;
		     currentMarker = getNextInSequence(currentMarker)) {
			if (!getUniqueness(getNode(currentMarker))) {
				continue;
			} else if (getNodeStatus(getNode(currentMarker))) {
				if (getStartOffset(currentMarker) >
				    2 * getWordLength(graph))
					break;
				for (newList = list; newList != NULL;
				     newList = newList->next) {
					if (newList->node ==
					    getNode(currentMarker)) {
						newList->multiplicity++;
						break;
					}
				}
				if (newList == NULL)
					abort();
				break;
			} else {
				if (getStartOffset(currentMarker) >
				    2 * getWordLength(graph))
					break;
				setSingleNodeStatus(getNode(currentMarker),
						    true);
				newList = allocateRBConnection();
				newList->node = getNode(currentMarker);
				newList->multiplicity = 1;
				newList->marker = startMarker;
				newList->next = list;
				list = newList;
				break;
			}
		}
	}

	while (list != NULL) {
		newList = list;
		list = newList->next;
		setSingleNodeStatus(newList->node, false);
		if (newList->multiplicity >= MULTIPLICITY_CUTOFF) {
			if (destination == NULL) {
				destination = newList->node;
				path = newList->marker;
			} else if (destination != newList->node)
				multipleHits = true;
		}
		deallocateRBConnection(newList);
	}

	if (multipleHits) {
		multCounter++;
		setUniqueness(startingNode, false);
		return false;
	}

	if (destination == NULL || destination == startingNode
	    || destination == getTwinNode(startingNode)) {
		nullCounter++;
		return false;
	}
	// Check for reciprocity
	for (startMarker = getMarker(getTwinNode(destination));
	     startMarker != NULL_IDX;
	     startMarker = getNextInNode(startMarker)) {
		if (getFinishOffset(startMarker) >
		    2 * getWordLength(graph))
			continue;

		for (currentMarker = getNextInSequence(startMarker);
		     currentMarker != NULL_IDX;
		     currentMarker = getNextInSequence(currentMarker)) {
			if (!getUniqueness(getNode(currentMarker))) {
				continue;
			} else if (getNodeStatus(getNode(currentMarker))) {
				if (getStartOffset(currentMarker) >
				    2 * getWordLength(graph))
					break;
				for (newList = list; newList != NULL;
				     newList = newList->next) {
					if (newList->node ==
					    getNode(currentMarker)) {
						newList->multiplicity++;
						break;
					}
				}
				if (newList == NULL)
					abort();
				break;
			} else {
				if (getStartOffset(currentMarker) >
				    2 * getWordLength(graph))
					break;
				setSingleNodeStatus(getNode(currentMarker),
						    true);
				newList = allocateRBConnection();
				newList->node = getNode(currentMarker);
				newList->multiplicity = 1;
				newList->next = list;
				list = newList;
				break;
			}
		}
	}

	while (list != NULL) {
		newList = list;
		list = newList->next;
		setSingleNodeStatus(newList->node, false);
		if (newList->multiplicity >= MULTIPLICITY_CUTOFF
		    && newList->node != getTwinNode(startingNode))
			multipleHits = true;
		deallocateRBConnection(newList);
	}

	if (multipleHits) {
		multCounter++;
		setUniqueness(destination, false);
		return false;
	}
	// Aligning long reads to each other:
	// TODO 

	// Merge pairwise alignments and produce consensus
	// TODO

	return true;
}

static boolean goesToNode(PassageMarkerI marker, Node * node)
{
	PassageMarkerI current;
	Node * start = getNode(marker);
	Node * twinStart = getTwinNode(start);
	Node * currentNode;

	for (current = getNextInSequence(marker); current != NULL_IDX;
	     current = getNextInSequence(current)) {
		currentNode = getNode(current);
		if (currentNode == start || currentNode == twinStart)
			return false;
		else if (currentNode == node)
			return true;
		else if (getUniqueness(currentNode))
			return false;
	}

	return false;
}

static void updateMembers(Node * bypass, Node * nextNode)
{
	PassageMarkerI marker, next, tmp;
	Coordinate nextLength = getNodeLength(nextNode);

	// Update  marker + arc info
	for (marker = getMarker(bypass); marker != NULL_IDX; marker = tmp) {
		tmp = getNextInNode(marker);

		if (!isTerminal(marker)
		    && getNode(getNextInSequence(marker)) == nextNode) {
			// Marker steps right into target
			next = getNextInSequence(marker);
			disconnectNextPassageMarker(marker, graph);
			destroyPassageMarker(next);
		} else if (getUniqueness(nextNode)
			   && goesToNode(marker, nextNode)) {
			// Marker goes indirectly to target
			while (getNode(getNextInSequence(marker)) !=
			       nextNode) {
				next = getNextInSequence(marker);
				disconnectNextPassageMarker(marker, graph);
				destroyPassageMarker(next);
			}

			next = getNextInSequence(marker);
			disconnectNextPassageMarker(marker, graph);
			destroyPassageMarker(next);
		} else if (!isTerminal(marker)
			   && getFinishOffset(marker) == 0) {
			// Marker goes somewhere else than to target
			next = getNextInSequence(marker);
			incrementFinishOffset(marker, nextLength);
		} else {
			// Marker goes nowhere
			incrementFinishOffset(marker, nextLength);
		}
	}
}

static void admitGroupies(Node * source, Node * bypass)
{
	PassageMarkerI marker, tmpMarker;

	for (marker = getMarker(source); marker != NULL_IDX;
	     marker = tmpMarker) {
		tmpMarker = getNextInNode(marker);
		extractPassageMarker(marker);
		transposePassageMarker(marker, bypass);
		incrementFinishOffset(getTwinMarker(marker),
				      getNodeLength(bypass));
	}

}

static void adjustShortReads(Node * target, PassageMarkerI pathMarker)
{
	ShortReadMarker *targetArray, *marker;
	IDnum targetLength, index;
	Coordinate position, nodeLength;

	if (!readStartsAreActivated(graph))
		return;

	targetArray = getNodeReads(getTwinNode(target), graph);
	targetLength = getNodeReadCount(getTwinNode(target), graph);

	nodeLength = getPassageMarkerLength(pathMarker);

	for (index = 0; index < targetLength; index++) {
		marker = getShortReadMarkerAtIndex(targetArray, index);
		position = getShortReadMarkerPosition(marker);
		position += nodeLength;
		setShortReadMarkerPosition(marker, position);
	}
}

Node *bypass()
{
	Node *bypass = getNode(path);
	Node *next = NULL;
	Arc *arc;
	PassageMarkerI nextMarker;

	// Remove unwanted arcs
	while (getArc(bypass) != NULL)
		destroyArc(getArc(bypass), graph);

	// Update extensive variables (length + descriptors + passage markers)
	while (!isTerminal(path)) {
		nextMarker = getNextInSequence(path);
		next = getNode(nextMarker);
		while (next == bypass) {
			disconnectNextPassageMarker(path, graph);
			destroyPassageMarker(nextMarker);
			nextMarker = getNextInSequence(path);
			next = getNode(nextMarker);
		}

		if (next == NULL)
			return bypass;

		// Overall node update 
		if (!getUniqueness(next)) {
			adjustShortReads(bypass, getNextInSequence(path));
			appendSequence(bypass, sequences,
				       getNextInSequence(path), graph);
		} else {
			concatenateReadStarts(bypass, next, graph);

#ifndef SINGLE_COV_CAT
			Category cat;
			for (cat = 0; cat < CATEGORIES; cat++) {
				// Update virtual coverage
				incrementVirtualCoverage(bypass, cat,
							 getVirtualCoverage(next, cat));
				// Update original virtual coverage
				incrementOriginalVirtualCoverage(bypass, cat,
								 getOriginalVirtualCoverage(next, cat));
			}
#else
			incrementVirtualCoverage(bypass, getVirtualCoverage(next));
#endif
			appendDescriptors(bypass, next);
		}

		// Members
		updateMembers(bypass, next);

		// Termination 
		if (isTerminal(path) || getUniqueness(next))
			break;
	}

	// Remove unique groupies from arrival 
	admitGroupies(next, bypass);

	// Copy destination arcs
	for (arc = getArc(next); arc != NULL; arc = getNextArc(arc)) {
		if (getDestination(arc) == next)
			continue;
		else if (getDestination(arc) == getTwinNode(next))
			createAnalogousArc(bypass, getTwinNode(bypass),
					   arc, graph);
		else
			createAnalogousArc(bypass, getDestination(arc),
					   arc, graph);
	}

	destroyNode(next, graph);

	return bypass;
}

static void trimLongReadTips()
{
	IDnum index;
	Node *node;
	PassageMarkerI marker, next;

	velvetLog("Trimming read tips\n");

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (getUniqueness(node))
			continue;

		for (marker = getMarker(node); marker != NULL_IDX;
		     marker = next) {
			next = getNextInNode(marker);

			if (!isInitial(marker) && !isTerminal(marker))
				continue;

			if (isTerminal(marker))
				marker = getTwinMarker(marker);

			while (!getUniqueness(getNode(marker))) {
				if (next != NULL_IDX
				    && (marker == next
					|| marker == getTwinMarker(next)))
					next = getNextInNode(next);
				if (getNextInSequence(marker) != NULL_IDX) {
					marker = getNextInSequence(marker);
					destroyPassageMarker
					    (getPreviousInSequence
					     (marker));
				} else {
					destroyPassageMarker(marker);
					break;
				}
			}
		}
	}
}

void readCoherentGraph(Graph * inGraph, boolean(*isUnique) (Node * node),
		       double coverage, ReadSet * reads)
{
	IDnum nodeIndex;
	Node *node;
	IDnum previousNodeCount = 0;

	graph = inGraph;
	listMemory = newRecycleBin(sizeof(PassageMarkerList), 100000);
	expected_coverage = coverage;
	sequences = reads->tSequences;

	velvetLog("Read coherency...\n");
	resetNodeStatus(graph);
	identifyUniqueNodes(isUnique);
	trimLongReadTips();

	previousNodeCount = 0;
	while (previousNodeCount != nodeCount(graph)) {

		previousNodeCount = nodeCount(graph);

		for (nodeIndex = 1; nodeIndex <= nodeCount(graph);
		     nodeIndex++) {

			node = getNodeInGraph(graph, nodeIndex);

			if (node == NULL || !getUniqueness(node))
				continue;

			while (uniqueNodesConnect(node))
				node = bypass();

			node = getTwinNode(node);

			while (uniqueNodesConnect(node))
				node = bypass();

		}

		renumberNodes(graph);
	}

	destroyRecycleBin(listMemory);
	destroyRecycleBin(nodeListMemory);

	velvetLog("Confronted to %li multiple hits and %li null over %li\n",
	       (long) multCounter, (long) nullCounter, (long) dbgCounter);

	velvetLog("Read coherency over!\n");
}

void setMultiplicityCutoff(int value)
{
	if (value < 0) {
		velvetLog("Negative long read multiplicity cutoff %i!\n",
		       value);
		velvetLog("Exiting...\n");
#ifdef DEBUG 
		abort();
#endif 
		exit(1);
	}
	MULTIPLICITY_CUTOFF = value;
}
