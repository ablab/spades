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

#define LONG_NODE_CUTOFF 50
#define LN2 0.693147
#define PROBABILITY_CUTOFF 5
#define MAX_READ_COUNT 100
#define MAX_READ_LENGTH 2000

static Graph *graph = NULL;
static PassageMarker *path = NULL;
static RecycleBin *listMemory = NULL;
static double expected_coverage = 1;
static TightString **sequences = NULL;
static int MULTIPLICITY_CUTOFF = 2;

static IDnum multCounter = 0;
static IDnum dbgCounter = 0;
static IDnum nullCounter = 0;

typedef struct rb_connection_st RBConnection;

struct rb_connection_st {
	Node *node;
	PassageMarker *marker;
	RBConnection *next;
	IDnum multiplicity;
};

static RecycleBin *nodeListMemory = NULL;

#define BLOCKSIZE 1000

static RBConnection *allocateRBConnection()
{
	if (nodeListMemory == NULL)
		nodeListMemory =
		    newRecycleBin(sizeof(RBConnection), BLOCKSIZE);

	return allocatePointer(nodeListMemory);
}

static void deallocateRBConnection(RBConnection * nodeList)
{
	deallocatePointer(nodeListMemory, nodeList);
}

void setBaseCoverage(double coverage)
{
	expected_coverage = coverage;
}

boolean isUniqueBasic(Node * node)
{
	if (getNodeLength(node) < LONG_NODE_CUTOFF) {
		return false;
	}
	if (readCoverage(node) / (double) getNodeLength(node) >
	    1.5 * expected_coverage) {
		return false;
	}

	return true;
}

boolean isUniqueSolexa(Node * node)
{

	Coordinate nodeLength = getNodeLength(node);
	Coordinate nodeCoverage =
	    (getVirtualCoverage(node, 0) + getVirtualCoverage(node, 1));
	double nodeDensity, probability;

	if (nodeLength == 0) {
		return false;
	}
	if (nodeLength > LONG_NODE_CUTOFF) {
		nodeDensity = nodeCoverage / (double) nodeLength;

		probability =
		    LN2 / 2 +
		    nodeLength / (2 * expected_coverage) *
		    (expected_coverage * expected_coverage -
		     nodeDensity * nodeDensity / 2);
		return probability > PROBABILITY_CUTOFF;
	} else {
		return false;
		probability =
		    expected_coverage * nodeLength - nodeCoverage / LN2;
		return probability > 0;
	}
}

static void identifyUniqueNodes(boolean(*isUniqueFunction) (Node *))
{
	IDnum index;
	Node *node;
	IDnum counter = 0;

	puts("Identifying unique nodes");

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (node == NULL)
			continue;

		setUniqueness(node, isUniqueFunction(node));

		if (getUniqueness(node))
			counter++;
	}

	printf("Done, %u unique nodes counted\n", counter);
}

static boolean uniqueNodesConnect(Node * startingNode)
{
	Node *destination = NULL;
	PassageMarker *startMarker, *currentMarker;
	RBConnection *newList;
	RBConnection *list = NULL;
	boolean multipleHits = false;

	if (arcCount(startingNode) == 0)
		return false;

	if (getMarker(startingNode) == NULL)
		return false;

	dbgCounter++;

	// Checking for multiple destinations
	for (startMarker = getMarker(startingNode); startMarker != NULL;
	     startMarker = getNextInNode(startMarker)) {
		if (getFinishOffset(startMarker) >
		    2 * getWordLength(graph))
			continue;

		for (currentMarker = getNextInSequence(startMarker);
		     currentMarker != NULL;
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
	     startMarker != NULL;
	     startMarker = getNextInNode(startMarker)) {
		if (getFinishOffset(startMarker) >
		    2 * getWordLength(graph))
			continue;

		for (currentMarker = getNextInSequence(startMarker);
		     currentMarker != NULL;
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

static boolean goesToNode(PassageMarker * marker, Node * node)
{
	PassageMarker *current;

	for (current = marker; current != NULL;
	     current = getNextInSequence(current))
		if (getNode(current) == node)
			return true;

	return false;
}

static void updateMembers(Node * bypass, Node * nextNode)
{
	PassageMarker *marker, *next, *tmp;
	Coordinate nextLength = getNodeLength(nextNode);

	// Update  marker + arc info
	for (marker = getMarker(bypass); marker != NULL; marker = tmp) {
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
	PassageMarker *marker, *tmpMarker;

	for (marker = getMarker(source); marker != NULL;
	     marker = tmpMarker) {
		tmpMarker = getNextInNode(marker);
		extractPassageMarker(marker);
		transposePassageMarker(marker, bypass);
		incrementFinishOffset(getTwinMarker(marker),
				      getNodeLength(bypass));
	}

}

static void adjustShortReads(Node * target, PassageMarker * pathMarker)
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

static Node *bypass()
{
	Node *bypass = getNode(path);
	Node *next = NULL;
	Arc *arc;
	Category cat;
	PassageMarker *nextMarker;

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
			// Update virtual coverage
			for (cat = 0; cat < CATEGORIES; cat++)
				incrementVirtualCoverage(bypass, cat,
							 getVirtualCoverage
							 (next, cat));

			// Update original virtual coverage
			for (cat = 0; cat < CATEGORIES; cat++)
				incrementOriginalVirtualCoverage(bypass,
								 cat,
								 getOriginalVirtualCoverage
								 (next,
								  cat));
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
	PassageMarker *marker, *next;

	printf("Trimming read tips\n");

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (getUniqueness(node))
			continue;

		for (marker = getMarker(node); marker != NULL;
		     marker = next) {
			next = getNextInNode(marker);

			if (!isInitial(marker) && !isTerminal(marker))
				continue;

			if (isTerminal(marker))
				marker = getTwinMarker(marker);

			while (!getUniqueness(getNode(marker))) {
				if (next != NULL
				    && (marker == next
					|| marker == getTwinMarker(next)))
					next = getNextInNode(next);
				if (getNextInSequence(marker) != NULL) {
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

	puts("Read coherency...");
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

	printf("Confronted to %i multiple hits and %i null over %i\n",
	       multCounter, nullCounter, dbgCounter);

	puts("Read coherency over!");
}

void setMultiplicityCutoff(int value)
{
	if (value < 0) {
		printf("Negative long read multiplicity cutoff %i!\n",
		       value);
		puts("Exiting...");
		exit(1);
	}
	MULTIPLICITY_CUTOFF = value;
}
