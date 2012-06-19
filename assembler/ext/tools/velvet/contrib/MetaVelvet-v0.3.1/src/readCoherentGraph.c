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
// Original
#include <string.h>
#include <math.h>
// Original

#include "globals.h"
#include "graph.h"
#include "recycleBin.h"
#include "passageMarker.h"
#include "graphStats.h"
#include "concatenatedGraph.h"
#include "readSet.h"
// Original
#include "readCoherentGraph.h"
#include "shortReadPairs.h"
// Original

#define LONG_NODE_CUTOFF 50
#define LN2 0.693147
#define PROBABILITY_CUTOFF 5
#define MAX_READ_COUNT 100
#define MAX_READ_LENGTH 2000
// Original
#define LEN_ARRAY_OUTIN 2
// Original

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

// Original
static void getStatusOfUniqueness(Node * node, double probability);
// Original

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
		    -1 * LN2 / 2 +
		    nodeLength / (2 * expected_coverage) *
		    (expected_coverage * expected_coverage -
		     nodeDensity * nodeDensity / 2);
 
		// Original
		getStatusOfUniqueness(node, probability);
		// Original
			
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
	// Original
	IDnum outputStartingNodeID = 0;
	IDnum outputDestinationNodeID = 0;
	// Original

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

	// Original
	// Reserving startingNode and destinationNode
	outputStartingNodeID = getNodeID(startingNode);
	outputDestinationNodeID = getNodeID(destination);
	// Original

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

	// Original
	if (outputStartingNodeID != 0 && outputDestinationNodeID != 0) {
                printf("RBConnection\tStarting : %d\tDestination : %d\n", 
		       outputStartingNodeID, outputDestinationNodeID);
	}
	// Original

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

// Original
static void getStatusOfUniqueness(Node * node, double probability)
{
	Node *twin = getTwinNode(node);
	char *repeatByProb = "Unique", *repeatByArcCount = "Unique";
	IDnum index;
	Node *arrayNode[] = {node, twin};
	Arc *tmpArc = NULL;
	Node *tmpNode = NULL;
	char *arcDirection = "";
	
	if (probability > PROBABILITY_CUTOFF) {
		repeatByProb = "Unique";
	} else {
		repeatByProb = "Repeat";
	}
	if ( (simpleArcCount(node) >= 2 && simpleArcCount(twin) >= 1)
	     || (simpleArcCount(twin) >= 2 && simpleArcCount(node) >= 1) )
		repeatByArcCount = "Repeat";
	else
		repeatByArcCount = "Unique";
	
	printf("Node : %d\tLen : %d\tCov : %6.2f\t",
	       getNodeID(node), getNodeLength(node), getNodeDensity(node));
	printf("F : %6.2f\tbyF : %s\tbyArc : %s\n",
	       probability, repeatByProb, repeatByArcCount);
	
	if ( !(strcmp(repeatByProb, "Repeat") == 0 || 
	       strcmp(repeatByArcCount, "Repeat") == 0) )
		return;

	for (index = 0; index <= 1; index++) {
		if (index == 0)
			arcDirection = "Out";
		else
			arcDirection = "In";
			
		tmpArc = getArc(arrayNode[index]);
		while (tmpArc != NULL) {
			tmpNode = getDestination(tmpArc);
			
			printf("\t%s\tNode : %d\t", 
			       arcDirection, getNodeID(tmpNode));
			printf("Len : %d\tCov : %6.2f\n", 
			       getNodeLength(tmpNode), getNodeDensity(node));
			
			tmpArc = getNextArc(tmpArc);
		}
	}
}

boolean isUniqueSolexaSubgraph(Node * node, double expCovSubgraph)
{
	int nodeLength = getNodeLength(node);
	double nodeDensity, probability;

	if (nodeLength == 0) {
		return false;
	}
	if (nodeLength > LONG_NODE_CUTOFF) {
		nodeDensity = getNodeDensity(node);

		probability =
		    -1 * LN2 / 2 +
		    nodeLength / (2 * expCovSubgraph) *
		    (expCovSubgraph * expCovSubgraph -
		     nodeDensity * nodeDensity / 2);

		return probability > PROBABILITY_CUTOFF;
	} else {
		return false;
	}
}

void identifyUniqueNodesSubgraph(Graph * graph, int * subgraphMask, 
				 boolean(*isUniqueSubgraph) (Node *, double), 
				 double expCovSubgraph)
{
        IDnum index;
        Node *node;

        for (index = 0; index < nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index + 1);
		if (node == NULL)
			continue;
		if (subgraphMask[index + 1 + nodeCount(graph)] == 2)
			setUniqueness(node, 
				      isUniqueSubgraph(node, expCovSubgraph));
		else
			setUniqueness(node, false);
	}
}

static double getNearestExpCov(double argCov, double * expCovMulti)
{
	int index;
	double difference, minDiff = 1000000.0, resultExpCov = 0.0;

	for (index = 0; index < 2; index++) {
		if (expCovMulti[index] == -1)
			break;
		difference = fabs(argCov - expCovMulti[index]);
		if (difference < minDiff) {
			minDiff = difference;
			resultExpCov = expCovMulti[index];
		}
	}
			
	return resultExpCov;
}

static boolean compareExpCovOutIn(Node * arrayOutInNode[][LEN_ARRAY_OUTIN], 
				  double * expCovMulti, 
				  double repeatNodeCov, double repeatNodeCovSD)
{
	int outIndex, inIndex, arrayIndex, sortIndex;
	double outExpCov, tmpExpCov, checkedOutExpCov[LEN_ARRAY_OUTIN];
	Node *recordedInNode[LEN_ARRAY_OUTIN];
	int checkIndex, checkedCount = 0;
	boolean flagChecked = false, flagOutInMatched = true;
	double arrayCovOutIn[LEN_ARRAY_OUTIN * 2], tmpCovOutIn;
	double aveCovOutInAll = 0.0, aveCovOutInPrimary = 0.0;
	
	// Initialize recordedInNode
	for (inIndex = 0; inIndex < LEN_ARRAY_OUTIN; inIndex++)
		recordedInNode[inIndex] = NULL;

	// Check Coverage of Repeat Node
	 // Prepare arrayCovOutIn
	arrayIndex = 0;
 	for (outIndex = 0; outIndex < 2; outIndex++) {
		for (inIndex = 0; inIndex < LEN_ARRAY_OUTIN; inIndex++) {
			tmpCovOutIn = getNodeDensity(arrayOutInNode[outIndex][inIndex]);
			arrayCovOutIn[arrayIndex++] = tmpCovOutIn;
			aveCovOutInAll += tmpCovOutIn;
		}
	}
	aveCovOutInAll /= 2.0;
	 // Sort
	for (sortIndex = 0; sortIndex < LEN_ARRAY_OUTIN * 2; sortIndex++) {
		for (arrayIndex = sortIndex + 1; arrayIndex < LEN_ARRAY_OUTIN * 2;
		     arrayIndex++) {
			if (arrayCovOutIn[sortIndex] < arrayCovOutIn[arrayIndex]) {
				tmpCovOutIn = arrayCovOutIn[arrayIndex];
				arrayCovOutIn[arrayIndex] = arrayCovOutIn[sortIndex];
				arrayCovOutIn[sortIndex] = tmpCovOutIn;
			}
		}
	}
	aveCovOutInPrimary = (arrayCovOutIn[0] + arrayCovOutIn[1]) / (double) 2;
	 // If Repeat Node Coverage is Too Small, "return false"
	//printf("CovRepeat = %f, aveCovAll = %f, aveCovPrimary = %f -> ",
	//       repeatNodeCov, aveCovOutInAll, aveCovOutInPrimary);
	if ((aveCovOutInAll * (1.0 + repeatNodeCovSD) < repeatNodeCov) 
	    || (aveCovOutInAll * (1.0 - repeatNodeCovSD) > repeatNodeCov)) {
		//printf("Dubious\n");
		return false;
	}
	//printf("Valid\n");


 	for (outIndex = 0; outIndex < LEN_ARRAY_OUTIN; outIndex++) {
		if (arrayOutInNode[0][outIndex] == NULL)
			break;

		if (!flagOutInMatched)
			return false;
		
		// Judge whether this OutCov has been checked or not
		flagChecked = false;
		tmpExpCov = getNodeDensity(arrayOutInNode[0][outIndex]);
		outExpCov = getNearestExpCov(tmpExpCov, expCovMulti);
		for (checkIndex = 0; checkIndex < checkedCount; checkIndex++) {
			if (outExpCov == checkedOutExpCov[checkIndex]) {
				flagChecked = true;
				break;
			}
		}
		if (flagChecked)
			continue;
		checkedOutExpCov[checkedCount++] = outExpCov;

		// Compare exp_cov between Out- and In-Nodes
		flagOutInMatched = false;
		for (inIndex = 0; inIndex < LEN_ARRAY_OUTIN; inIndex++) {
			if (arrayOutInNode[1][inIndex] == NULL)
				break;
			
			// Record corresponding In-Node
			tmpExpCov = getNodeDensity(arrayOutInNode[1][inIndex]);
			if (outExpCov == getNearestExpCov(tmpExpCov, expCovMulti)) {
				recordedInNode[outIndex] = arrayOutInNode[1][inIndex];
				flagOutInMatched = true;
				break;
			}
		}
	}

	for (inIndex = 0; inIndex < LEN_ARRAY_OUTIN; inIndex++) {
		if (arrayOutInNode[1][inIndex] == NULL)
			break;
		arrayOutInNode[1][inIndex] = recordedInNode[inIndex];
	}

	if (flagOutInMatched && checkedCount > 1)
		return true;
	else
		return false;
}

static boolean isInterRepeat(Node * node, double * expCovMulti,
			     Node * arrayOutInNode[2][LEN_ARRAY_OUTIN],
			     double repeatNodeCovSD)
{
	Node *twin = getTwinNode(node);
	char *isRepeat = "Unique";
	IDnum nodeIndex, arcIndex;
	Node *arrayNode[] = {node, twin};
	Arc *tmpArc = NULL;
	Node *tmpNode = NULL;
	char *arcDirection = "";
	boolean flagOutputNodeInfo = true;
	boolean resultCompareECOI = false;
	
	// Judge Unique or Repeat by ArcCount
	if (simpleArcCount(node) == 2 && simpleArcCount(twin) == 2)
		isRepeat = "Repeat";
	else
		isRepeat = "Unique";
	
	// Not InterRepeat but Unique
	if (strcmp(isRepeat, "Unique") == 0)
		return false;

	// Record Out- and In-Nodes
	for (nodeIndex = 0; nodeIndex <= 1; nodeIndex++) {
		arcIndex = 0;
		tmpArc = getArc(arrayNode[nodeIndex]);
		while (tmpArc != NULL) {
			tmpNode = getDestination(tmpArc);
			
			// Record to array
			if (nodeIndex == 0) {
				arrayOutInNode[nodeIndex][arcIndex++] 
					= tmpNode;
			}
			else {
				arrayOutInNode[nodeIndex][arcIndex++]
					= getTwinNode(tmpNode);
			}
			
			tmpArc = getNextArc(tmpArc);
		}
	}

	// Judge Intra- or Inter-Repeat by Out- and In-Nodes
	resultCompareECOI = compareExpCovOutIn(arrayOutInNode, expCovMulti,
					       getNodeDensity(node), repeatNodeCovSD);

	if (resultCompareECOI) {
		// Output Out- and In-Nodes
		for (nodeIndex = 0; nodeIndex <= 1; nodeIndex++) {
			arcIndex = 0;
			tmpArc = getArc(arrayNode[nodeIndex]);
			while (tmpArc != NULL) {
				tmpNode = getDestination(tmpArc);
				
				if (nodeIndex == 0)
					arcDirection = "Out";
				else
					arcDirection = "In";
				
				// Node Information
				if (flagOutputNodeInfo) {
					printf("Node : %d\tLen : %d \t",
					       getNodeID(node), getNodeLength(node));
					printf("Cov : %6.2f\tbyArc : %s\n", 
					       getNodeDensity(node), isRepeat);
					flagOutputNodeInfo = false;
				}
				// Connecting Information 
				printf("\t%s\tNode : %d\t", 
				       arcDirection, getNodeID(tmpNode));
				printf("Len : %d \tCov : %6.2f\n", 
				       getNodeLength(tmpNode), 
				       getNodeDensity(tmpNode));
				
				tmpArc = getNextArc(tmpArc);
			}
		}
	}

	return resultCompareECOI;
}

int identifyAndSeparateInterRepeats(Graph * argGraph, double * expCovMulti,
				    double repeatNodeCovSD)
{
	int graphIndex, nodeIndex, arcIndex;
	Graph *graph = argGraph;
	Node *node, *outNode, *inNode;
	Node *arrayOutInNode[2][LEN_ARRAY_OUTIN];
	int numInterRepeat = 0;

	puts("\nIdentifying and Separating InterRepeats");
	
	// Reset NodeStatus and Uniqueness
	resetNodeStatus(graph); resetUniqueness(graph);

	for (graphIndex = 0; graphIndex < nodeCount(graph); graphIndex++) {
		node = getNodeInGraph(graph, graphIndex + 1);
		
		if (getNodeID(node) == 0)
			continue;

		// Initialize arrayOutIn and arrayOutInNode
		for (arcIndex = 0; arcIndex < LEN_ARRAY_OUTIN; arcIndex++) {
			for (nodeIndex = 0; nodeIndex <= 1; nodeIndex++) {
				arrayOutInNode[nodeIndex][arcIndex] = NULL;
			}
		}
				
		// Identify InterRepeats
		if (!isInterRepeat(node, expCovMulti, arrayOutInNode, repeatNodeCovSD))
			continue;
		printf("Identified InterRepeat Node %d\n", getNodeID(node));
		numInterRepeat++;
		
		// Separate the InterRepeat
		for (arcIndex = 0; arcIndex < LEN_ARRAY_OUTIN; arcIndex++) {
			if (arrayOutInNode[1][arcIndex] == NULL)
				break;

			inNode = arrayOutInNode[1][arcIndex];
			outNode = arrayOutInNode[0][arcIndex];

			setNodeStatus(node, true); setUniqueness(node, false);
			setNodeStatus(inNode, true); setUniqueness(inNode, true);
			setNodeStatus(outNode, true); setUniqueness(outNode, true);
		       			
			if (!pushNeighboursInterRepeat(inNode, node, outNode, graph)) {
				printf("Error!! Separating Failed at Node %d",
					getNodeID(node));
				printf(" -- In : %d Out : %d\n\n", 
					getNodeID(inNode), getNodeID(outNode));
				exit(1);
			}
		}
	}

	// Reset NodeStatus and Uniqueness
	resetNodeStatus(graph); resetUniqueness(graph);

	// Return the number of InterRepeats
	return numInterRepeat;
}

static boolean trimLongReadTipsSubgraph(int * subgraphMask)
{
	IDnum index;
	Node *node;
	PassageMarker *marker, *next;
	boolean flagLongRead = false;

	//printf("Trimming read tips\n");

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (node == NULL || subgraphMask[index + nodeCount(graph)] != 1)
			continue;

		if (getUniqueness(node))
			continue;

		for (marker = getMarker(node); marker != NULL; marker = next) {
			next = getNextInNode(marker);

			flagLongRead = true;

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

	return flagLongRead;
}

static int computeNodeCount(Graph * argGraph)
{
	int index;
	int count = 0;
	Node *node;
	Graph *graph = argGraph;

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node != NULL)
			count++;
	}

	return count;
}

void readCoherentSubgraph(Graph * inGraph, double expCovSubgraph, 
			  ReadSet * reads, int * subgraphMask)
{
	IDnum nodeIndex;
	Node *node;
	IDnum previousNodeCount = 0;
	int checkModified = -1;
	
	graph = inGraph;
	listMemory = newRecycleBin(sizeof(PassageMarkerList), 100000);
	expected_coverage = expCovSubgraph;
	sequences = reads->tSequences;
	
	//puts("Read coherency...");

	if (!trimLongReadTipsSubgraph(subgraphMask)) {
		destroyRecycleBin(listMemory);
		//puts("Read Coherency didn't work. No Long Reads in the Subgraph");
		return;
	}

	while (previousNodeCount != computeNodeCount(graph)) {

		previousNodeCount = computeNodeCount(graph);
			
		for (nodeIndex = 1; nodeIndex <= nodeCount(graph); nodeIndex++) {
			node = getNodeInGraph(graph, nodeIndex);
			if (node == NULL || !getUniqueness(node))
				continue;

			while (uniqueNodesConnect(node))
				node = bypass();
			
			node = getTwinNode(node);
			while (uniqueNodesConnect(node))
				node = bypass();
		}

		checkModified++;
	}
	
	destroyRecycleBin(listMemory);

	//printf("readCoherentSubgraph checkModified = %d\n", checkModified);
 
	//puts("Read coherency over!");
}
// Original
