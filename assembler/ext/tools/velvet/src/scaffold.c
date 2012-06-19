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
#include <math.h>
#include <sys/time.h>
 
#ifdef _OPENMP
#include <omp.h>
#endif

#include "globals.h"
#include "graph.h"
#include "concatenatedGraph.h"
#include "recycleBin.h"
#include "locallyCorrectedGraph.h"
#include "passageMarker.h"
#include "readSet.h"
#include "utility.h"
#include "scaffold.h"

#define BLOCK_SIZE  100000
#define LN2 1.4

static int PEBBLE_ROUND_NUM = 0;

typedef struct readOccurence_st ReadOccurence;
static double paired_exp_fraction = 0.1;

struct connection_st {
	Node *destination;
	Connection *right;
	Connection *left;
	Connection *twin;
	float distance;
	float variance;
	IDnum direct_count;
	IDnum paired_count;
	unsigned char clean;
}  ATTRIBUTE_PACKED;

struct readOccurence_st {
	IDnum position;
	IDnum offset;
	IDnum nodeID;
}  ATTRIBUTE_PACKED;

// Global params
static IDnum UNRELIABLE_CONNECTION_CUTOFF = 5;

// Global pointers
static Graph *graph;
static Connection **scaffold = NULL;
static RecycleBin *connectionMemory = NULL;
static boolean estimated[CATEGORIES + 1];

#ifdef _OPENMP

#define READS_PER_LOCK 32


/* Array of reads locks */
static omp_lock_t *readsLocks = NULL;
/* Array of per-node locks */
static omp_lock_t *nodeLocks = NULL;

static void
createReadsLocks()
{
	Coordinate nbLocks;
	Coordinate lockIndex;

	if (readsLocks)
		free (readsLocks);
	nbLocks = 1 + sequenceCount(graph) / READS_PER_LOCK;
	readsLocks = mallocOrExit(nbLocks, omp_lock_t);

	#pragma omp parallel for
	for (lockIndex = 0; lockIndex < nbLocks; lockIndex++)
		omp_init_lock(readsLocks + lockIndex);
}

static inline void lockRead(IDnum readID)
{
	omp_set_lock (readsLocks + readID / READS_PER_LOCK);
}

static inline void unLockRead(IDnum readID)
{
	omp_unset_lock (readsLocks + readID / READS_PER_LOCK);
}

static void
createNodeLocks(Graph *graph)
{
	IDnum nbNodes;
	IDnum nodeIndex;

	nbNodes = nodeCount(graph) + 1;
	if (nodeLocks)
		free (nodeLocks);
	nodeLocks = mallocOrExit(nbNodes, omp_lock_t);

	#pragma omp parallel for
	for (nodeIndex = 0; nodeIndex < nbNodes; nodeIndex++)
		omp_init_lock(nodeLocks + nodeIndex);
}

/* Tries to avoid deadlocking */
static inline void lockTwoNodes(IDnum nodeID, IDnum node2ID)
{
	if (nodeID < 0)
		nodeID = -nodeID;
	if (node2ID < 0)
		node2ID = -node2ID;

	/* Lock lowest ID first to avoid deadlocks */
	if (nodeID < node2ID)
	{
		omp_set_lock (nodeLocks + nodeID);
		omp_set_lock (nodeLocks + node2ID);
	}
	else
	{
		omp_set_lock (nodeLocks + node2ID);
		omp_set_lock (nodeLocks + nodeID);
	}
}

static inline void unLockTwoNodes(IDnum nodeID, IDnum node2ID)
{
	if (nodeID < 0)
		nodeID = -nodeID;
	if (node2ID < 0)
		node2ID = -node2ID;

	omp_unset_lock (nodeLocks + nodeID);
	omp_unset_lock (nodeLocks + node2ID);
}
#endif

static Connection *allocateConnection()
{
	Connection *connect;
#ifdef _OPENMP
#pragma omp critical
	{
#endif
	if (connectionMemory == NULL)
		connectionMemory =
		    newRecycleBin(sizeof(Connection), BLOCK_SIZE);

	connect = allocatePointer(connectionMemory);
#ifdef _OPENMP
	}
#endif
	connect->destination = NULL;
	connect->clean = false;
	return connect;
}

static void deallocateConnection(Connection * connect)
{
	deallocatePointer(connectionMemory, connect);
}

Node * getConnectionDestination(Connection * connect) {
	return connect->destination;
}

Connection * getNextConnection(Connection * connect) {
	return connect->right;
}

Connection * getTwinConnection(Connection * connect) {
	return connect->twin;
}

Coordinate getConnectionDistance(Connection * connect) {
	return (Coordinate) connect->distance;
}

double getConnectionVariance(Connection * connect) {
	return connect->variance;
}

IDnum getConnectionDirectCount(Connection * connect) {
	return connect->direct_count;
}

IDnum getConnectionPairedCount(Connection * connect) {
	return connect->paired_count;
}

Connection * getConnection(Node * node) {
	return scaffold[getNodeID(node) + nodeCount(graph)];
}

void incrementConnectionDistance(Connection * connect, Coordinate increment) {
	connect->distance += increment;
}

static double norm(double X)
{
	return 0.4 * exp(-X * X / 2);
}

static double normInt(double X, double Y)
{
	return (erf(0.7 * Y) - erf(0.7 * X)) / 2;
}

static IDnum expectedNumberOfConnections(IDnum IDA, Connection * connect,
					 IDnum ** counts, Category cat)
{
	Node *A = getNodeInGraph(graph, IDA);
	Node *B = connect->destination;
	double left, middle, right;
	Coordinate longLength, shortLength, D;
	IDnum longCount;
	double M, N, O, P;
	Coordinate mu = getInsertLength(graph, cat);
	double sigma = sqrt(getInsertLength_var(graph, cat));
	double result;

	if (mu <= 0)
		return 0;

	if (getNodeLength(A) < getNodeLength(B)) {
		longLength = getNodeLength(B);
		shortLength = getNodeLength(A);
		longCount = counts[cat][getNodeID(B) + nodeCount(graph)];
	} else {
		longLength = getNodeLength(A);
		shortLength = getNodeLength(B);
		longCount = counts[cat][IDA + nodeCount(graph)];
	}

	D = getConnectionDistance(connect) - (longLength + shortLength) / 2;

	M = (D - mu) / sigma;
	N = (D + shortLength - mu) / sigma;
	O = (D + longLength - mu) / sigma;
	P = (D + shortLength + longLength - mu) / sigma;

	left = ((norm(M) - norm(N)) - M * normInt(M, N)) * sigma;
	middle = shortLength * normInt(N, O);
	right = ((norm(O) - norm(P)) - P * normInt(O, P)) * (-sigma);

	result = (longCount * (left + middle + right)) / longLength;

	if (result > 0)
		return (IDnum) result;
	else
		return 0;
}

void destroyConnection(Connection * connect, IDnum nodeID)
{
	Connection *previous, *next;

	//velvetLog("Destroying connection from %li to %li\n", nodeID, getNodeID(connect->destination));

	if (connect == NULL)
		return;

	previous = connect->left;
	next = connect->right;

	if (previous != NULL)
		previous->right = next;
	if (next != NULL)
		next->left = previous;

	if (scaffold[nodeID + nodeCount(graph)] == connect)
		scaffold[nodeID + nodeCount(graph)] = next;

	if (connect->twin != NULL) {
		connect->twin->twin = NULL;
		destroyConnection(connect->twin,
				  getNodeID(connect->destination));
	}

	deallocateConnection(connect);
}

static boolean testConnection(IDnum IDA,
			      Connection *connect,
			      IDnum **counts,
			      boolean *shadows)
{
	IDnum total = 0;
	Category cat;

	// Spare unique -> undetermined node connections
	if (!getUniqueness(connect->destination))
		return true;

	// Destroy tenuous connections
	if (connect->paired_count + connect->direct_count <
	    UNRELIABLE_CONNECTION_CUTOFF)
		return false;

	for (cat = 0; cat < CATEGORIES; cat++)
		if (!shadows[cat] || cat <= PEBBLE_ROUND_NUM)
			total += expectedNumberOfConnections(IDA, connect, counts, cat);

	// Remove inconsistent connections
	return connect->paired_count >= total * paired_exp_fraction;
}

static IDnum *computeReadToNodeCounts(Coordinate *totalCount)
{
	IDnum nodeIndex;
	IDnum maxNodeIndex = 2 * nodeCount(graph) + 1;
	IDnum maxReadIndex = sequenceCount(graph) + 1;
	IDnum *readNodeCounts = callocOrExit(maxReadIndex, IDnum);
	unsigned char *readMarker = callocOrExit(1 + maxReadIndex / 8, unsigned char);
	Coordinate total = 0;

	velvetLog("Computing read to node mapping array sizes\n");

#ifdef _OPENMP
	#pragma omp parallel for reduction(+:total)
#endif
	for (nodeIndex = 0; nodeIndex < maxNodeIndex; nodeIndex++) {
		Node *node;
		ShortReadMarker *nodeArray;
		IDnum nodeReadCount;
		IDnum readIndex;

		node = getNodeInGraph(graph, nodeIndex - nodeCount(graph));
		if (node == NULL)
			continue;
		nodeArray = getNodeReads(node, graph);
		nodeReadCount = getNodeReadCount(node, graph);

		// Short reads
		for (readIndex = 0; readIndex < nodeReadCount; readIndex++) {
			ShortReadMarker *shortMarker;
			IDnum readID;

			shortMarker = getShortReadMarkerAtIndex(nodeArray,
								readIndex);
			readID = getShortReadMarkerID(shortMarker);
#ifdef _OPENMP
			#pragma omp atomic
#endif
			readNodeCounts[readID]++;
			total++;
		}
	}

	for (nodeIndex = 0; nodeIndex < maxNodeIndex; nodeIndex++) {
		Node *node;
		PassageMarkerI marker;

		node = getNodeInGraph(graph, nodeIndex - nodeCount(graph));
		if (node == NULL)
			continue;
		// Long reads
		for (marker = getMarker(node); marker != NULL_IDX;
		     marker = getNextInNode(marker)) {
			IDnum readIndex = getPassageMarkerSequenceID(marker);;

			if (readIndex < 0)
				continue;

			const unsigned int idx = readIndex / 8;
			const unsigned int mask = 1 << (readIndex & 7);
			if (readMarker[idx] & mask)
				continue;

			readNodeCounts[readIndex]++;
			total++;
			readMarker[idx] |= mask;
		}

		// Clean up marker array
		for (marker = getMarker(node); marker != NULL_IDX;
		     marker = getNextInNode(marker)) {
			IDnum readIndex = getPassageMarkerSequenceID(marker);
			if (readIndex > 0)
				// No need to go bit-wise
				readMarker[readIndex / 8] = 0;
		}
	}

	*totalCount = total;
	free(readMarker);
	return readNodeCounts;
}

static ReadOccurence **allocateReadToNodeTables(IDnum * readNodeCounts,
						Coordinate totalCount,
						ReadOccurence **readNodesArray)
{
	Coordinate offset = 0;
	IDnum readIndex;
	IDnum maxReadIndex = sequenceCount(graph) + 1;
	ReadOccurence **readNodes = callocOrExit(maxReadIndex, ReadOccurence *);
	*readNodesArray = callocOrExit(totalCount, ReadOccurence);

	for (readIndex = 1; readIndex < maxReadIndex; readIndex++) {
		if (readNodeCounts[readIndex] != 0) {
			readNodes[readIndex] = *readNodesArray + offset;
			offset += readNodeCounts[readIndex];
			readNodeCounts[readIndex] = 0;
		}
	}

	return readNodes;
}

static void computePartialReadToNodeMappingShort(IDnum nodeID,
						 ReadOccurence ** readNodes,
						 IDnum * readNodeCounts)
{
	ShortReadMarker *shortMarker;
	IDnum index, readIndex;
	ReadOccurence *readArray, *readOccurence;
	Node *node = getNodeInGraph(graph, nodeID);
	ShortReadMarker *nodeArray = getNodeReads(node, graph);
	IDnum nodeReadCount = getNodeReadCount(node, graph);

	for (index = 0; index < nodeReadCount; index++) {
		shortMarker = getShortReadMarkerAtIndex(nodeArray, index);
		readIndex = getShortReadMarkerID(shortMarker);
		readArray = readNodes[readIndex];
#ifdef _OPENMP
		lockRead(readIndex);
#endif
		readOccurence = &readArray[readNodeCounts[readIndex]];
		readOccurence->nodeID = nodeID;
		readOccurence->position =
		    getShortReadMarkerPosition(shortMarker);
		readOccurence->offset =
		    getShortReadMarkerOffset(shortMarker);
		readNodeCounts[readIndex]++;
#ifdef _OPENMP
		unLockRead(readIndex);
#endif
	}
}

static void computePartialReadToNodeMappingLong(IDnum nodeID,
						ReadOccurence ** readNodes,
						IDnum * readNodeCounts,
						unsigned char *readMarker,
						ReadSet * reads)
{
	IDnum readIndex;
	ReadOccurence *readArray, *readOccurence;
	Node *node = getNodeInGraph(graph, nodeID);
	PassageMarkerI marker;

	for (marker = getMarker(node); marker != NULL_IDX;
	     marker = getNextInNode(marker)) {
		readIndex = getPassageMarkerSequenceID(marker);
		if (readIndex <= 0 || reads->categories[readIndex - 1] == REFERENCE)
			continue;

		const unsigned int idx = readIndex / 8;
		const unsigned int mask = 1 << (readIndex & 7);
		if (readMarker[idx] & mask) {
			readArray = readNodes[readIndex];
			readOccurence =
			    &readArray[readNodeCounts[readIndex] - 1];
			readOccurence->position = -1;
			readOccurence->offset = -1;
		} else {
			readArray = readNodes[readIndex];
			readOccurence =
			    &readArray[readNodeCounts[readIndex]];
			readOccurence->nodeID = nodeID;
			readOccurence->position = getStartOffset(marker);
			readOccurence->offset =
			    getPassageMarkerStart(marker);
			readNodeCounts[readIndex]++;
			readMarker[idx] |= mask;
		}
	}

	for (marker = getMarker(node); marker != NULL_IDX;
	     marker = getNextInNode(marker)) {
		readIndex = getPassageMarkerSequenceID(marker);
		if (readIndex > 0)
			// No need to go bit-wise
			readMarker[readIndex / 8] = 0;
	}
}

static ReadOccurence **computeReadToNodeMappings(IDnum * readNodeCounts,
						 ReadSet * reads,
						 Coordinate totalCount,
						 ReadOccurence **readNodesArray)
{
	unsigned char *readMarker;
	IDnum nodeID;
	IDnum nodes = nodeCount(graph);
	ReadOccurence **readNodes = allocateReadToNodeTables(readNodeCounts,
							     totalCount,
							     readNodesArray);

	velvetLog("Computing read to node mappings\n");

#ifdef _OPENMP
	createReadsLocks();
	#pragma omp parallel for
#endif
	for (nodeID = -nodes; nodeID <= nodes; nodeID++)
		if (nodeID != 0 && getNodeInGraph(graph, nodeID))
			computePartialReadToNodeMappingShort(nodeID, readNodes,
							     readNodeCounts);

#ifdef _OPENMP
	free(readsLocks);
	readsLocks = NULL;
#endif

	readMarker = callocOrExit(1 + sequenceCount(graph) / 8, unsigned char);
	for (nodeID = -nodes; nodeID <= nodes; nodeID++)
		if (nodeID != 0 && getNodeInGraph(graph, nodeID))
			computePartialReadToNodeMappingLong(nodeID, readNodes,
							    readNodeCounts,
							    readMarker,
							    reads);

	free(readMarker);
	return readNodes;
}

static unsigned char * countCoOccurences(IDnum * coOccurencesCount,
					 ReadOccurence ** readNodes,
					 IDnum * readNodeCounts,
					 IDnum * readPairs,
					 Category * cats)
{
	IDnum readIndex, readPairIndex;
	IDnum readNodeCount;
	IDnum readOccurenceIndex, readPairOccurenceIndex;
	ReadOccurence * readOccurence, *readPairOccurence;
	unsigned char *interestingReads = callocOrExit(1 + sequenceCount(graph) / 8, unsigned char);
	Category libID;

	for (libID = 0; libID < CATEGORIES + 1; libID++)
		coOccurencesCount[libID] = 0;

	for (readIndex = 0; readIndex < sequenceCount(graph); readIndex++) {
		// Eliminating dodgy, unpaired, already counted or user-specified reads
		if ( readPairs[readIndex] < readIndex
		    || getInsertLength(graph, cats[readIndex]) > -1)
			continue;

		// Check for co-occurence
		// We know that for each read the read occurences are ordered by increasing node ID
		// Therefore one list is followed by increasing index, whereas the other is followed 
		// by decreasing index
		libID = cats[readIndex] / 2;
		readPairIndex = readPairs[readIndex];
		
		readOccurenceIndex = 0;
		readOccurence = readNodes[readIndex + 1];
		readNodeCount = readNodeCounts[readIndex + 1];

		readPairOccurenceIndex = readNodeCounts[readPairIndex + 1] - 1;
		readPairOccurence = &(readNodes[readPairIndex + 1][readPairOccurenceIndex]);

		while (readOccurenceIndex < readNodeCount && readPairOccurenceIndex >= 0) {
			if (readOccurence->nodeID == -readPairOccurence->nodeID) {
				if (readOccurence->position > 0 && readPairOccurence->position > 0) {
					coOccurencesCount[libID]++;
					interestingReads[readIndex / 8] |= 1 << (readIndex & 7);
					break;
				} else {
					readOccurence++;
					readOccurenceIndex++;	
					readPairOccurence--;
					readPairOccurenceIndex--;	
				}
			} else if (readOccurence->nodeID < -readPairOccurence->nodeID) {
				readOccurence++;
				readOccurenceIndex++;	
			} else {
				readPairOccurence--;
				readPairOccurenceIndex--;	
			}
		}
	}

	return interestingReads;
}

static void measureCoOccurences(IDnum ** coOccurences,
				unsigned char * interestingReads,
				ReadOccurence ** readNodes,
				IDnum * readNodeCounts,
				IDnum * readPairs,
				Category * cats)
{
	IDnum coOccurencesIndex[CATEGORIES + 1];
	IDnum observationIndex;
	IDnum readIndex, readPairIndex;
	IDnum readNodeCount;
	IDnum readOccurenceIndex, readPairOccurenceIndex;
	ReadOccurence * readOccurence, *readPairOccurence;
	Category libID;

	for (libID = 0; libID < CATEGORIES + 1; libID++)
		coOccurencesIndex[libID] = 0;

	for (readIndex = 0; readIndex < sequenceCount(graph); readIndex++) {
		// Eliminating dodgy, unpaired, already counted or user-specified reads
		if (!(interestingReads[readIndex / 8] & (1 << (readIndex & 7))))
			continue;
		
		// Find co-occurence
		// We know that for each read the read occurences are ordered by increasing node ID
		libID = cats[readIndex]/2;
		readPairIndex = readPairs[readIndex];	
		observationIndex = coOccurencesIndex[libID];
		
		readOccurence = readNodes[readIndex + 1];
		readOccurenceIndex = 0;
		readNodeCount = readNodeCounts[readIndex + 1];

		readPairOccurenceIndex = readNodeCounts[readPairIndex + 1] - 1;
		readPairOccurence = &(readNodes[readPairIndex + 1][readPairOccurenceIndex]);

		while (readOccurenceIndex < readNodeCount && readPairOccurenceIndex >= 0) {
			if (readOccurence->nodeID == -readPairOccurence->nodeID) {
				if (readOccurence->position > 0 && readPairOccurence->position > 0) {
					coOccurences[libID][observationIndex] = 
					      getNodeLength(getNodeInGraph(graph, readOccurence->nodeID))
					      + getWordLength(graph) - 1
					      - (readOccurence->position - readOccurence->offset)	
					      - (readPairOccurence->position - readPairOccurence->offset);
					coOccurencesIndex[libID]++;
					break;
				} else {
					readOccurence++;
					readOccurenceIndex++;	
					readPairOccurence--;
					readPairOccurenceIndex--;	
				}
			} else if (readOccurence->nodeID < -readPairOccurence->nodeID) {
				readOccurence++;
				readOccurenceIndex++;	
			} else {
				readPairOccurence--;
				readPairOccurenceIndex--;	
			}
		}
	}
}

int compareReadOccurences(const void *A, const void * B) {
	IDnum * cA = (IDnum *) A;
	IDnum * cB = (IDnum *) B;

	if (*cA > *cB)
		return 1;
	if (*cA == *cB)
		return 0;
	return -1;	
}

static void estimateLibraryInsertLength(IDnum * coOccurences, IDnum coOccurencesCount, Category libID) {
	Coordinate median, variance;
	IDnum index;
	int counter = 0;
	qsort(coOccurences, coOccurencesCount, sizeof(IDnum), compareReadOccurences);

	median = coOccurences[coOccurencesCount / 2];

	// Modified variance around the median (proxy for expected value) 
	// interval censoring
	variance = 0;
	for (index = 0; index < coOccurencesCount; index++) {
		if (coOccurences[index] > 0 && coOccurences[index] < 5 * median) {
			variance += (coOccurences[index] - median) * (coOccurences[index] - median);
			counter++;
		}
	}
	if (counter) 
		variance /= counter;
	else {
		variance = 0;
		for (index = 0; index < coOccurencesCount; index++)
			variance += (coOccurences[index] - median) * (coOccurences[index] - median);
		variance /= coOccurencesCount;
	}
	
	// To avoid subsequent divisions by zero
	if (variance == 0)
		variance = 1;

	velvetLog("Paired-end library %i has length: %lli, sample standard deviation: %lli\n", libID + 1, (long long) median, (long long) sqrt(variance));
	setInsertLengths(graph, libID, median, sqrt(variance));
	estimated[libID] = true;
}

static void estimateLibraryInsertLengths(IDnum ** coOccurences, IDnum * coOccurencesCounts) {
	Category libID;

	for (libID = 0; libID < CATEGORIES + 1; libID++)
		estimated[libID] = false;

	for (libID = 0; libID < CATEGORIES + 1; libID++)
		if (coOccurencesCounts[libID] > 0)
			estimateLibraryInsertLength(coOccurences[libID], coOccurencesCounts[libID], libID);
}

static void estimateMissingInsertLengths(ReadOccurence ** readNodes, IDnum * readNodeCounts, IDnum * readPairs, Category * cats) {
	IDnum * coOccurences[CATEGORIES + 1];
	IDnum coOccurencesCounts[CATEGORIES + 1]; 
	Category libID;

	velvetLog("Estimating library insert lengths...\n");

	unsigned char * interestingReads = countCoOccurences(coOccurencesCounts, readNodes, readNodeCounts, readPairs, cats);

	for (libID = 0; libID < CATEGORIES + 1; libID++)
		coOccurences[libID] = callocOrExit(coOccurencesCounts[libID], IDnum);

	measureCoOccurences(coOccurences, interestingReads, readNodes, readNodeCounts, readPairs, cats);
	estimateLibraryInsertLengths(coOccurences, coOccurencesCounts);

	for (libID = 0; libID < CATEGORIES + 1; libID++)
		free(coOccurences[libID]);
	
	free(interestingReads);

	velvetLog("Done\n");
}

static void createTwinConnection(IDnum nodeID, IDnum node2ID,
				 Connection * connect)
{
	Connection *newConnection = allocateConnection();
	IDnum nodeIndex = nodeID + nodeCount(graph);

	// Fill in
	newConnection->distance = connect->distance;
	newConnection->variance = connect->variance;
	newConnection->direct_count = connect->direct_count;
	newConnection->paired_count = connect->paired_count;
	newConnection->destination = getNodeInGraph(graph, node2ID);

	// Batch to twin
	newConnection->twin = connect;
	connect->twin = newConnection;

	// Insert in scaffold
	newConnection->left = NULL;
	newConnection->right = scaffold[nodeIndex];
	if (scaffold[nodeIndex] != NULL)
		scaffold[nodeIndex]->left = newConnection;
	scaffold[nodeIndex] = newConnection;
}

Connection *createNewConnection(IDnum nodeID, IDnum node2ID,
				       IDnum direct_count,
				       IDnum paired_count,
				       Coordinate distance,
				       double variance)
{
	Node *destination = getNodeInGraph(graph, node2ID);
	IDnum nodeIndex = nodeID + nodeCount(graph);
	Connection *connect = allocateConnection();

	// Fill in 
	connect->destination = destination;
	connect->direct_count = direct_count;
	connect->paired_count = paired_count;
	connect->distance = (double) distance;
	connect->variance = variance;

	// Insert in scaffold
	connect->left = NULL;
	connect->right = scaffold[nodeIndex];
	if (scaffold[nodeIndex] != NULL)
		scaffold[nodeIndex]->left = connect;
	scaffold[nodeIndex] = connect;

	// Event. pair up to twin
	if (getUniqueness(destination))
		createTwinConnection(node2ID, nodeID, connect);
	else
		connect->twin = NULL;

	return connect;
}

void readjustConnection(Connection * connect, Coordinate distance,
			       double variance, IDnum direct_count,
			       IDnum paired_count)
{
	connect->direct_count += direct_count;
	connect->paired_count += paired_count;

	connect->distance =
	    (variance * connect->distance +
	     distance * connect->variance) / (variance +
					      connect->variance);
	connect->variance =
	    (variance *
	     connect->variance) / (variance + connect->variance);

	if (connect->twin != NULL) {
		connect->twin->distance = connect->distance;
		connect->twin->variance = connect->variance;
		connect->twin->direct_count = connect->direct_count;
		connect->twin->paired_count = connect->paired_count;
	}
}

//////////////////////////////////////
// Splay tree function for Connections
//////////////////////////////////////

/* This function can be called only if K2 has a left child */
/* Perform a rotate between a node (K2) and its left child */
/* Update heights, then return new root */

static Connection *connectionSingleRotateWithLeft(Connection * K2)
{
	Connection *K1;

	K1 = K2->left;
	K2->left = K1->right;
	K1->right = K2;

	return K1;		/* New root */
}

/* This function can be called only if K1 has a right child */
/* Perform a rotate between a node (K1) and its right child */
/* Update heights, then return new root */

static Connection *connectionSingleRotateWithRight(Connection * K1)
{
	Connection *K2;

	K2 = K1->right;
	K1->right = K2->left;
	K2->left = K1;

	return K2;		/* New root */
}

/* Top-down splay procedure, */
/* not requiring destination to be in tree */

static Connection *splayConnection(Connection * T, IDnum nodeID)
{
	Connection Header;
	Connection *LeftTreeMax, *RightTreeMin;

	if (T == NULL)
		return NULL;

	Header.left = Header.right = NULL;
	LeftTreeMax = RightTreeMin = &Header;

	while (nodeID != getNodeID(T->destination))
	{
		if (nodeID < getNodeID(T->destination))
		{
			if (T->left == NULL)
				break;
			if (nodeID  < getNodeID(T->left->destination))
				T = connectionSingleRotateWithLeft(T);
			if (T->left == NULL)
				break;
			/* Link right */
			RightTreeMin->left = T;
			RightTreeMin = T;
			T = T->left;
		}
		else
		{
			if (T->right == NULL)
				break;
			if (nodeID > getNodeID(T->right->destination))
				T = connectionSingleRotateWithRight(T);
			if (T->right == NULL)
				break;
			/* Link left */
			LeftTreeMax->right = T;
			LeftTreeMax = T;
			T = T->right;
		}
	} /* while nodeID != T->destination */

	/* Reassemble */
	LeftTreeMax->right = T->left;
	RightTreeMin->left = T->right;
	T->left = Header.right;
	T->right = Header.left;

	return T;
}

static Connection* findOrCreateConnection(IDnum nodeID,
					  IDnum node2ID)
{
	Connection **T;
	Connection *newConnection;
	IDnum nodeIndex;

	nodeIndex = nodeID + nodeCount(graph);
	T = scaffold + nodeIndex;

	if (*T == NULL)
	{
		newConnection = allocateConnection();

		newConnection->left = NULL;
		newConnection->right = NULL;
		*T = newConnection;
	}
	else
	{
		IDnum destID;

		*T = splayConnection(*T, node2ID);
		destID = getNodeID((*T)->destination);
		if (destID == node2ID)
			newConnection = *T;
		else
		{
			newConnection = allocateConnection();
			if (node2ID < destID)
			{
				newConnection->left = (*T)->left;
				newConnection->right = *T;
				(*T)->left = NULL;
			}
			else if (node2ID > destID)
			{
				newConnection->right = (*T)->right;
				newConnection->left = *T;
				(*T)->right = NULL;
			}
			*T = newConnection;
		}
	}

	return newConnection;
}

static Connection* findConnection(IDnum nodeID,
				  IDnum node2ID)
{
	Connection **T;
	IDnum nodeIndex;

	nodeIndex = nodeID + nodeCount(graph);
	T = scaffold + nodeIndex;

	if (*T == NULL)
		return NULL;
	else
	{
		IDnum destID;

		*T = splayConnection(*T, node2ID);
		destID = getNodeID((*T)->destination);
		if (destID == node2ID)
			return *T;
	}
	return NULL;
}

RecycleBin *connectionStackMemory = NULL;

typedef struct ConnectionStack_st ConnectionStack;

struct ConnectionStack_st
{
	Connection *connection;
	ConnectionStack *next;
};

#ifdef _OPENMP
static void initConnectionStackMemory(void)
{
	int n = omp_get_max_threads();

	#pragma omp critical
	{
		if (connectionStackMemory == NULL)
			connectionStackMemory = newRecycleBinArray(n, sizeof(ConnectionStack), BLOCK_SIZE);
	}
}
#endif

static ConnectionStack *allocateConnectionStack(void)
{
#ifdef _OPENMP
#ifdef DEBUG
	if (connectionStackMemory == NULL)
	{
		velvetLog("The memory for connection stack seems uninitialised, "
			  "this is probably a bug, aborting.\n");
		abort();
	}
#endif
	return allocatePointer(getRecycleBinInArray(connectionStackMemory,
				omp_get_thread_num()));
#else
	if (connectionStackMemory == NULL)
		connectionStackMemory =
		    newRecycleBin(sizeof(ConnectionStack), BLOCK_SIZE);

	return allocatePointer(connectionStackMemory);
#endif
}

static void deallocateConnectionStack(ConnectionStack *stack)
{
#ifdef _OPENMP
	deallocatePointer(getRecycleBinInArray(connectionStackMemory,
					       omp_get_thread_num()),
			  stack);
#else
	deallocatePointer(connectionStackMemory, stack);
#endif
}

static void destroyConnectionStackMemory(void)
{
#ifdef _OPENMP
	destroyRecycleBinArray(connectionStackMemory);
#else
	destroyRecycleBin(connectionStackMemory);
#endif
	connectionStackMemory = NULL;
}

static void pushConnectionStack(ConnectionStack **stack, Connection *connection)
{
	ConnectionStack *newElement;

	newElement = allocateConnectionStack();
	newElement->connection = connection;
	newElement->next = *stack;
	*stack = newElement;
}

static Connection *popConnectionStack(ConnectionStack **stack)
{
	ConnectionStack *nextElement;
	Connection *connection;

	if (*stack == NULL)
		return NULL;

	nextElement = (*stack)->next;
	connection = (*stack)->connection;
	deallocateConnectionStack(*stack);
	*stack = nextElement;

	return connection;
}

static void splayToList(Connection **connection)
{
	ConnectionStack *stack = NULL;
	Connection *current;
	Connection *list = NULL;

	if (*connection == NULL)
		return;

	for (current = *connection; current != NULL; current = popConnectionStack(&stack))
	{
		Connection *right;
		Connection *left;

		right = current->right;
		if (right != NULL)
			pushConnectionStack(&stack, right);
		left = current->left;
		if (left != NULL)
			pushConnectionStack(&stack, left);
		if (list != NULL)
			list->left = current;
		current->right = list;
		list = current;
	}
	list->left = NULL;
	*connection = list;
}

static void setAllConnectionsClean(void)
{
	IDnum nodeID;
	IDnum nodes = nodeCount(graph);

#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (nodeID = 2 * nodes; nodeID >= 0; nodeID--)
	{
		ConnectionStack *stack = NULL;
		Connection **connect;
		Connection *current;

		connect = scaffold + nodeID;
		if (*connect == NULL)
			continue;

		for (current = *connect; current != NULL; current = popConnectionStack(&stack))
		{
			Connection *right;
			Connection *left;

			current->clean = true;
			right = current->right;
			if (right != NULL)
				pushConnectionStack(&stack, right);
			left = current->left;
			if (left != NULL)
				pushConnectionStack(&stack, left);
		}
	}
}

static void fillNewConnectionInTree(Connection *connect,
				    Node *destination,
				    IDnum direct_count,
				    IDnum paired_count,
				    Coordinate distance,
				    double variance)
{
	connect->destination = destination;
	connect->direct_count = direct_count;
	connect->paired_count = paired_count;
	connect->distance = (double)distance;
	connect->variance = variance;
}

static void readjustConnectionInTree(Connection *connect,
				     IDnum direct_count,
				     IDnum paired_count,
				     Coordinate distance,
				     double variance)
{
	connect->direct_count += direct_count;
	connect->paired_count += paired_count;
	connect->distance = (variance * connect->distance + distance * connect->variance) /
			    (variance + connect->variance);
	connect->variance = (variance * connect->variance) / (variance + connect->variance);

	if (connect->twin != NULL)
	{
		connect->twin->direct_count = connect->direct_count;
		connect->twin->paired_count = connect->paired_count;
		connect->twin->distance = connect->distance;
		connect->twin->variance = connect->variance;
	}
}

static void createTwinConnectionInTree(IDnum nodeID,
				       IDnum node2ID,
				       Connection *connect)
{
	Connection *newConnection;

	newConnection = findOrCreateConnection(nodeID, node2ID);
	if (newConnection->destination == NULL)
	{
		fillNewConnectionInTree(newConnection,
					getNodeInGraph(graph, node2ID),
					connect->direct_count,
					connect->paired_count,
					(Coordinate)connect->distance,
					connect->variance);
		// Batch to twin
		newConnection->twin = connect;
		connect->twin = newConnection;
	}
	else
		readjustConnectionInTree(newConnection,
					 connect->direct_count,
					 connect->paired_count,
					 (Coordinate)connect->distance,
					 connect->variance);
}

static void createConnection(IDnum nodeID,
			     IDnum node2ID,
			     IDnum direct_count,
			     IDnum paired_count,
			     Coordinate distance,
			     double variance)
{
	Connection *connect;

	if (getUniqueness(getNodeInGraph(graph, node2ID)) && node2ID < nodeID) {
		return;
	}	

#ifdef _OPENMP
	lockTwoNodes(nodeID, node2ID);
#endif
	connect = findOrCreateConnection(nodeID, node2ID);
	if (connect->destination == NULL)
	{
		Node *destination = getNodeInGraph(graph, node2ID);
		fillNewConnectionInTree(connect,
					destination,
					direct_count,
					paired_count,
					distance,
					variance);

		if (getUniqueness(destination))
			createTwinConnectionInTree(node2ID, nodeID, connect);
		else
			connect->twin = NULL;
	}
	else
		readjustConnectionInTree(connect,
					 direct_count,
					 paired_count,
					 distance,
					 variance);

#ifdef _OPENMP
	unLockTwoNodes(nodeID, node2ID);
#endif
}

static void projectFromSingleRead(Node * node,
				  ReadOccurence * readOccurence,
				  Coordinate position,
				  Coordinate offset, Coordinate length)
{
	Coordinate distance = 0;
	Node *target = getNodeInGraph(graph, -readOccurence->nodeID);
	double variance = 1;

	if (target == getTwinNode(node) || target == node)
		return;

	if (position < 0) {
		variance += getNodeLength(node) * getNodeLength(node) / 16;
		// distance += 0;
	} else {
		// variance += 0;
		distance += position - getNodeLength(node) / 2;
	}

	if (readOccurence->position < 0) {
		variance +=
		    getNodeLength(target) * getNodeLength(target) / 16;
		//distance += 0;
	} else {
		// variance += 0;
		distance +=
		    -readOccurence->position + getNodeLength(target) / 2;
	}

	if (readOccurence->offset < 0 || offset < 0) { 
		variance += length * length / 16;
		//distance += 0;
	} else {
		// variance += 0;
		distance += readOccurence->offset - offset;
	}

	// Relative ordering
	if (offset > 0 && readOccurence->offset > 0) {
		if (offset < readOccurence->offset) {
			if (distance - getNodeLength(node)/2 - getNodeLength(target)/2 < -10)
				;
			else if (distance < getNodeLength(node)/2 + getNodeLength(target)/2)
				createConnection(getNodeID(node), getNodeID(target), 1, 0,
						 getNodeLength(node)/2 + getNodeLength(target)/2, variance);
			else
				createConnection(getNodeID(node), getNodeID(target), 1, 0,
						 distance, variance);
		} else if (offset > readOccurence->offset) {
			if (-distance - getNodeLength(node)/2 - getNodeLength(target)/2 < -10)
				;
			else if (-distance < getNodeLength(node)/2 + getNodeLength(target)/2)
				createConnection(-getNodeID(node), -getNodeID(target), 1,
						 0, getNodeLength(node)/2 + getNodeLength(target)/2 , variance);
			else 
				createConnection(-getNodeID(node), -getNodeID(target), 1,
						 0, -distance, variance);
		}
	} else if (offset > 0 && position > 0) {
		if (distance - offset > -getNodeLength(node)/2 && distance - offset + length > getNodeLength(node)/2)
			createConnection(getNodeID(node), getNodeID(target), 1, 0,
					 getNodeLength(node)/2 + getNodeLength(target)/2, variance);
		else if (distance - offset < -getNodeLength(node)/2 && distance - offset + length < getNodeLength(node)/2)
			createConnection(-getNodeID(node), -getNodeID(target), 1, 0,
					 getNodeLength(node)/2 + getNodeLength(target)/2, variance);
		else {
			createConnection(getNodeID(node), getNodeID(target), 1, 0,
					 getNodeLength(node)/2 + getNodeLength(target)/2, variance);
			createConnection(-getNodeID(node), -getNodeID(target), 1, 0,
					 getNodeLength(node)/2 + getNodeLength(target)/2, variance);
		}
	} else if (readOccurence->offset > 0 && readOccurence->position > 0) {
		if (-distance - readOccurence->offset > -getNodeLength(target)/2 && -distance - readOccurence->offset + length > getNodeLength(target)/2)
			createConnection(-getNodeID(node), -getNodeID(target), 1, 0,
					 getNodeLength(node)/2 + getNodeLength(target)/2, variance);
		if (-distance - readOccurence->offset < -getNodeLength(target)/2 && -distance - readOccurence->offset + length < getNodeLength(target)/2)
			createConnection(getNodeID(node), getNodeID(target), 1, 0,
					 getNodeLength(node)/2 + getNodeLength(target)/2, variance);
		else {
			createConnection(getNodeID(node), getNodeID(target), 1, 0,
					 getNodeLength(node)/2 + getNodeLength(target)/2, variance);
			createConnection(-getNodeID(node), -getNodeID(target), 1, 0,
					 getNodeLength(node)/2 + getNodeLength(target)/2, variance);
		}
	} else {
		createConnection(getNodeID(node), getNodeID(target), 1, 0,
				 getNodeLength(node)/2 + getNodeLength(target)/2, variance);
		createConnection(-getNodeID(node), -getNodeID(target), 1, 0,
				 getNodeLength(node)/2 + getNodeLength(target)/2, variance);
	}
}

static void projectFromReadPair(Node * node, ReadOccurence * readOccurence,
				Coordinate position, Coordinate offset,
				Coordinate insertLength,
				double insertVariance,
				boolean doMatePairs)
{
	Coordinate distance = insertLength;
	Coordinate variance = insertVariance;
	Node *target = getNodeInGraph(graph, readOccurence->nodeID);
	IDnum nodeID;
	IDnum node2ID;

	if (target == getTwinNode(node) || target == node)
		return;

	nodeID  = getNodeID(node);
	node2ID = getNodeID(target);

	if (getUniqueness(target) && node2ID < nodeID)
		return;

	// Check if a conflicting PE (or MP from a smaller size lib) connection
	// already exists
	if (doMatePairs) {
		Connection *reverseConnect;

#ifdef _OPENMP
		lockTwoNodes(nodeID, node2ID);
#endif
		reverseConnect = findConnection(-nodeID, -node2ID);
#ifdef _OPENMP
		unLockTwoNodes(nodeID, node2ID);
#endif

		if (reverseConnect != NULL &&
		    reverseConnect->clean &&
		    reverseConnect->paired_count +
		    reverseConnect->direct_count >= UNRELIABLE_CONNECTION_CUTOFF)
			return;
	}

	if (position < 0) {
		variance += getNodeLength(node) * getNodeLength(node) / 16;
		// distance += 0;
	} else {
		// variance += 0;
		distance += position - offset - getNodeLength(node) / 2;
	}

	if (readOccurence->position < 0) {
		variance +=
		    getNodeLength(target) * getNodeLength(target) / 16;
		//distance += 0;
	} else {
		// variance += 0;
		distance +=
		    readOccurence->position - readOccurence->offset -
		    getNodeLength(target) / 2;
	}

	if (distance - getNodeLength(node)/2 - getNodeLength(target)/2 < -6 * sqrt(insertVariance))
		return;
	else if (distance < getNodeLength(node)/2 + getNodeLength(target)/2)
		distance = getNodeLength(node)/2 + getNodeLength(target)/2;

	createConnection(nodeID, node2ID, 0, 1,
			 distance, variance);
}

static void projectFromShortRead(Node * node,
				 ShortReadMarker * shortMarker,
				 IDnum * readPairs, Category * cats,
				 ReadOccurence ** readNodes,
				 IDnum * readNodeCounts,
				 ShortLength * lengths,
				 boolean * shadows,
				 boolean doMatePairs,
				 Category thisCat)
{
	IDnum index;
	IDnum readIndex = getShortReadMarkerID(shortMarker);
	ReadOccurence *readArray;
	IDnum readPairIndex;
	Category cat;
	Coordinate position = getShortReadMarkerPosition(shortMarker);
	Coordinate offset = getShortReadMarkerOffset(shortMarker);
	Coordinate length = lengths[getShortReadMarkerID(shortMarker) - 1];
	Coordinate insertLength;
	double insertVariance;

	// Going through single-read information
	if (!doMatePairs && readNodeCounts[readIndex] > 1) {
		readArray = readNodes[readIndex];
		for (index = 0; index < readNodeCounts[readIndex]; index++)
			projectFromSingleRead(node, &readArray[index],
					      position, offset, length);
	}
	// Going through paired read information
	if (readPairs == NULL)
		return;

	readPairIndex = readPairs[readIndex - 1] + 1;

	if (readPairIndex == 0)
		return;

	cat = cats[readIndex - 1];
	insertLength = getInsertLength(graph, cat);
	insertVariance = getInsertLength_var(graph, cat);
	cat /= 2;
	if (shadows[cat] && cat > PEBBLE_ROUND_NUM)
		return;

	if (!shadows[cat] && !doMatePairs) {
		readArray = readNodes[readPairIndex];
		for (index = 0; index < readNodeCounts[readPairIndex]; index++)
			projectFromReadPair(node, &readArray[index], position,
					offset, insertLength, insertVariance, false);
	}
	else if (shadows[cat] && doMatePairs && cat == thisCat) {
		readArray = readNodes[readPairIndex];
		for (index = 0; index < readNodeCounts[readPairIndex]; index++)
			projectFromReadPair(node, &readArray[index], position,
					offset, insertLength, insertVariance, true);
	}

}

static void projectFromLongRead(Node * node, PassageMarkerI marker,
				IDnum * readPairs, Category * cats,
				ReadOccurence ** readNodes,
				IDnum * readNodeCounts,
				ShortLength * lengths)
{
	IDnum index;
	IDnum readIndex = getPassageMarkerSequenceID(marker);
	ReadOccurence *readArray;
	IDnum readPairIndex;
	Category cat;
	Coordinate position = getStartOffset(marker);
	Coordinate offset = getPassageMarkerStart(marker);
	Coordinate length =
	    lengths[getPassageMarkerSequenceID(marker) - 1];
	Coordinate insertLength;
	double insertVariance;

	// Going through single-read information
	if (readNodeCounts[readIndex] > 1 && position > 0) {
		readArray = readNodes[readIndex];
		for (index = 0; index < readNodeCounts[readIndex]; index++)
			projectFromSingleRead(node, &readArray[index],
					      position, offset, length);
	}
	// Going through paired read information
	if (readPairs == NULL)
		return;

	readPairIndex = readPairs[readIndex - 1] + 1;

	if (readPairIndex == 0)
		return;

	cat = cats[readIndex - 1];
	insertLength = getInsertLength(graph, cat);
	insertVariance = getInsertLength_var(graph, cat);

	readArray = readNodes[readPairIndex];
	for (index = 0; index < readNodeCounts[readPairIndex]; index++)
		projectFromReadPair(node, &readArray[index], position,
				    offset, insertLength, insertVariance, false);

}

static void projectFromNode(IDnum nodeID,
			    ReadOccurence ** readNodes,
			    IDnum * readNodeCounts,
			    IDnum * readPairs, Category * cats,
			    boolean * dubious, ShortLength * lengths,
			    boolean * shadows,
			    boolean doMatePairs,
			    Category thisCat)
{
	IDnum index;
	ShortReadMarker *nodeArray, *shortMarker;
	PassageMarkerI marker;
	Node *node;
	IDnum nodeReadCount;

	node = getNodeInGraph(graph, nodeID);

	if (node == NULL || !getUniqueness(node))
		return;

	nodeArray = getNodeReads(node, graph);
	nodeReadCount = getNodeReadCount(node, graph);
	for (index = 0; index < nodeReadCount; index++) {
		shortMarker = getShortReadMarkerAtIndex(nodeArray, index);
		if (dubious[getShortReadMarkerID(shortMarker) - 1])
			continue;
		projectFromShortRead(node, shortMarker, readPairs, cats,
				     readNodes, readNodeCounts, lengths,
				     shadows,
				     doMatePairs,
				     thisCat);
	}

	if (!doMatePairs)
		for (marker = getMarker(node); marker != NULL_IDX;
				marker = getNextInNode(marker)) {
			if (getPassageMarkerSequenceID(marker) > 0)
				projectFromLongRead(node, marker, readPairs, cats,
						readNodes, readNodeCounts,
						lengths);
		}
}

static Connection **computeNodeToNodeMappings(ReadOccurence ** readNodes,
					      IDnum * readNodeCounts,
					      IDnum * readPairs,
					      Category * cats,
					      boolean * dubious,
					      boolean * shadows,
					      ShortLength * lengths)
{
	IDnum nodeID;
	IDnum nodes = nodeCount(graph);
	struct timeval start, end, diff;
	Category cat;
	boolean hasShadow;

	scaffold = callocOrExit(2 * nodes + 1, Connection *);

	velvetLog("Computing direct node to node mappings\n");

	gettimeofday(&start, NULL);
#ifdef _OPENMP
	createNodeLocks(graph);

	int threads = omp_get_max_threads();
	if (threads > 32)
		threads = 32;

	#pragma omp parallel for num_threads(threads)
#endif 
	for (nodeID = -nodes; nodeID <= nodes; nodeID++)
	{
		if (nodeID % 10000 == 0)
			velvetLog("Scaffolding node %li\n", (long) nodeID);

		projectFromNode(nodeID, readNodes, readNodeCounts,
				readPairs, cats, dubious, lengths, shadows, false, 0);
	}

#ifdef _OPENMP
	initConnectionStackMemory();
#endif

	hasShadow = false;
	for (cat = 0; cat < CATEGORIES; cat++)
		if (shadows[cat])
		{
			hasShadow = true;
			break;
		}

	if (hasShadow)
	{
		for (cat = 0; cat < CATEGORIES; cat++)
		{
			setAllConnectionsClean();
			if (!shadows[cat])
				continue;
			velvetLog("Scaffolding MP library %i\n", cat);
#ifdef _OPENMP
			#pragma omp parallel for
#endif 
			for (nodeID = -nodes; nodeID <= nodes; nodeID++)
				projectFromNode(nodeID, readNodes, readNodeCounts,
						readPairs, cats, dubious, lengths,
						shadows, true, cat);
		}
	}
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (nodeID = 2 * nodes; nodeID >= 0; nodeID--)
		splayToList(scaffold + nodeID);

	destroyConnectionStackMemory();

#ifdef _OPENMP
	free(nodeLocks);
	nodeLocks = NULL;
#endif
	gettimeofday(&end, NULL);
	timersub(&end, &start, &diff);
	velvetLog(" === Nodes Scaffolded in %ld.%06ld s\n", (long) diff.tv_sec, (long) diff.tv_usec);

	PEBBLE_ROUND_NUM++;

	return scaffold;
}

static IDnum **countShortReads(Graph * graph, ReadSet * reads)
{
	IDnum **counts = callocOrExit(CATEGORIES + 1, IDnum *);
	Category cat;
	IDnum nodeIndex;
	IDnum nodes = nodeCount(graph);
	Node *node;
	ShortReadMarker *array, *marker;
	IDnum readCount, readIndex, readID;

	// Allocate memory where needed
	for (cat = 0; cat <= CATEGORIES; cat++)
		if (getInsertLength(graph, cat) > 0)
			counts[cat] =
			    callocOrExit(2 * nodeCount(graph) + 1,
				   IDnum);

	// Start fillin'
	for (nodeIndex = 0; nodeIndex < 2 * nodes + 1; nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex - nodes);

		if (node == NULL || !getUniqueness(node))
			continue;

		array = getNodeReads(node, graph);
		readCount = getNodeReadCount(node, graph);
		for (readIndex = 0; readIndex < readCount; readIndex++) {
			marker =
			    getShortReadMarkerAtIndex(array, readIndex);
			readID = getShortReadMarkerID(marker);
			cat = reads->categories[readID - 1];
			if (cat % 2 == 1 && counts[cat / 2] != NULL)
				counts[cat / 2][nodeIndex]++;
		}
	}

	return counts;
}

static void removeUnreliableConnections(ReadSet * reads, boolean *shadows)
{
	IDnum maxNodeIndex = nodeCount(graph) * 2 + 1;
	IDnum index;
	Connection *connect, *next;
	Category cat;
	IDnum **counts = countShortReads(graph, reads);
	IDnum nodes = nodeCount(graph);

	for (index = 0; index < maxNodeIndex; index++) {
		for (connect = scaffold[index]; connect != NULL;
		     connect = next) {
			next = connect->right;
			if (!testConnection(index - nodes, connect, counts, shadows))
				destroyConnection(connect, index - nodes);
		}
	}

	// Free memory
	for (cat = 0; cat <= CATEGORIES; cat++)
		if (counts[cat])
			free(counts[cat]);
	free(counts);
}

void printConnections(ReadSet * reads, boolean * shadows)
{
	IDnum maxNodeIndex = nodeCount(graph) * 2 + 1;
	IDnum index;
	Connection *connect, *next;
	Node *node;
	IDnum **counts = countShortReads(graph, reads);
	IDnum nodes = nodeCount(graph);
	Category cat;

	puts("CONNECT IDA IDB dcount pcount dist lengthA lengthB var countA countB coordA coordB real exp distance test");

	for (index = 0; index < maxNodeIndex; index++) {
		node = getNodeInGraph(graph, index - nodeCount(graph));
		for (connect = scaffold[index]; connect != NULL;
		     connect = next) {
			next = getNextConnection(connect);
			printf
			    ("CONNECT %ld %ld %ld %ld %lld %lld %lld %f %ld %ld",
			     (long) index - nodeCount(graph),
			     (long) getNodeID(connect->destination),
			     (long) connect->direct_count,
			     (long) connect->paired_count,
			     (long long) getConnectionDistance(connect),
			     (long long) getNodeLength(node), (long long)
			     getNodeLength(connect->destination),
			     connect->variance,
			     (long) getNodeReadCount(node, graph),
			     (long) getNodeReadCount(connect->destination,
						     graph));
			if (markerCount(node) == 1
			    && markerCount(connect->destination) == 1)
				printf(" %lld %lld %lld", (long long)
				       getPassageMarkerFinish(getMarker
							      (node)),
				       (long long)
				       getPassageMarkerFinish(getMarker
							      (connect->
							       destination)),
				       (long
					long) (getPassageMarkerFinish
					       (getMarker(node)) -
					       getPassageMarkerFinish
					       (getMarker
						(connect->destination))));
			else
				printf(" ? ? ?");
			printf(" %ld",
			       (long) expectedNumberOfConnections(index -
								  nodeCount
								  (graph),
								  connect,
								  counts,
								  0));
			printf(" %lld",
			       (long long) (getConnectionDistance(connect)
					    - (getNodeLength(node) +
					       getNodeLength
					       (connect->destination)) /
					    2));
			if (testConnection(index - nodes, connect, counts, shadows))
				puts(" OK");
			else
				puts(" NG");
		}
	}

	for (cat = 0; cat <= CATEGORIES; cat++)
		if (counts[cat])
			free(counts[cat]);
	free(counts);
}

void buildScaffold(Graph * argGraph,
		   ReadSet * reads,
		   boolean * dubious,
		   boolean * shadows)
{
	IDnum *readPairs;
	Category *cats;
	IDnum *readNodeCounts;
	ReadOccurence **readNodes;
	ReadOccurence *readNodesArray = NULL;
	ShortLength *lengths = getSequenceLengths(reads, getWordLength(argGraph));
	Coordinate totalCount = 0;

	graph = argGraph;
	readPairs = reads->mateReads;
	cats = reads->categories;

	// Prepare primary scaffold
	readNodeCounts = computeReadToNodeCounts(&totalCount);
	readNodes = computeReadToNodeMappings(readNodeCounts, reads, totalCount, &readNodesArray);

	estimateMissingInsertLengths(readNodes, readNodeCounts, readPairs, cats);

	scaffold = computeNodeToNodeMappings(readNodes, readNodeCounts,
				      readPairs, cats, dubious, shadows, lengths);
	removeUnreliableConnections(reads, shadows);

	free(readNodesArray);
	free(readNodes);
	free(readNodeCounts);
	free(lengths);
}

//DEBUG
void printScaffold(Graph * argGraph,
		   ReadSet * reads,
		   boolean * dubious,
		   boolean * shadows)
{
	IDnum *readPairs;
	Category *cats;
	IDnum *readNodeCounts;
	ReadOccurence **readNodes;
	ReadOccurence *readNodesArray = NULL;
	ShortLength *lengths = getSequenceLengths(reads, getWordLength(argGraph));
	Coordinate totalCount = 0;

	graph = argGraph;
	readPairs = reads->mateReads;
	cats = reads->categories;

	// Prepare primary scaffold
	readNodeCounts = computeReadToNodeCounts(&totalCount);
	readNodes = computeReadToNodeMappings(readNodeCounts, reads, totalCount, &readNodesArray);

	estimateMissingInsertLengths(readNodes, readNodeCounts, readPairs, cats);

	scaffold = computeNodeToNodeMappings(readNodes, readNodeCounts,
				      readPairs, cats, dubious, shadows, lengths);
	printConnections(reads, shadows);

	free(readNodesArray);
	free(readNodes);
	free(readNodeCounts);
	free(lengths);
	cleanScaffoldMemory();
}

void setUnreliableConnectionCutoff(int val)
{
	UNRELIABLE_CONNECTION_CUTOFF = (IDnum) val;
}

void cleanScaffoldMemory() {
	Category libID;

	for (libID = 0; libID < CATEGORIES + 1; libID++)
		if (estimated[libID])
			setInsertLengths(graph, libID, -1, -1);

	destroyRecycleBin(connectionMemory);
	free(scaffold);
	connectionMemory = NULL;
}

void setPairedExpFraction(double x) {
	paired_exp_fraction = x;
}
