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

typedef struct readOccurence_st ReadOccurence;

struct connection_st {
	Node *destination;
	Connection *next;
	Connection *previous;
	Connection *twin;
	double distance;
	double variance;
	IDnum direct_count;
	IDnum paired_count;
};

struct readOccurence_st {
	Coordinate position;
	Coordinate offset;
	IDnum nodeID;
};

// Global params
static IDnum UNRELIABLE_CONNECTION_CUTOFF = 5;

// Global pointers
static Graph *graph;
static Connection **scaffold = NULL;
static RecycleBin *connectionMemory = NULL;
static boolean estimated[CATEGORIES + 1];

static Connection *allocateConnection()
{
	if (connectionMemory == NULL)
		connectionMemory =
		    newRecycleBin(sizeof(Connection), BLOCK_SIZE);

	return allocatePointer(connectionMemory);
}

static void deallocateConnection(Connection * connect)
{
	deallocatePointer(connectionMemory, connect);
}

Node * getConnectionDestination(Connection * connect) {
	return connect->destination;
}

Connection * getNextConnection(Connection * connect) {
	return connect->next;
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

	//printf("Destroying connection from %li to %li\n", nodeID, getNodeID(connect->destination));

	if (connect == NULL)
		return;

	previous = connect->previous;
	next = connect->next;

	if (previous != NULL)
		previous->next = next;
	if (next != NULL)
		next->previous = previous;

	if (scaffold[nodeID + nodeCount(graph)] == connect)
		scaffold[nodeID + nodeCount(graph)] = next;

	if (connect->twin != NULL) {
		connect->twin->twin = NULL;
		destroyConnection(connect->twin,
				  getNodeID(connect->destination));
	}

	deallocateConnection(connect);
}

static boolean testConnection(IDnum IDA, Connection * connect,
			      IDnum ** counts)
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

	for (cat = 0; cat <= CATEGORIES; cat++)
		total +=
		    expectedNumberOfConnections(IDA, connect, counts, cat);

	// Remove inconsistent connections
	return connect->paired_count >= total / 10;
}

static IDnum *computeReadToNodeCounts()
{
	IDnum readIndex, nodeIndex;
	IDnum maxNodeIndex = 2 * nodeCount(graph) + 1;
	IDnum maxReadIndex = sequenceCount(graph) + 1;
	IDnum *readNodeCounts = callocOrExit(maxReadIndex, IDnum);
	boolean *readMarker = callocOrExit(maxReadIndex, boolean);
	ShortReadMarker *nodeArray, *shortMarker;
	PassageMarker *marker;
	Node *node;
	IDnum nodeReadCount;

	// Original
	/*
	puts("Computing read to node mapping array sizes");
	*/
	// Original

	for (nodeIndex = 0; nodeIndex < maxNodeIndex; nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex - nodeCount(graph));
		if (node == NULL)
			continue;
		nodeArray = getNodeReads(node, graph);
		nodeReadCount = getNodeReadCount(node, graph);

		// Short reads
		for (readIndex = 0; readIndex < nodeReadCount; readIndex++) {
			shortMarker =
			    getShortReadMarkerAtIndex(nodeArray,
						      readIndex);
			readNodeCounts[getShortReadMarkerID
				       (shortMarker)]++;
		}

		// Long reads
		for (marker = getMarker(node); marker != NULL;
		     marker = getNextInNode(marker)) {
			readIndex = getPassageMarkerSequenceID(marker);
			if (readIndex < 0)
				continue;

			if (readMarker[readIndex])
				continue;

			readNodeCounts[readIndex]++;
			readMarker[readIndex] = true;
		}

		// Clean up marker array
		for (marker = getMarker(node); marker != NULL;
		     marker = getNextInNode(marker)) {
			readIndex = getPassageMarkerSequenceID(marker);
			if (readIndex > 0)
				readMarker[readIndex] = false;
		}
	}

	free(readMarker);
	return readNodeCounts;
}

static ReadOccurence **allocateReadToNodeTables(IDnum * readNodeCounts)
{
	IDnum readIndex;
	IDnum maxReadIndex = sequenceCount(graph) + 1;
	ReadOccurence **readNodes =
	    callocOrExit(maxReadIndex, ReadOccurence *);

	for (readIndex = 1; readIndex < maxReadIndex; readIndex++) {
		if (readNodeCounts[readIndex] != 0) {
			readNodes[readIndex] =
			    callocOrExit(readNodeCounts[readIndex],
				   ReadOccurence);
			readNodeCounts[readIndex] = 0;
		}
	}

	return readNodes;
}

static void computePartialReadToNodeMapping(IDnum nodeID,
					    ReadOccurence ** readNodes,
					    IDnum * readNodeCounts,
					    boolean * readMarker)
{
	ShortReadMarker *shortMarker;
	IDnum index, readIndex;
	ReadOccurence *readArray, *readOccurence;
	Node *node = getNodeInGraph(graph, nodeID);
	ShortReadMarker *nodeArray = getNodeReads(node, graph);
	IDnum nodeReadCount = getNodeReadCount(node, graph);
	PassageMarker *marker;

	for (index = 0; index < nodeReadCount; index++) {
		shortMarker = getShortReadMarkerAtIndex(nodeArray, index);
		readIndex = getShortReadMarkerID(shortMarker);
		readArray = readNodes[readIndex];
		readOccurence = &readArray[readNodeCounts[readIndex]];
		readOccurence->nodeID = nodeID;
		readOccurence->position =
		    getShortReadMarkerPosition(shortMarker);
		readOccurence->offset =
		    getShortReadMarkerOffset(shortMarker);
		readNodeCounts[readIndex]++;
	}

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker)) {
		readIndex = getPassageMarkerSequenceID(marker);
		if (readIndex < 0)
			continue;

		if (!readMarker[readIndex]) {
			readArray = readNodes[readIndex];
			readOccurence =
			    &readArray[readNodeCounts[readIndex]];
			readOccurence->nodeID = nodeID;
			readOccurence->position = getStartOffset(marker);
			readOccurence->offset =
			    getPassageMarkerStart(marker);
			readNodeCounts[readIndex]++;
			readMarker[readIndex] = true;
		} else {
			readArray = readNodes[readIndex];
			readOccurence =
			    &readArray[readNodeCounts[readIndex] - 1];
			readOccurence->position = -1;
			readOccurence->offset = -1;
		}
	}

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker)) {
		readIndex = getPassageMarkerSequenceID(marker);
		if (readIndex > 0)
			readMarker[readIndex] = false;
	}
}

static ReadOccurence **computeReadToNodeMappings(IDnum * readNodeCounts)
{
	IDnum nodeID;
	IDnum nodes = nodeCount(graph);
	ReadOccurence **readNodes =
	    allocateReadToNodeTables(readNodeCounts);
	boolean *readMarker =
	    callocOrExit(sequenceCount(graph) + 1, boolean);

	// Original
	/*
	puts("Computing read to node mappings");
	*/
	// Original

	for (nodeID = -nodes; nodeID <= nodes; nodeID++)
		if (nodeID != 0 && getNodeInGraph(graph, nodeID))
			computePartialReadToNodeMapping(nodeID, readNodes,
							readNodeCounts,
							readMarker);

	free(readMarker);
	return readNodes;
}

static boolean * countCoOccurences(IDnum * coOccurencesCount, ReadOccurence ** readNodes, IDnum * readNodeCounts, IDnum * readPairs, Category * cats) {
	IDnum readIndex, readPairIndex;
	IDnum readNodeCount;
	IDnum readOccurenceIndex, readPairOccurenceIndex;
	ReadOccurence * readOccurence, *readPairOccurence;
	boolean * interestingReads = callocOrExit(sequenceCount(graph), boolean);
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
		libID = cats[readIndex]/2;
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
					interestingReads[readIndex] = true;
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

static void measureCoOccurences(Coordinate ** coOccurences, boolean * interestingReads, ReadOccurence ** readNodes, IDnum * readNodeCounts, IDnum * readPairs, Category * cats) {
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
		if (!interestingReads[readIndex])
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
	Coordinate * cA = (Coordinate *) A;
	Coordinate * cB = (Coordinate *) B;

	if (*cA > *cB)
		return 1;
	if (*cA == *cB)
		return 0;
	return -1;	
}

static void estimateLibraryInsertLength(Coordinate * coOccurences, IDnum coOccurencesCount, Category libID) {
	Coordinate median, variance;
	IDnum index;
	int counter = 0;
	qsort(coOccurences, coOccurencesCount, sizeof(Coordinate), compareReadOccurences);

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

	printf("Paired-end library %i has length: %lli, sample standard deviation: %lli\n", libID + 1, (long long) median, (long long) sqrt(variance));
	setInsertLengths(graph, libID, median, sqrt(variance));
	estimated[libID] = true;
}

static void estimateLibraryInsertLengths(Coordinate ** coOccurences, IDnum * coOccurencesCounts) {
	Category libID;

	for (libID = 0; libID < CATEGORIES + 1; libID++)
		estimated[libID] = false;

	for (libID = 0; libID < CATEGORIES + 1; libID++)
		if (coOccurencesCounts[libID] > 0)
			estimateLibraryInsertLength(coOccurences[libID], coOccurencesCounts[libID], libID);
}

static void estimateMissingInsertLengths(ReadOccurence ** readNodes, IDnum * readNodeCounts, IDnum * readPairs, Category * cats) {
	Coordinate * coOccurences[CATEGORIES + 1];
	IDnum coOccurencesCounts[CATEGORIES + 1]; 
	Category libID;

	// Original
	/*
	puts("Estimating library insert lengths...");
	*/
	// Original

	boolean * interestingReads = countCoOccurences(coOccurencesCounts, readNodes, readNodeCounts, readPairs, cats);

	for (libID = 0; libID < CATEGORIES + 1; libID++)
		coOccurences[libID] = callocOrExit(coOccurencesCounts[libID], Coordinate);

	measureCoOccurences(coOccurences, interestingReads, readNodes, readNodeCounts, readPairs, cats);
	estimateLibraryInsertLengths(coOccurences, coOccurencesCounts);

	for (libID = 0; libID < CATEGORIES + 1; libID++)
		free(coOccurences[libID]);
	
	free(interestingReads);

	// Original
	/*
	puts("Done");
	*/
	// Original
}

static Connection *findConnection(IDnum nodeID, IDnum node2ID)
{
	Node *node2 = getNodeInGraph(graph, node2ID);
	Connection *connect;

	if (node2 == NULL)
		return NULL;

	for (connect = scaffold[nodeID + nodeCount(graph)];
	     connect != NULL; connect = connect->next)
		if (connect->destination == node2)
			break;

	return connect;
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
	newConnection->previous = NULL;
	newConnection->next = scaffold[nodeIndex];
	if (scaffold[nodeIndex] != NULL)
		scaffold[nodeIndex]->previous = newConnection;
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
	connect->previous = NULL;
	connect->next = scaffold[nodeIndex];
	if (scaffold[nodeIndex] != NULL)
		scaffold[nodeIndex]->previous = connect;
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

static void createConnection(IDnum nodeID, IDnum node2ID,
			     IDnum direct_count,
			     IDnum paired_count,
			     Coordinate distance, double variance)
{
	Connection *connect = findConnection(nodeID, node2ID);

	if (connect != NULL)
		readjustConnection(connect, distance, variance,
				   direct_count, paired_count);
	else
		createNewConnection(nodeID, node2ID, direct_count,
				    paired_count, distance, variance);
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
				double insertVariance)
{
	Coordinate distance = insertLength;
	Coordinate variance = insertVariance;
	Node *target = getNodeInGraph(graph, readOccurence->nodeID);

	if (target == getTwinNode(node) || target == node)
		return;

	if (getUniqueness(target) && getNodeID(target) < getNodeID(node))
		return;

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

	createConnection(getNodeID(node), getNodeID(target), 0, 1,
			 distance, variance);
}

static void projectFromShortRead(Node * node,
				 ShortReadMarker * shortMarker,
				 IDnum * readPairs, Category * cats,
				 ReadOccurence ** readNodes,
				 IDnum * readNodeCounts,
				 Coordinate * lengths)
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
	if (readNodeCounts[readIndex] > 1) {
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
				    offset, insertLength, insertVariance);

}

static void projectFromLongRead(Node * node, PassageMarker * marker,
				IDnum * readPairs, Category * cats,
				ReadOccurence ** readNodes,
				IDnum * readNodeCounts,
				Coordinate * lengths)
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
				    offset, insertLength, insertVariance);

}

static void projectFromNode(IDnum nodeID,
			    ReadOccurence ** readNodes,
			    IDnum * readNodeCounts,
			    IDnum * readPairs, Category * cats,
			    boolean * dubious, Coordinate * lengths)
{
	IDnum index;
	ShortReadMarker *nodeArray, *shortMarker;
	PassageMarker *marker;
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
				     readNodes, readNodeCounts, lengths);
	}

	for (marker = getMarker(node); marker != NULL;
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
					      Coordinate * lengths)
{
	IDnum nodeID;
	IDnum nodes = nodeCount(graph);
	scaffold = callocOrExit(2 * nodes + 1, Connection *);

	// Original
	/*
	puts("Computing direct node to node mappings");
	*/
	// Original

	for (nodeID = -nodes; nodeID <= nodes; nodeID++) {
		// Original
		/*
		if (nodeID % 10000 == 0)
			printf("Scaffolding node %d\n", nodeID);
		*/
		// Original

		projectFromNode(nodeID, readNodes, readNodeCounts,
				readPairs, cats, dubious, lengths);
	}

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

void printConnections(ReadSet * reads)
{
	IDnum maxNodeIndex = nodeCount(graph) * 2 + 1;
	IDnum index;
	Connection *connect, *next;
	Node *node;
	IDnum **counts = countShortReads(graph, reads);
	IDnum nodes = nodeCount(graph);

	puts("CONNECT IDA IDB dcount pcount dist lengthA lengthB var countA countB coordA coordB real exp distance test");

	for (index = 0; index < maxNodeIndex; index++) {
		node = getNodeInGraph(graph, index - nodeCount(graph));
		for (connect = scaffold[index]; connect != NULL;
		     connect = next) {
			next = connect->next;
			if (getUniqueness(connect->destination)) {
				printf
				    ("CONNECT %ld %ld %ld %ld %lld %lld %lld %f %ld %ld",
				     (long) index - nodeCount(graph),
				     (long) getNodeID(connect->destination),
				     (long) connect->direct_count,
				     (long) connect->paired_count,
				     (long long) getConnectionDistance(connect),
				     (long long) getNodeLength(node),
				     (long long) getNodeLength(connect->destination),
				     connect->variance,
				     (long) getNodeReadCount(node, graph),
				     (long) getNodeReadCount(connect->destination,
						      graph));
				if (markerCount(node) == 1
				    && markerCount(connect->destination) ==
				    1)
					printf(" %lld %lld %lld",
					       (long long) getPassageMarkerFinish
					       (getMarker(node)),
					       (long long) getPassageMarkerFinish
					       (getMarker
						(connect->destination)),
					       (long long) (getPassageMarkerFinish
					       (getMarker(node)) - 
					       getPassageMarkerFinish
					       (getMarker
						(connect->destination))));
				else
					printf(" ? ?");
				printf(" %ld", (long) expectedNumberOfConnections(index-nodeCount(graph), connect, counts, 0));
				printf(" %lld", (long long) (getConnectionDistance(connect) - (getNodeLength(node) + getNodeLength(connect->destination))/2));
				if (testConnection
				    (index - nodes, connect, counts))
					puts(" OK");
				else
					puts(" NG");
			}
		}
	}
}

static void removeUnreliableConnections(ReadSet * reads)
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
			next = connect->next;
			if (!testConnection
			    (index - nodes, connect, counts))
				destroyConnection(connect, index - nodes);
		}
	}

	// Free memory
	for (cat = 0; cat <= CATEGORIES; cat++)
		if (counts[cat])
			free(counts[cat]);
	free(counts);
}

void buildScaffold(Graph * argGraph, ReadSet * reads, boolean * dubious) {
	IDnum *readPairs;
	Category *cats;
	IDnum *readNodeCounts;
	ReadOccurence **readNodes;
	Coordinate *lengths =
	    getSequenceLengths(reads, getWordLength(argGraph));
	IDnum index;

	graph = argGraph;
	readPairs = reads->mateReads;
	cats = reads->categories;

	// Prepare primary scaffold
	readNodeCounts = computeReadToNodeCounts();
	readNodes = computeReadToNodeMappings(readNodeCounts);

	estimateMissingInsertLengths(readNodes, readNodeCounts, readPairs, cats);

	scaffold = computeNodeToNodeMappings(readNodes, readNodeCounts,
				      readPairs, cats, dubious, lengths);
	removeUnreliableConnections(reads);

	// Clean up memory
	for (index = 1; index <= sequenceCount(graph); index++)
		free(readNodes[index]);

	free(readNodes);
	free(readNodeCounts);
	free(lengths);
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
