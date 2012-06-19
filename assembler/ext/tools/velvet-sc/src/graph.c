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
#include <string.h>

#include "globals.h"
#include "graph.h"
#include "recycleBin.h"
#include "tightString.h"
#include "passageMarker.h"
#include "utility.h"
#include "kmer.h"

#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

struct arc_st {
	Arc *twinArc;		// 64
	Arc *next;		// 64
	Arc *previous;		// 64
	Arc *nextInLookupTable;	// 64
	Node *destination;	// 64
	IDnum multiplicity;	// 32
};				// 352 Total

struct node_st {
	Node *twinNode;		// 64
	Arc *arc;		// 64
	Descriptor *descriptor;	// 64
	PassageMarker *marker;	// 64
	Coordinate length;	// 32
	Coordinate virtualCoverage[CATEGORIES];	// 32 * 2
	Coordinate originalVirtualCoverage[CATEGORIES];	// 32 * 2
	IDnum ID;		// 32
	IDnum arcCount;		// 32
	boolean status;		// 1
	boolean uniqueness;	// 1
};				// 418 Total

struct shortReadMarker_st {
	Coordinate position;
	IDnum readID;
	ShortLength offset;
};

struct gapMarker_st {
	GapMarker *next;
	Coordinate position;
	ShortLength length;
};

struct graph_st {
	IDnum sequenceCount;
	IDnum nodeCount;
	Node **nodes;
	Arc **arcLookupTable;
	ShortReadMarker **nodeReads;
	IDnum *nodeReadCounts;
	Coordinate insertLengths[CATEGORIES + 1];
	double insertLengths_var[CATEGORIES + 1];
	int wordLength;
	GapMarker **gapMarkers;
};

static RecycleBin *arcMemory = NULL;
static RecycleBin *nodeMemory = NULL;
static RecycleBin *gapMarkerMemory = NULL;

#define BLOCKSIZE 50
#define GAPBLOCKSIZE 10000

Arc *allocateArc()
{
	if (arcMemory == NULL)
		arcMemory = newRecycleBin(sizeof(Arc), BLOCKSIZE);

	return allocatePointer(arcMemory);
}

void deallocateArc(Arc * arc)
{
	deallocatePointer(arcMemory, arc);
}

Node *allocateNode()
{
	if (nodeMemory == NULL)
		nodeMemory = newRecycleBin(sizeof(Node), BLOCKSIZE);

	return (Node *) allocatePointer(nodeMemory);
}

void deallocateNode(Node * node)
{
	deallocatePointer(nodeMemory, node);
}

// Returns the twin node of a given node
Node *getTwinNode(Node * node)
{
	return node->twinNode;
}

// Inserts new passage marker in the marker list of destination node
void insertPassageMarker(PassageMarker * marker, Node * destination)
{
	setTopOfTheNode(marker);
	setNextInNode(marker, destination->marker);
	destination->marker = marker;
}

// Returns the length of the node's descriptor list
Coordinate getNodeLength(Node * node)
{
	return node->length;
}

// Returns the number of nodes in the graph
IDnum nodeCount(Graph * graph)
{
	return graph->nodeCount;
}

// returns the number of sequences used to buid the graph
IDnum sequenceCount(Graph * graph)
{
	return graph->sequenceCount;
}

// Creates an arc from node origin to node destination.
// If this arc already exists, increments its multiplicity by 1.
Arc *createArc(Node * originNode, Node * destinationNode, Graph * graph)
{
	Arc *arc, *twinArc;
	Node *destinationTwin;
	IDnum lookupIndex;

	if (originNode == NULL || destinationNode == NULL)
		return NULL;

//      printf("Connecting nodes %i -> %i\n", originNode->ID, destinationNode->ID);

	arc = getArcBetweenNodes(originNode, destinationNode, graph);

	if (arc != NULL) {
		arc->multiplicity++;
		arc->twinArc->multiplicity++;
		return arc;
	}
	// If not found
	arc = allocateArc();
	arc->destination = destinationNode;
	arc->multiplicity = 1;
	arc->previous = NULL;
	arc->next = originNode->arc;
	if (originNode->arc != NULL)
		originNode->arc->previous = arc;
	originNode->arc = arc;
	originNode->arcCount++;

	destinationTwin = destinationNode->twinNode;

	// Hairpin case
	if (destinationTwin == originNode) {
		arc->multiplicity++;
		arc->twinArc = arc;
		if (graph->arcLookupTable != NULL) {
			lookupIndex =
			    2 * originNode->ID + destinationNode->ID +
			    3 * graph->nodeCount;
			arc->nextInLookupTable =
			    graph->arcLookupTable[lookupIndex];
			graph->arcLookupTable[lookupIndex] = arc;
		}
		return arc;
	}

	twinArc = allocateArc();
	twinArc->destination = originNode->twinNode;
	twinArc->multiplicity = 1;
	twinArc->previous = NULL;
	twinArc->next = destinationTwin->arc;
	if (destinationTwin->arc != NULL)
		destinationTwin->arc->previous = twinArc;
	destinationTwin->arc = twinArc;
	destinationTwin->arcCount++;

	arc->twinArc = twinArc;
	twinArc->twinArc = arc;

	if (graph->arcLookupTable != NULL) {
		lookupIndex =
		    2 * originNode->ID + destinationNode->ID +
		    3 * graph->nodeCount;
		arc->nextInLookupTable =
		    graph->arcLookupTable[lookupIndex];
		graph->arcLookupTable[lookupIndex] = arc;

		lookupIndex =
		    -2 * destinationNode->ID - originNode->ID +
		    3 * graph->nodeCount;
		twinArc->nextInLookupTable =
		    graph->arcLookupTable[lookupIndex];
		graph->arcLookupTable[lookupIndex] = twinArc;
	}
	return arc;
}

void createAnalogousArc(Node * originNode, Node * destinationNode,
			Arc * refArc, Graph * graph)
{
	Arc *arc, *twinArc;
	Node *destinationTwin;
	IDnum lookupIndex;

	if (originNode == NULL || destinationNode == NULL)
		return;

//      printf("Connecting nodes %i -> %i\n", originNode->ID, destinationNode->ID);

	arc = getArcBetweenNodes(originNode, destinationNode, graph);

	if (arc != NULL) {
		if (refArc->twinArc != refArc) {
			arc->multiplicity += getMultiplicity(refArc);
			arc->twinArc->multiplicity +=
			    getMultiplicity(refArc);
		} else {
			arc->multiplicity += getMultiplicity(refArc) / 2;
			arc->twinArc->multiplicity +=
			    getMultiplicity(refArc) / 2;
		}
		return;
	}
	// If not found
	arc = allocateArc();
	arc->destination = destinationNode;
	arc->multiplicity = getMultiplicity(refArc);
	arc->previous = NULL;
	arc->next = originNode->arc;
	if (originNode->arc != NULL)
		originNode->arc->previous = arc;
	originNode->arc = arc;
	originNode->arcCount++;

	destinationTwin = destinationNode->twinNode;

	// Hairpin case
	if (destinationTwin == originNode) {
		arc->twinArc = arc;
		if (refArc->twinArc != refArc)
			arc->multiplicity *= 2;

		if (graph->arcLookupTable != NULL) {
			lookupIndex =
			    2 * originNode->ID + destinationNode->ID
			    + 3 * graph->nodeCount;
			arc->nextInLookupTable =
			    graph->arcLookupTable[lookupIndex];
			graph->arcLookupTable[lookupIndex] = arc;
		}
		return;
	}

	twinArc = allocateArc();
	twinArc->destination = originNode->twinNode;
	twinArc->multiplicity = getMultiplicity(refArc);
	twinArc->previous = NULL;
	twinArc->next = destinationTwin->arc;
	if (destinationTwin->arc != NULL)
		destinationTwin->arc->previous = twinArc;
	destinationTwin->arc = twinArc;
	destinationTwin->arcCount++;

	arc->twinArc = twinArc;
	twinArc->twinArc = arc;

	if (graph->arcLookupTable != NULL) {
		lookupIndex =
		    2 * originNode->ID + destinationNode->ID +
		    3 * graph->nodeCount;
		arc->nextInLookupTable =
		    graph->arcLookupTable[lookupIndex];
		graph->arcLookupTable[lookupIndex] = arc;

		lookupIndex =
		    -2 * destinationNode->ID - originNode->ID +
		    3 * graph->nodeCount;
		twinArc->nextInLookupTable =
		    graph->arcLookupTable[lookupIndex];
		graph->arcLookupTable[lookupIndex] = twinArc;
	}
}

void changeMultiplicity(Arc * arc, IDnum variation)
{
	if (arc == NULL)
		return;
	arc->multiplicity += variation;
	arc->twinArc->multiplicity += variation;
}

Arc *getArcBetweenNodes(Node * originNode, Node * destinationNode,
			Graph * graph)
{
	Arc *arc;
	Node *twinDestination, *twinOrigin;

	if (originNode == NULL || destinationNode == NULL)
		return NULL;

	if (graph->arcLookupTable != NULL) {
		for (arc =
		     graph->arcLookupTable[2 * originNode->ID +
					   destinationNode->ID +
					   3 * graph->nodeCount];
		     arc != NULL; arc = arc->nextInLookupTable) {
			if (arc->destination == destinationNode) {
				return arc;
			}
		}
		return NULL;
	}

	twinDestination = destinationNode->twinNode;
	if (originNode->arcCount <= twinDestination->arcCount) {
		for (arc = originNode->arc; arc != NULL; arc = arc->next)
			if (arc->destination == destinationNode)
				return arc;
		return NULL;
	}

	twinOrigin = originNode->twinNode;
	for (arc = twinDestination->arc; arc != NULL; arc = arc->next)
		if (arc->destination == twinOrigin)
			return arc->twinArc;
	return NULL;
}

void destroyArc(Arc * arc, Graph * graph)
{
	Node *origin, *destination;
	Arc *twinArc;
	Arc *currentArc;
	IDnum lookupIndex;

	if (arc == NULL)
		return;

	twinArc = arc->twinArc;
	origin = twinArc->destination->twinNode;
	destination = arc->destination->twinNode;

	//printf("Destroying arc %p\n", arc);

	// Removing arc from list
	if (origin->arc == arc) {
		origin->arc = arc->next;
		if (origin->arc != NULL)
			origin->arc->previous = NULL;
	} else {
		arc->previous->next = arc->next;
		if (arc->next != NULL)
			arc->next->previous = arc->previous;
	}

	origin->arcCount--;

	if (destination == origin) {
		if (graph->arcLookupTable != NULL) {
			lookupIndex =
			    2 * origin->ID - destination->ID +
			    3 * graph->nodeCount;
			currentArc = graph->arcLookupTable[lookupIndex];
			if (currentArc == arc)
				graph->arcLookupTable[lookupIndex] =
				    arc->nextInLookupTable;
			else {
				while (currentArc->nextInLookupTable !=
				       arc)
					currentArc =
					    currentArc->nextInLookupTable;

				currentArc->nextInLookupTable =
				    twinArc->nextInLookupTable;
			}
		}

		deallocateArc(arc);
		return;
	}
	// Removing arc's twin from list
	if (destination->arc == twinArc) {
		destination->arc = twinArc->next;
		if (destination->arc != NULL)
			destination->arc->previous = NULL;
	} else {
		twinArc->previous->next = twinArc->next;
		if (twinArc->next != NULL)
			twinArc->next->previous = twinArc->previous;
	}

	destination->arcCount--;

	if (graph->arcLookupTable != NULL) {
		lookupIndex =
		    2 * origin->ID - destination->ID +
		    3 * graph->nodeCount;
		currentArc = graph->arcLookupTable[lookupIndex];
		if (currentArc == arc)
			graph->arcLookupTable[lookupIndex] =
			    arc->nextInLookupTable;
		else {
			while (currentArc->nextInLookupTable != arc)
				currentArc = currentArc->nextInLookupTable;

			currentArc->nextInLookupTable =
			    arc->nextInLookupTable;
		}

		lookupIndex =
		    2 * destination->ID - origin->ID +
		    3 * graph->nodeCount;
		currentArc = graph->arcLookupTable[lookupIndex];
		if (currentArc == twinArc)
			graph->arcLookupTable[lookupIndex] =
			    twinArc->nextInLookupTable;
		else {
			while (currentArc->nextInLookupTable != twinArc)
				currentArc = currentArc->nextInLookupTable;

			currentArc->nextInLookupTable =
			    twinArc->nextInLookupTable;
		}
	}
	// Freeing memory
	deallocateArc(arc);
	deallocateArc(twinArc);
}

void destroyNode(Node * node, Graph * graph)
{
	Node *twin = node->twinNode;
	IDnum ID = node->ID;
	IDnum index;

	//printf("Destroying %d\n and twin %d\n", getNodeID(node), getNodeID(twin));

	if (ID < 0)
		ID = -ID;

	// Node arcs:
	while (node->arc != NULL)
		destroyArc(node->arc, graph);
	while (twin->arc != NULL)
		destroyArc(twin->arc, graph);

	// Descriptors
	free(node->descriptor);
	free(twin->descriptor);

	// Passage markers
	while (node->marker != NULL)
		destroyPassageMarker(node->marker);

	// Reads starts
	if (graph->nodeReads != NULL) {
		index = ID + graph->nodeCount;
		free(graph->nodeReads[index]);
		graph->nodeReads[index] = NULL;
		graph->nodeReadCounts[index] = 0;

		index = -ID + graph->nodeCount;
		free(graph->nodeReads[index]);
		graph->nodeReads[index] = NULL;
		graph->nodeReadCounts[index] = 0;
	}

	graph->nodes[ID] = NULL;
	deallocateNode(node);
	deallocateNode(twin);
}

int outDegree(Node * node)
{
	int result = 0;
	Arc *arc = node->arc;
	while (arc != NULL) {
		result += arc->multiplicity;
		arc = arc->next;
	}

	return result;
}

int simpleArcCount(Node * node)
{
	return node->arcCount;
}

int arcCount(Node * node)
{
	int result = 0;
	Arc *arc;

	if (node == NULL)
		return result;

	arc = node->arc;
	while (arc != NULL) {
		result++;
		if (arc->destination == node->twinNode)
			result++;
		arc = arc->next;
	}

	return result;

}

static Nucleotide getNucleotideInDescriptor(Descriptor * descriptor,
					    Coordinate i)
{
	Descriptor *fourMer = descriptor + i / 4;

	switch (i % 4) {
	case 0:
		return (*fourMer & 3);
	case 1:
		return (*fourMer & 12) >> 2;
	case 2:
		return (*fourMer & 48) >> 4;
	case 3:
		return (*fourMer & 192) >> 6;
	}
	return 0;
}

Nucleotide getNucleotideInNode(Node * node, Coordinate index) {
        return getNucleotideInDescriptor(node->descriptor, index);
}

char *readNode(Node * node)
{
	char *s = callocOrExit(1000000, char);
	char tmpString[100000];
	Descriptor *descriptor = node->descriptor;
	Nucleotide nucleotide;
	Coordinate i;

	sprintf(s, "NODE %d :", node->ID);

	for (i = 0; i < node->length; i++) {
		nucleotide = getNucleotideInDescriptor(descriptor, i);
		switch (nucleotide) {
		case ADENINE:
			tmpString[i] = 'A';
			break;
		case CYTOSINE:
			tmpString[i] = 'C';
			break;
		case GUANINE:
			tmpString[i] = 'G';
			break;
		case THYMINE:
			tmpString[i] = 'T';
			break;
		}
	}

	tmpString[i] = '\0';
	strcat(s, tmpString);

	/*
	   while (arc != NULL) {
	   sprintf(tmpString, " %d(%dx);", arc->destination->ID,
	   arc->multiplicity);
	   strcat(s, tmpString);
	   arc = arc->next;
	   }

	   sprintf(tmpString, " lgth: %d", node->length);
	   strcat(s, tmpString);
	 */

	return s;
}

void displayGraph(Graph * graph)
{
	Node *currentNode;
	IDnum nodeIndex;

	printf("%d sequences\n", graph->sequenceCount);
	printf("%d*2 nodes\n", graph->nodeCount);

	for (nodeIndex = 1; nodeIndex <= graph->nodeCount; nodeIndex++) {
		currentNode = graph->nodes[nodeIndex];
		printf("%s\n", readNode(currentNode));
		printf("%s\n", readNode(currentNode->twinNode));
	}
}

PassageMarker *getMarker(Node * node)
{
	return node->marker;
}

void setMarker(Node * node, PassageMarker * marker)
{
	if (node == NULL)
		return;

	if (marker == NULL) {
		node->marker = NULL;
		node->twinNode->marker = NULL;
		return;
	}

	node->marker = marker;
	setTopOfTheNode(marker);
	node->twinNode->marker = getTwinMarker(marker);
	setTopOfTheNode(getTwinMarker(marker));
}

void setNodeStatus(Node * node, boolean status)
{
	node->status = status;
	node->twinNode->status = status;
}

void setSingleNodeStatus(Node * node, boolean status)
{
	node->status = status;
}

boolean getNodeStatus(Node * node)
{
	if (node == NULL)
		return false;
	return node->status;
}

IDnum getNodeID(Node * node)
{
	if (node == NULL)
		return 0;

	return node->ID;
}

void resetNodeStatus(Graph * graph)
{
	IDnum nodeIndex;
	Node *node;

	for (nodeIndex = 1; nodeIndex <= graph->nodeCount; nodeIndex++) {
		node = graph->nodes[nodeIndex];
		if (node == NULL)
			continue;

		node->status = false;
		node->twinNode->status = false;
	}
}

void resetPassageMarkersStatus(Graph * graph)
{
	IDnum nodeIndex;
	Node *node;
	PassageMarker *marker;

	for (nodeIndex = 1; nodeIndex <= graph->nodeCount; nodeIndex++) {
		node = graph->nodes[nodeIndex];
		if (node == NULL)
			continue;

		for (marker = node->marker; marker != NULL;
		     marker = getNextInNode(marker))
			setPassageMarkerStatus(marker, false);
	}
}

Node *getNodeInGraph(Graph * graph, IDnum nodeID)
{
	if (nodeID == 0)
		return NULL;
	else if (nodeID > 0)
		return graph->nodes[nodeID];
	else if (graph->nodes[-nodeID] == NULL)
		return NULL;
	else
		return graph->nodes[-nodeID]->twinNode;
}

Arc *getArc(Node * node)
{
	return node->arc;
}

Arc *getNextArc(Arc * arc)
{
	return arc->next;
}

IDnum getMultiplicity(Arc * arc)
{
	if (arc == NULL)
		return 0;

	return arc->multiplicity;
}

Node *getOrigin(Arc * arc)
{
	if (arc == NULL)
		return NULL;

	return arc->twinArc->destination->twinNode;
}

Node *getDestination(Arc * arc)
{
	if (arc == NULL)
		return NULL;

	return arc->destination;
}

IDnum markerCount(Node * node)
{
	IDnum count = 0;
	PassageMarker *marker;

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker))
		count++;

	return count;
}

void appendNodeSequence(Node * node, TightString * sequence,
			Coordinate writeIndex)
{
	Coordinate i;
	Nucleotide nucleotide;

	//printf("Getting sequence from node %d of length %d (%d)\n", getNodeID(node), getNodeLength(node), getLength(nodeLabel));

	for (i = 0; i < getNodeLength(node); i++) {
		nucleotide =
		    getNucleotideInDescriptor(node->descriptor, i);
		writeNucleotideAtPosition(nucleotide, i + writeIndex,
					  sequence);
	}
}

static void writeNucleotideInDescriptor(Nucleotide nucleotide,
					Descriptor * descriptor,
					Coordinate i)
{
	Descriptor *fourMer = descriptor + i / 4;
	switch (i % 4) {
	case 3:
		*fourMer &= 63;
		*fourMer += nucleotide << 6;
		return;
	case 2:
		*fourMer &= 207;
		*fourMer += nucleotide << 4;
		return;
	case 1:
		*fourMer &= 243;
		*fourMer += nucleotide << 2;
		return;
	case 0:
		*fourMer &= 252;
		*fourMer += nucleotide;
	}
}

static inline Descriptor *mergeDescriptors(Descriptor * descr,
					   Coordinate destinationLength,
					   Descriptor * copy,
					   Coordinate sourceLength,
					   size_t arrayLength)
{
	Descriptor *readPtr, *writePtr;
	Descriptor readCopy;
	int readOffset, writeOffset;
	Descriptor *new = callocOrExit(arrayLength, Descriptor);
	Coordinate index;

	readPtr = descr;
	readCopy = *readPtr;
	writePtr = new;
	writeOffset = 0;
	for (index = 0; index < destinationLength; index++) {
		(*writePtr) >>= 2;
		(*writePtr) += (readCopy & 3) << 6;
		readCopy >>= 2;

		writeOffset++;
		if (writeOffset == 4) {
			writePtr++;
			readPtr++;
			if (index < destinationLength - 1)
				readCopy = *readPtr;
			writeOffset = 0;
		}
	}

	readPtr = copy;
	readCopy = *readPtr;
	readOffset = 0;
	for (index = 0; index < sourceLength; index++) {
		(*writePtr) >>= 2;
		(*writePtr) += (readCopy & 3) << 6;
		readCopy >>= 2;

		writeOffset++;
		if (writeOffset == 4) {
			writePtr++;
			writeOffset = 0;
		}

		readOffset++;
		if (readOffset == 4) {
			readPtr++;
			if (index < sourceLength - 1)
				readCopy = *readPtr;
			readOffset = 0;
		}
	}

	if (writeOffset != 0) {
		while (writeOffset != 4) {
			(*writePtr) >>= 2;
			writeOffset++;
		}
	}

	return new;
}

static void addBufferToDescriptor(Node * node, Coordinate length)
{
	Descriptor *descr;
	Descriptor *twinDescr;
	Coordinate newLength;
	size_t arrayLength;
	Node *twinNode;
	Coordinate index;
	Descriptor *old_descriptor;

	if (node == NULL)
		return;

	twinNode = node->twinNode;
	descr = node->descriptor;
	twinDescr = twinNode->descriptor;

	// Amendments for empty descriptors
	if (descr == NULL) {
		arrayLength = length / 4;
		if (length % 4 != 0)
			arrayLength++;

		node->descriptor = callocOrExit(arrayLength, Descriptor);
		node->length = length;
		twinNode->descriptor =
		    callocOrExit(arrayLength, Descriptor);
		twinNode->length = length;
		return;
	}

	newLength = node->length + length;
	arrayLength = newLength / 4;
	if (newLength % 4 != 0)
		arrayLength++;

	// Merging forward descriptors
	node->descriptor =
	    reallocOrExit(node->descriptor, arrayLength, Descriptor);

	for (index = node->length; index < newLength; index++)
		writeNucleotideInDescriptor(ADENINE, node->descriptor,
					    index);
	node->length = newLength;

	// Merging reverse descriptors
	old_descriptor = twinNode->descriptor;
	twinNode->descriptor = callocOrExit(arrayLength, Descriptor);
	for (index = 0; index < twinNode->length; index++)
		writeNucleotideInDescriptor(getNucleotideInDescriptor
					    (old_descriptor, index),
					    twinNode->descriptor,
					    index + length);
	for (index = 0; index < length; index++)
		writeNucleotideInDescriptor(THYMINE, twinNode->descriptor,
					    index);
	free(old_descriptor);
	twinNode->length = newLength;
}

void appendDescriptors(Node * destination, Node * source)
{
	Descriptor *copy;
	Descriptor *twinCopy;
	Descriptor *descr;
	Descriptor *twinDescr;
	Coordinate newLength, destinationLength, sourceLength;
	size_t arrayLength;
	Descriptor *new;
	Node *twinDestination;

	if (source == NULL || destination == NULL)
		return;

	twinDestination = destination->twinNode;
	descr = destination->descriptor;
	twinDescr = twinDestination->descriptor;
	copy = source->descriptor;
	twinCopy = source->twinNode->descriptor;

	// Amendments for empty descriptors
	if (getNodeLength(source) == 0)
		return;
	if (getNodeLength(destination) == 0) {
		destination->descriptor = copy;
		twinDestination->descriptor = twinCopy;
		source->descriptor = NULL;
		source->twinNode->descriptor = NULL;
		destination->length = source->length;
		destination->twinNode->length = source->length;
		source->length = 0;
		source->twinNode->length = 0;
		return;
	}

	destinationLength = destination->length;
	sourceLength = source->length;
	newLength = destinationLength + sourceLength;
	arrayLength = newLength / 4;
	if (newLength % 4 != 0)
		arrayLength++;

	// Merging forward descriptors
	new =
	    mergeDescriptors(descr, destinationLength, copy, sourceLength,
			     arrayLength);
	free(descr);
	destination->descriptor = new;
	destination->length = newLength;

	// Merging reverse descriptors
	new =
	    mergeDescriptors(twinCopy, sourceLength, twinDescr,
			     destinationLength, arrayLength);
	free(twinDescr);
	twinDestination->descriptor = new;
	twinDestination->length = newLength;
}

static void catDescriptors(Descriptor * descr, Coordinate destinationLength, Descriptor * copy, Coordinate sourceLength) 
{
	Coordinate index;
	Nucleotide nucleotide;

	for (index = 0; index < sourceLength; index++) {
		nucleotide = getNucleotideInDescriptor(copy, index);
		writeNucleotideInDescriptor(nucleotide, descr, index + destinationLength);
	}
}

static void reverseCatDescriptors(Descriptor * descr, Coordinate destinationLength, Descriptor * copy, Coordinate sourceLength, Coordinate totalLength) 
{
	Coordinate shift = totalLength - destinationLength - sourceLength;
	Coordinate index;
	Nucleotide nucleotide;

	for (index = 0; index < sourceLength; index++) {
		nucleotide = getNucleotideInDescriptor(copy, index);
		writeNucleotideInDescriptor(nucleotide, descr, index + shift);
	}
}

void directlyAppendDescriptors(Node * destination, Node * source, Coordinate totalLength)
{
	Descriptor *copy;
	Descriptor *twinCopy;
	Descriptor *descr;
	Descriptor *twinDescr;
	Coordinate destinationLength, sourceLength;

	if (source == NULL || destination == NULL)
		return;

	descr = destination->descriptor;
	twinDescr = destination->twinNode->descriptor;
	copy = source->descriptor;
	twinCopy = source->twinNode->descriptor;

	// Amendments for empty descriptors
	if (getNodeLength(source) == 0)
		return;

	destinationLength = destination->length;
	sourceLength = source->length;

	// Merging forward descriptors
	catDescriptors(descr, destinationLength, copy, sourceLength);

	// Merging reverse descriptors
	reverseCatDescriptors(twinDescr, destinationLength, twinCopy, sourceLength, totalLength);

	destination->length += source->length;
	destination->twinNode->length += source->length;
}

static void copyDownDescriptor(Descriptor ** writePtr, int *writeOffset,
			       Descriptor * source, Coordinate length)
{
	Descriptor *readPtr = source;
	Descriptor readCopy = *readPtr;
	int readOffset = 0;
	Coordinate index;

	for (index = 0; index < length; index++) {
		(**writePtr) >>= 2;
		(**writePtr) += (readCopy & 3) << 6;
		readCopy >>= 2;

		(*writeOffset)++;
		if (*writeOffset == 4) {
			(*writePtr)++;
			*writeOffset = 0;
		}

		readOffset++;
		if (readOffset == 4) {
			readPtr++;
			if (index < length - 1)
				readCopy = *readPtr;
			readOffset = 0;
		}
	}
}

static void copyDownSequence(Descriptor ** writePtr, int *writeOffset,
			     TightString * sequence, Coordinate start,
			     Coordinate finish, int WORDLENGTH)
{
	boolean forward = (start < finish);
	Coordinate sourceLength = finish - start;
	Coordinate index;
	Nucleotide nucleotide;

	if (!forward)
		sourceLength *= -1;

	for (index = 0; index < sourceLength; index++) {
		if (forward)
			nucleotide =
			    getNucleotide(start + WORDLENGTH - 1 + index,
					  sequence);
		else
			nucleotide =
#ifndef COLOR
			    3 - getNucleotide(start - index - 1, sequence);
#else
			    getNucleotide(start - index - 1, sequence);
#endif

		(**writePtr) >>= 2;
		(**writePtr) += nucleotide << 6;

		(*writeOffset)++;
		if (*writeOffset == 4) {
			(*writePtr)++;
			*writeOffset = 0;
		}
	}
}

static Descriptor *appendSequenceToDescriptor(Descriptor * descr,
					      Coordinate nodeLength,
					      PassageMarker * marker,
					      TightString ** sequences,
					      int WORDLENGTH,
					      size_t arrayLength,
					      boolean downStream)
{
	int writeOffset = 0;
	Descriptor *new = callocOrExit(arrayLength, Descriptor);
	Descriptor *writePtr = new;
	TightString *sequence;
	IDnum sequenceID = getPassageMarkerSequenceID(marker);
	Coordinate start = getPassageMarkerStart(marker);
	Coordinate finish = getPassageMarkerFinish(marker);

	if (sequenceID > 0)
		sequence = sequences[sequenceID - 1];
	else
		sequence = sequences[-sequenceID - 1];

	if (downStream)
		copyDownDescriptor(&writePtr, &writeOffset, descr,
				   nodeLength);

	copyDownSequence(&writePtr, &writeOffset, sequence, start, finish,
			 WORDLENGTH);

	if (!downStream)
		copyDownDescriptor(&writePtr, &writeOffset, descr,
				   nodeLength);

	if (writeOffset != 0) {
		while (writeOffset != 4) {
			(*writePtr) >>= 2;
			writeOffset++;
		}
	}

	return new;
}

void appendSequence(Node * node, TightString ** reads,
		    PassageMarker * guide, Graph * graph)
{
	Descriptor *descr;
	Descriptor *twinDescr;
	Coordinate newLength, nodeLength, sourceLength;
	size_t arrayLength;
	Descriptor *new;
	Node *twinNode;

	if (node == NULL)
		return;

	twinNode = node->twinNode;
	descr = node->descriptor;
	twinDescr = twinNode->descriptor;
	nodeLength = node->length;
	sourceLength = getPassageMarkerLength(guide);

	// Amendments for empty descriptors
	if (sourceLength == 0)
		return;

	newLength = nodeLength + sourceLength;
	arrayLength = newLength / 4;
	if (newLength % 4 != 0)
		arrayLength++;

	// Merging forward descriptors
	new =
	    appendSequenceToDescriptor(descr, nodeLength, guide, reads,
				       getWordLength(graph), arrayLength,
				       true);
	free(descr);
	node->descriptor = new;
	node->length = newLength;

	// Merging reverse descriptors
	new =
	    appendSequenceToDescriptor(twinDescr, nodeLength,
				       getTwinMarker(guide), reads,
				       getWordLength(graph), arrayLength,
				       false);
	free(twinDescr);
	twinNode->descriptor = new;
	twinNode->length = newLength;
}

void setMultiplicity(Arc * arc, IDnum mult)
{
	arc->multiplicity = mult;
	arc->twinArc->multiplicity = mult;
}

// Reshuffles the graph->nodes array to remove NULL pointers
// Beware that node IDs are accordingly reshuffled (all pointers remain valid though)
void renumberNodes(Graph * graph)
{
	IDnum nodeIndex;
	Node *currentNode;
	IDnum counter = 0;
	IDnum nodes = graph->nodeCount;
	IDnum newIndex;

	puts("Renumbering nodes");
	printf("Initial node count %d\n", graph->nodeCount);

	for (nodeIndex = 1; nodeIndex <= nodes; nodeIndex++) {
		currentNode = getNodeInGraph(graph, nodeIndex);

		if (currentNode == NULL)
			counter++;
		else if (counter != 0) {
			newIndex = nodeIndex - counter;
			currentNode->ID = newIndex;
			currentNode->twinNode->ID = -newIndex;
			graph->nodes[newIndex] = currentNode;

			if (graph->nodeReads != NULL) {
				graph->nodeReads[newIndex + nodes] =
				    graph->nodeReads[nodeIndex + nodes];
				graph->nodeReadCounts[newIndex + nodes] =
				    graph->nodeReadCounts[nodeIndex +
							  nodes];

				graph->nodeReads[nodeIndex + nodes] = NULL;
				graph->nodeReadCounts[nodeIndex + nodes] =
				    0;

				graph->nodeReads[-newIndex + nodes] =
				    graph->nodeReads[-nodeIndex + nodes];
				graph->nodeReadCounts[-newIndex + nodes] =
				    graph->nodeReadCounts[-nodeIndex +
							  nodes];

				graph->nodeReads[-nodeIndex + nodes] =
				    NULL;
				graph->nodeReadCounts[-nodeIndex + nodes] =
				    0;
			}

			if (graph->gapMarkers != NULL) {
				graph->gapMarkers[newIndex] =
				    graph->gapMarkers[nodeIndex];
				graph->gapMarkers[nodeIndex] = NULL;
			}
		}
	}

	// Shitfting array to the left
	if (graph->nodeReads != NULL && counter != 0) {
		for (nodeIndex = counter; nodeIndex <= 2 * nodes - counter;
		     nodeIndex++) {
			graph->nodeReads[nodeIndex - counter] =
			    graph->nodeReads[nodeIndex];
			graph->nodeReadCounts[nodeIndex - counter] =
			    graph->nodeReadCounts[nodeIndex];
		}
	}

	// Rellocating node space
	graph->nodeCount -= counter;
	graph->nodes =
	    reallocOrExit(graph->nodes, graph->nodeCount + 1, Node *);

	// Reallocating short read marker arrays
	if (graph->nodeReads != NULL) {
		graph->nodeReads =
		    reallocOrExit(graph->nodeReads,
			    2 * graph->nodeCount +
			     1, ShortReadMarker *);
		graph->nodeReadCounts =
		    reallocOrExit(graph->nodeReadCounts,
			    2 * graph->nodeCount + 1, IDnum);
	}

	// Reallocating gap marker table
	if (graph->gapMarkers != NULL)
		graph->gapMarkers = reallocOrExit(graph->gapMarkers,
					    graph->nodeCount +
					     1, GapMarker *);

	printf("Removed %d null nodes\n", counter);
}

void splitNodeDescriptor(Node * source, Node * target, Coordinate offset)
{
	Coordinate originalLength = source->length;
	Coordinate backLength = originalLength - offset;
	Coordinate index;
	Descriptor *descriptor, *new;
	size_t arrayLength;
	Nucleotide nucleotide;

	source->length = offset;
	source->twinNode->length = offset;

	if (target != NULL) {
		target->length = backLength;
		target->twinNode->length = backLength;
		free(target->descriptor);
		free(target->twinNode->descriptor);
		target->descriptor = NULL;
		target->twinNode->descriptor = NULL;
	}

	if (backLength == 0)
		return;

	descriptor = source->descriptor;

	arrayLength = backLength / 4;
	if (backLength % 4 > 0)
		arrayLength++;

	if (target != NULL) {
		// Target node .. forwards
		new = mallocOrExit(arrayLength, Descriptor);
		target->descriptor = new;
		for (index = 0; index < backLength; index++) {
			nucleotide =
			    getNucleotideInDescriptor(descriptor, index);
			writeNucleotideInDescriptor(nucleotide, new,
						    index);
		}
	}
	// Source node
	for (index = backLength; index < originalLength; index++) {
		nucleotide = getNucleotideInDescriptor(descriptor, index);
		writeNucleotideInDescriptor(nucleotide, descriptor,
					    index - backLength);
	}

	if (target == NULL)
		return;

	// target node other way
	descriptor = source->twinNode->descriptor;
	new = mallocOrExit(arrayLength, Descriptor);
	target->twinNode->descriptor = new;

	for (index = offset; index < originalLength; index++) {
		nucleotide = getNucleotideInDescriptor(descriptor, index);
		writeNucleotideInDescriptor(nucleotide, new,
					    index - offset);
	}
}

void reduceNode(Node * node)
{
	free(node->descriptor);
	node->descriptor = NULL;
	node->length = 0;

	free(node->twinNode->descriptor);
	node->twinNode->descriptor = NULL;
	node->twinNode->length = 0;
}

void checkPassageMarkersStatus(Graph * graph)
{
	IDnum nodeIndex;
	Node *node;
	PassageMarker *marker;

	for (nodeIndex = 1; nodeIndex <= graph->nodeCount; nodeIndex++) {
		node = graph->nodes[nodeIndex];
		if (node == NULL)
			continue;

		for (marker = node->marker; marker != NULL;
		     marker = getNextInNode(marker)) {
			if (getPassageMarkerStatus(marker)) {
				printf("TRUE marker %s\n",
				       readPassageMarker(marker));
				exit(-1);
			}

			if (getNextInSequence(marker) != NULL
			    && getArcBetweenNodes(node,
						  getNode(getNextInSequence
							  (marker)),
						  graph) == NULL) {
				printf
				    ("Missing arc %d -> %d (for %d)\n",
				     getNodeID(node),
				     getNodeID(getNode
					       (getNextInSequence
						(marker))),
				     getPassageMarkerSequenceID(marker));
				abort();
			}
			if (getPreviousInSequence(marker) != NULL
			    &&
			    getArcBetweenNodes(getNode
					       (getPreviousInSequence
						(marker)), node,
					       graph) == NULL) {
				printf
				    ("Missing arc %d -> %d (for %d)\n",
				     getNodeID(getNode
					       (getNextInSequence
						(marker))),
				     getNodeID(node),
				     getPassageMarkerSequenceID(marker));
				abort();
			}
		}
	}
}

void reassessArcMultiplicities(Graph * graph)
{
	IDnum index;
	Node *node, *twin;
	Arc *arc;
	PassageMarker *marker;

	for (index = 1; index <= graph->nodeCount; index++) {
		node = getNodeInGraph(graph, index);

		if (node == NULL)
			continue;

		for (arc = getArc(node); arc != NULL;
		     arc = getNextArc(arc))
			setMultiplicity(arc, 0);
		for (arc = getArc(getTwinNode(node)); arc != NULL;
		     arc = getNextArc(arc))
			setMultiplicity(arc, 0);
	}

	for (index = 1; index <= graph->nodeCount; index++) {
		node = getNodeInGraph(graph, index);

		if (node == NULL)
			continue;

		twin = getTwinNode(node);

		for (marker = getMarker(node); marker != NULL;
		     marker = getNextInNode(marker)) {
			if (getPassageMarkerSequenceID(marker) > 0
			    && !isTerminal(marker)) {
				arc = getArcBetweenNodes(node,
							 getNode
							 (getNextInSequence
							  (marker)),
							 graph);
				if (arc != NULL)
					changeMultiplicity(arc, 1);
			} else if (getPassageMarkerSequenceID(marker) < 0
				   && !isInitial(marker)) {
				arc = getArcBetweenNodes(twin,
							 getTwinNode
							 (getNode
							  (getPreviousInSequence
							   (marker))),
							 graph);
				if (arc != NULL)
					changeMultiplicity(arc, 1);
			}
		}

	}
}

// Allocate memory for an empty graph created with sequenceCount different sequences
Graph *emptyGraph(IDnum sequenceCount, int wordLength)
{
	Graph *newGraph = mallocOrExit(1, Graph);
	newGraph->sequenceCount = sequenceCount;
	newGraph->arcLookupTable = NULL;
	newGraph->nodeReads = NULL;
	newGraph->nodeReadCounts = NULL;
	newGraph->wordLength = wordLength;
	newGraph->gapMarkers = NULL;
	return newGraph;
}

static Descriptor *newPositiveDescriptor(IDnum sequenceID,
					 Coordinate start,
					 Coordinate finish,
					 TightString ** sequences,
					 int WORDLENGTH)
{
	Coordinate index;
	Nucleotide nucleotide;
	TightString *tString = sequences[sequenceID - 1];
	Coordinate length = finish - start;
	Descriptor *res;
	size_t arrayLength = length / 4;

	if (length % 4 > 0)
		arrayLength++;

	res = mallocOrExit(arrayLength, Descriptor);

	for (index = 0; index < length; index++) {
		nucleotide =
		    getNucleotide(start + index + WORDLENGTH - 1, tString);
		writeNucleotideInDescriptor(nucleotide, res, index);
	}

	return res;

}

static Descriptor *newNegativeDescriptor(IDnum sequenceID,
					 Coordinate start,
					 Coordinate finish,
					 TightString ** sequences,
					 int WORDLENGTH)
{
	Coordinate index;
	Nucleotide nucleotide;
	TightString *tString = sequences[-sequenceID - 1];
	Coordinate length = start - finish;
	Descriptor *res;
	size_t arrayLength = length / 4;

	if (length % 4 > 0)
		arrayLength++;

	res = mallocOrExit(arrayLength, Descriptor);

	for (index = 0; index < length; index++) {
		nucleotide = getNucleotide(start - index, tString);
#ifndef COLOR
		writeNucleotideInDescriptor(3 - nucleotide, res, index);
#else
		writeNucleotideInDescriptor(nucleotide, res, index);
#endif
	}

	return res;

}

static Descriptor *newDescriptor(IDnum sequenceID, Coordinate start,
				 Coordinate finish,
				 TightString ** sequences, int WORDLENGTH)
{
	if (sequenceID > 0)
		return newPositiveDescriptor(sequenceID, start, finish,
					     sequences, WORDLENGTH);
	else
		return newNegativeDescriptor(sequenceID, start, finish,
					     sequences, WORDLENGTH);
}

// Constructor
// Memory allocated
Node *newNode(IDnum sequenceID, Coordinate start, Coordinate finish,
	      Coordinate offset, IDnum ID, TightString ** sequences,
	      int WORDLENGTH)
{
	Node *newnd = allocateNode();
	Node *antiNode = allocateNode();
	Category cat;

	newnd->ID = ID;
	newnd->descriptor =
	    newDescriptor(sequenceID, start + offset, finish + offset,
			  sequences, WORDLENGTH);
	newnd->arc = NULL;
	newnd->arcCount = 0;
	newnd->marker = NULL;
	newnd->status = false;
	for (cat = 0; cat < CATEGORIES; cat++) {
		newnd->virtualCoverage[cat] = 0;
		newnd->originalVirtualCoverage[cat] = 0;
	}

	antiNode->ID = -ID;
	antiNode->descriptor =
	    newDescriptor(-sequenceID, finish + offset - 1,
			  start + offset - 1, sequences, WORDLENGTH);
	antiNode->arc = NULL;
	antiNode->arcCount = 0;
	antiNode->marker = NULL;
	antiNode->status = false;
	for (cat = 0; cat < CATEGORIES; cat++) {
		antiNode->virtualCoverage[cat] = 0;
		antiNode->originalVirtualCoverage[cat] = 0;
	}

	newnd->twinNode = antiNode;
	antiNode->twinNode = newnd;

	if (sequenceID > 0) {
		newnd->length = finish - start;
		antiNode->length = finish - start;
	} else {
		newnd->length = start - finish;
		antiNode->length = start - finish;
	}

	return newnd;
}

void allocateNodeSpace(Graph * graph, IDnum nodeCount)
{
	graph->nodes = callocOrExit(nodeCount + 1, Node *);
	graph->nodeCount = nodeCount;
}

void addNodeToGraph(Graph * graph, Node * node)
{
	graph->nodes[node->ID] = node;
}

boolean getUniqueness(Node * node)
{
	return node->uniqueness;
}

void setUniqueness(Node * node, boolean value)
{
	node->uniqueness = value;
	node->twinNode->uniqueness = value;
}

Node *emptyNode()
{
	Node *newnd = allocateNode();
	Node *antiNode = allocateNode();
	Category cat;

	newnd->ID = 0;
	newnd->descriptor = NULL;
	newnd->arc = NULL;
	newnd->arcCount = 0;
	newnd->marker = NULL;
	newnd->length = 0;
	newnd->uniqueness = false;
	for (cat = 0; cat < CATEGORIES; cat++) {
		newnd->virtualCoverage[cat] = 0;
		newnd->originalVirtualCoverage[cat] = 0;
	}

	antiNode->ID = 0;
	antiNode->descriptor = NULL;
	antiNode->arc = NULL;
	antiNode->arcCount = 0;
	antiNode->marker = NULL;
	antiNode->length = 0;
	antiNode->uniqueness = false;
	for (cat = 0; cat < CATEGORIES; cat++) {
		antiNode->virtualCoverage[cat] = 0;
		antiNode->originalVirtualCoverage[cat] = 0;
	}

	newnd->twinNode = antiNode;
	antiNode->twinNode = newnd;

	return newnd;

}

Node *addEmptyNodeToGraph(Graph * graph, IDnum ID)
{
	Node *newnd = emptyNode();

	newnd->ID = ID;
	newnd->twinNode->ID = -ID;

	graph->nodes[ID] = newnd;

	return newnd;

}

void setVirtualCoverage(Node * node, Category category,
			Coordinate coverage)
{
	node->virtualCoverage[category] = coverage;
	node->twinNode->virtualCoverage[category] =
	    node->virtualCoverage[category];
}

void incrementVirtualCoverage(Node * node, Category category,
			      Coordinate coverage)
{
	node->virtualCoverage[category] += coverage;
	node->twinNode->virtualCoverage[category] =
	    node->virtualCoverage[category];

}

Coordinate getVirtualCoverage(Node * node, Category category)
{
	return node->virtualCoverage[category];
}

void setOriginalVirtualCoverage(Node * node, Category category,
				Coordinate coverage)
{
	node->originalVirtualCoverage[category] = coverage;
	node->twinNode->originalVirtualCoverage[category] =
	    node->originalVirtualCoverage[category];
}

void incrementOriginalVirtualCoverage(Node * node, Category category,
				      Coordinate coverage)
{
	node->originalVirtualCoverage[category] += coverage;
	node->twinNode->originalVirtualCoverage[category] =
	    node->originalVirtualCoverage[category];
}

Coordinate getOriginalVirtualCoverage(Node * node, Category category)
{
	return node->originalVirtualCoverage[category];
}

void clipNodeLength(Node * node, Coordinate startClip,
		    Coordinate finishClip)
{
	Descriptor *descriptor;
	Coordinate finalLength =
	    getNodeLength(node) - startClip - finishClip;
	Coordinate index;
	Node *twin = getTwinNode(node);
	Nucleotide nucleotide;

	if (finalLength < 0) {
		puts("Can't clip node that much!!");
		exit(-1);
	}

	if (getNodeLength(node) == 0) {
		puts("Short enough as is");
		exit(-1);
	}
	// One way
	descriptor = node->descriptor;
	for (index = 0; index < finalLength; index++) {
		nucleotide =
		    getNucleotideInDescriptor(descriptor,
					      index + startClip);
		writeNucleotideInDescriptor(nucleotide, descriptor, index);
	}

	// Same thing in the other direction
	descriptor = twin->descriptor;
	for (index = 0; index < finalLength; index++) {
		nucleotide =
		    getNucleotideInDescriptor(descriptor,
					      index + finishClip);
		writeNucleotideInDescriptor(nucleotide, descriptor, index);
	}

	// Length
	node->length = finalLength;
	node->twinNode->length = node->length;
}

boolean hasSingleArc(Node * node)
{
	return node->arcCount == 1;
}

void activateArcLookupTable(Graph * graph)
{
	IDnum index;
	Node *node;
	Arc *arc;
	IDnum nodes = graph->nodeCount;
	IDnum twinOriginID, destinationID, hash;
	Arc **table;

	puts("Activating arc lookup table");

	graph->arcLookupTable = callocOrExit(6 * nodes + 1, Arc *);

	table = graph->arcLookupTable;

	for (index = -nodes; index <= nodes; index++) {
		if (index == 0)
			continue;

		node = getNodeInGraph(graph, index);
		if (node == 0)
			continue;

		for (arc = getArc(node); arc != NULL;
		     arc = getNextArc(arc)) {
			twinOriginID = arc->twinArc->destination->ID;
			destinationID = arc->destination->ID;
			hash =
			    3 * nodes - 2 * twinOriginID + destinationID;
			arc->nextInLookupTable = table[hash];
			table[hash] = arc;
		}
	}

	puts("Done activating arc lookup table");
}

void deactivateArcLookupTable(Graph * graph)
{
	free(graph->arcLookupTable);
	graph->arcLookupTable = NULL;
}

static void exportNode(FILE * outfile, Node * node, void *withSequence)
{
	Coordinate index;
	Nucleotide nucleotide;
	Category cat;

	if (node == NULL)
		return;

	fprintf(outfile, "NODE\t%ld\t%lld", (long) node->ID, (long long) node->length);
	for (cat = 0; cat < CATEGORIES; cat++)
		fprintf(outfile, "\t%lld\t%lld", (long long) node->virtualCoverage[cat],
			(long long) node->originalVirtualCoverage[cat]);
	fprintf(outfile, "\n");

	if (withSequence == NULL)
		return;

	for (index = 0; index < node->length; index++) {
		nucleotide =
		    getNucleotideInDescriptor(node->descriptor, index);
		switch (nucleotide) {
		case ADENINE:
			fprintf(outfile, "A");
			break;
		case CYTOSINE:
			fprintf(outfile, "C");
			break;
		case GUANINE:
			fprintf(outfile, "G");
			break;
		case THYMINE:
			fprintf(outfile, "T");
			break;
		}
	}
	fprintf(outfile, "\n");

	for (index = 0; index < node->length; index++) {
		nucleotide =
		    getNucleotideInDescriptor(node->twinNode->descriptor,
					      index);
		switch (nucleotide) {
		case ADENINE:
			fprintf(outfile, "A");
			break;
		case CYTOSINE:
			fprintf(outfile, "C");
			break;
		case GUANINE:
			fprintf(outfile, "G");
			break;
		case THYMINE:
			fprintf(outfile, "T");
			break;
		}
	}
	fprintf(outfile, "\n");
}

static void exportArc(FILE * outfile, Arc * arc)
{
	IDnum originID, destinationID;
	IDnum absOriginID, absDestinationID;

	if (arc == NULL)
		return;

	absOriginID = originID = -arc->twinArc->destination->ID;
	absDestinationID = destinationID = arc->destination->ID;

	if (absOriginID < 0)
		absOriginID = -absOriginID;
	if (absDestinationID < 0)
		absDestinationID = -absDestinationID;

	if (absDestinationID < absOriginID)
		return;

	if (originID == destinationID && originID < 0)
		return;

	fprintf(outfile, "ARC\t%d\t%d\t%d\n", originID, destinationID,
		arc->multiplicity);
}

// Merges two lists of annotations in order of increasing position (used in mergeSort mainly)
static Arc *mergeArcLists(Arc * left, Arc * right)
{
	Arc *mergedList = NULL;
	Arc *tail = NULL;

	// Choose first element:
	if (left->destination->ID <= right->destination->ID) {
		mergedList = left;
		tail = left;
		left = left->next;
	} else {
		mergedList = right;
		tail = right;
		right = right->next;
	}

	// Iterate while both lists are still non empty
	while (left != NULL && right != NULL) {
		if (left->destination->ID <= right->destination->ID) {
			tail->next = left;
			left->previous = tail;
			left = left->next;
		} else {
			tail->next = right;
			right->previous = tail;
			right = right->next;
		}

		tail = tail->next;
	}

	// Concatenate the remaining list at the end of the merged list
	if (left != NULL) {
		tail->next = left;
		left->previous = tail;
	}

	if (right != NULL) {
		tail->next = right;
		right->previous = tail;
	}

	return mergedList;
}

static void arcMergeSort(Arc ** arcPtr, IDnum count)
{

	IDnum half = count / 2;
	Arc *left = *arcPtr;
	Arc *ptr = left;
	Arc *right;
	IDnum index;

	if (count == 0 || count == 1)
		return;

	if (count == 2) {
		if ((*arcPtr)->destination->ID >
		    (*arcPtr)->next->destination->ID) {
			(*arcPtr)->next->next = *arcPtr;
			(*arcPtr)->previous = (*arcPtr)->next;
			*arcPtr = (*arcPtr)->next;
			(*arcPtr)->next->next = NULL;
			(*arcPtr)->previous = NULL;
		}
		return;
	}

	for (index = 0; index < half - 1; index++) {
		ptr = ptr->next;
		if (ptr == NULL)
			return;
	}

	right = ptr->next;
	ptr->next = NULL;
	right->previous = NULL;

	arcMergeSort(&left, half);
	arcMergeSort(&right, count - half);
	*arcPtr = mergeArcLists(left, right);
}

static void sortNodeArcs(Node * node)
{
	Arc *arc;
	IDnum count = 0;

	for (arc = getArc(node); arc != NULL; arc = getNextArc(arc))
		count++;

	if (count == 0)
		return;

	arc = getArc(node);
	arcMergeSort(&arc, count);

	node->arc = arc;
}

// Merges two lists of annotations in order of increasing position (used in mergeSort mainly)
static GapMarker *mergeGapMarkerLists(GapMarker * left, GapMarker * right)
{
	GapMarker *mergedList = NULL;
	GapMarker *tail = NULL;

	// Choose first element:
	if (left->position <= right->position) {
		mergedList = left;
		tail = left;
		left = left->next;
	} else {
		mergedList = right;
		tail = right;
		right = right->next;
	}

	// Iterate while both lists are still non empty
	while (left != NULL && right != NULL) {
		if (left->position <= right->position) {
			tail->next = left;
			left = left->next;
		} else {
			tail->next = right;
			right = right->next;
		}

		tail = tail->next;
	}

	// Concatenate the remaining list at the end of the merged list
	if (left != NULL)
		tail->next = left;

	if (right != NULL)
		tail->next = right;

	return mergedList;
}

static void gapMergeSort(GapMarker ** gapPtr, IDnum count)
{

	IDnum half = count / 2;
	GapMarker *left = *gapPtr;
	GapMarker *ptr = left;
	GapMarker *right;
	IDnum index;

	if (count == 0 || count == 1)
		return;

	if (count == 2) {
		if ((*gapPtr)->position > (*gapPtr)->next->position) {
			(*gapPtr)->next->next = *gapPtr;
			*gapPtr = (*gapPtr)->next;
			(*gapPtr)->next->next = NULL;
		}
		return;
	}

	for (index = 0; index < half - 1; index++) {
		ptr = ptr->next;
		if (ptr == NULL)
			return;
	}

	right = ptr->next;
	ptr->next = NULL;

	gapMergeSort(&left, half);
	gapMergeSort(&right, count - half);
	*gapPtr = mergeGapMarkerLists(left, right);
}

static void sortNodeGapMarkers(Node * node, Graph * graph)
{
	GapMarker *gap;
	IDnum count = 0;
	IDnum ID = getNodeID(node);

	if (ID < 0)
		ID = -ID;

	for (gap = graph->gapMarkers[ID]; gap != NULL; gap = gap->next)
		count++;

	if (count == 0)
		return;

	gap = graph->gapMarkers[ID];
	gapMergeSort(&gap, count);

	graph->gapMarkers[ID] = gap;
}

void sortGapMarkers(Graph * graph)
{
	IDnum index;
	Node *node;

	if (graph->gapMarkers == NULL)
		return;

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node)
			sortNodeGapMarkers(node, graph);
	}
}

void exportGraph(char *filename, Graph * graph, TightString ** sequences)
{
	IDnum index;
	FILE *outfile;
	Node *node;
	Arc *arc;
	PassageMarker *marker;
	ShortReadMarker *reads;
	IDnum readCount, readIndex;

	if (graph == NULL) {
		return;
	}

	outfile = fopen(filename, "w");
	if (outfile == NULL) {
		puts("Couldn't open file, sorry");
		return;
	} else
		printf("Writing into graph file %s...\n", filename);

	// General data
	fprintf(outfile, "%d\t%d\t%i\n", graph->nodeCount,
		graph->sequenceCount, graph->wordLength);

	// Node info
	for (index = 1; index <= graph->nodeCount; index++) {
		node = getNodeInGraph(graph, index);
		exportNode(outfile, node, (void *) sequences);
	}

	// Arc info
	for (index = 1; index <= graph->nodeCount; index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL)
			continue;

		sortNodeArcs(node);
		sortNodeArcs(getTwinNode(node));

		for (arc = node->arc; arc != NULL; arc = arc->next)
			exportArc(outfile, arc);
		for (arc = node->twinNode->arc; arc != NULL;
		     arc = arc->next)
			exportArc(outfile, arc);
	}

	// Sequence info
	for (index = 1; index <= graph->nodeCount; index++) {
		node = getNodeInGraph(graph, index);
		if (node == NULL)
			continue;
		for (marker = node->marker; marker != NULL;
		     marker = getNextInNode(marker))
			exportMarker(outfile, marker, sequences,
				     graph->wordLength);
	}

	// Node reads
	if (readStartsAreActivated(graph)) {
		for (index = 0; index <= graph->nodeCount * 2; index++) {
			readCount = graph->nodeReadCounts[index];
			if (readCount == 0)
				continue;

			fprintf(outfile, "NR\t%d\t%d\n",
				index - graph->nodeCount, readCount);

			reads = graph->nodeReads[index];
			for (readIndex = 0; readIndex < readCount;
			     readIndex++)
				fprintf(outfile, "%ld\t%lld\t%d\n",
					(long) reads[readIndex].readID,
					(long long) reads[readIndex].position,
					(int) reads[readIndex].offset);
		}
	}

	fclose(outfile);
}

Graph *importGraph(char *filename)
{
	FILE *file = fopen(filename, "r");
	const int maxline = MAXLINE;
	char line[MAXLINE];
	Graph *graph;
	Coordinate coverage, originalCoverage;
	IDnum nodeCounter, sequenceCount;
	Node *node, *twin;
	Arc *arc;
	IDnum originID, destinationID, multiplicity;
	PassageMarker *newMarker, *marker;
	IDnum nodeID, seqID;
	Coordinate index;
	Coordinate start, finish;
	Coordinate startOffset, finishOffset;
	boolean finished = false;
	size_t arrayLength;
	IDnum readCount;
	ShortReadMarker *array;
	int wordLength, sCount;
	ShortLength length;
	Category cat;
	long long_var, long_var2, long_var3;
	long long longlong_var, longlong_var2, longlong_var3, longlong_var4;
	short short_var;
	char c;

	if (file == NULL) 
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	printf("Reading graph file %s\n", filename);

	// First  line
	if (!fgets(line, maxline, file))
		exitErrorf(EXIT_FAILURE, true, "Graph file incomplete");
	sscanf(line, "%ld\t%ld\t%i\n", &long_var, &long_var2,
	       &wordLength);
	nodeCounter = (IDnum) long_var;
	sequenceCount = (IDnum) long_var2;
	graph = emptyGraph(sequenceCount, wordLength);
	resetWordFilter(wordLength);
	allocateNodeSpace(graph, nodeCounter);
	printf("Graph has %ld nodes and %ld sequences\n", (long) nodeCounter,
	       (long) sequenceCount);

	// Read nodes
	if (!fgets(line, maxline, file))
		exitErrorf(EXIT_FAILURE, true, "Graph file incomplete");
	while (!finished && strncmp(line, "NODE", 4) == 0) {
		strtok(line, "\t\n");
		sscanf(strtok(NULL, "\t\n"), "%ld", &long_var);
		nodeID = (IDnum) long_var;
		node = addEmptyNodeToGraph(graph, nodeID);
		sscanf(strtok(NULL, "\t\n"), "%lld", &longlong_var);
		node->length = (Coordinate) longlong_var;
		for (cat = 0; cat < CATEGORIES; cat++) {
			sscanf(strtok(NULL, "\t\n"), "%lld", &longlong_var);
			coverage = (Coordinate) longlong_var;
			setVirtualCoverage(node, cat, coverage);
			sscanf(strtok(NULL, "\t\n"), "%lld",
			       &longlong_var);
			originalCoverage = (Coordinate) longlong_var;
			setOriginalVirtualCoverage(node, cat,
						   originalCoverage);
		}

		arrayLength = node->length / 4;
		if (node->length % 4 > 0)
			arrayLength++;
		node->descriptor =
		    callocOrExit(arrayLength, Descriptor);

		index = 0;
		while ((c = fgetc(file)) != '\n' && c != EOF) {
			if (c == 'A')
				writeNucleotideInDescriptor(ADENINE,
							    node->
							    descriptor,
							    index++);
			else if (c == 'C')
				writeNucleotideInDescriptor(CYTOSINE,
							    node->
							    descriptor,
							    index++);
			else if (c == 'G')
				writeNucleotideInDescriptor(GUANINE,
							    node->
							    descriptor,
							    index++);
			else if (c == 'T')
				writeNucleotideInDescriptor(THYMINE,
							    node->
							    descriptor,
							    index++);
		}

		twin = node->twinNode;
		twin->length = node->length;
		twin->descriptor =
		    callocOrExit(arrayLength, Descriptor);
		index = 0;
		while ((c = fgetc(file)) != '\n' && c != EOF) {
			if (c == 'A')
				writeNucleotideInDescriptor(ADENINE,
							    twin->
							    descriptor,
							    index++);
			else if (c == 'C')
				writeNucleotideInDescriptor(CYTOSINE,
							    twin->
							    descriptor,
							    index++);
			else if (c == 'G')
				writeNucleotideInDescriptor(GUANINE,
							    twin->
							    descriptor,
							    index++);
			else if (c == 'T')
				writeNucleotideInDescriptor(THYMINE,
							    twin->
							    descriptor,
							    index++);
		}

		if (fgets(line, maxline, file) == NULL)
			finished = true;
	}

	// Read arcs
	while (!finished && line[0] == 'A') {
		sscanf(line, "ARC\t%ld\t%ld\t%ld\n", &long_var,
		       &long_var2, &long_var3);
		originID = (IDnum) long_var;
		destinationID = (IDnum) long_var2;
		multiplicity = (IDnum) long_var3;
		arc =
		    createArc(getNodeInGraph(graph, originID),
			      getNodeInGraph(graph, destinationID), graph);
		setMultiplicity(arc, multiplicity);
		if (fgets(line, maxline, file) == NULL)
			finished = true;
	}

	// Read sequences
	while (!finished && line[0] != 'N') {
		sscanf(line, "SEQ\t%ld\n", &long_var);
		seqID = (IDnum) long_var;
		marker = NULL;
		if (!fgets(line, maxline, file))
			exitErrorf(EXIT_FAILURE, true, "Graph file incomplete");

		while (!finished && line[0] != 'N' && line[0] != 'S') {
			sCount =
			    sscanf(line, "%ld\t%lld\t%lld\t%lld\t%lld\n",
				   &long_var, &longlong_var, &longlong_var2, &longlong_var3,
				   &longlong_var4);
			nodeID = (IDnum) long_var;
			startOffset = (Coordinate) longlong_var;
			start = (Coordinate) longlong_var2;
			finish = (Coordinate) longlong_var3;
			finishOffset = (Coordinate) longlong_var4;
			if (sCount != 5) {
				printf
				    ("ERROR: reading in graph - only %d items read for line '%s'",
				     sCount, line);
				exit(1);
			}
			newMarker =
			    newPassageMarker(seqID, start, finish,
					     startOffset, finishOffset);
			transposePassageMarker(newMarker,
					       getNodeInGraph(graph,
							      nodeID));
			connectPassageMarkers(marker, newMarker, graph);
			marker = newMarker;
			if (fgets(line, maxline, file) == NULL)
				finished = true;
		}
	}

	// Node reads
	while (!finished) {
		sscanf(line, "NR\t%ld\t%ld\n", &long_var, &long_var2);
		nodeID = (IDnum) long_var;
		readCount = (IDnum) long_var2;
		if (!readStartsAreActivated(graph))
			activateReadStarts(graph);

		graph->nodeReadCounts[nodeID + graph->nodeCount] =
		    readCount;
		array = mallocOrExit(readCount, ShortReadMarker);
		graph->nodeReads[nodeID + graph->nodeCount] = array;

		readCount = 0;
		if (!fgets(line, maxline, file))
			exitErrorf(EXIT_FAILURE, true, "Graph file incomplete");
		while (!finished && line[0] != 'N') {
			sscanf(line, "%ld\t%lld\t%hd\n", &long_var,
			       &longlong_var, &short_var);
			seqID = (IDnum) long_var;
			startOffset = (Coordinate) longlong_var;
			length = (ShortLength) short_var;
			array[readCount].readID = seqID;
			array[readCount].position = startOffset;
			array[readCount].offset = length;
			readCount++;
			if (fgets(line, maxline, file) == NULL)
				finished = true;
		}
	}

	//printf("New graph has %d nodes\n", graph->nodeCount);

	fclose(file);
	//puts("Done, exiting");
	return graph;
}

Graph *importSimplifiedGraph(char *filename)
{
	FILE *file = fopen(filename, "r");
	const int maxline = MAXLINE;
	char line[MAXLINE];
	Graph *graph;
	Coordinate coverage, originalCoverage;
	IDnum nodeCounter, sequenceCount;
	Node *node, *twin;
	PassageMarker *newMarker, *marker;
	IDnum nodeID, seqID;
	Coordinate index;
	Coordinate start, finish;
	Coordinate startOffset, finishOffset;
	boolean finished = false;
	size_t arrayLength;
	IDnum readCount;
	ShortReadMarker *array = NULL;
	int wordLength, sCount;
	ShortLength length;
	Category cat;
	long long_var, long_var2;
	long long longlong_var, longlong_var2, longlong_var3, longlong_var4;
	short short_var;
	char c;

	if (file == NULL) 
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", filename);

	printf("Reading graph file %s\n", filename);

	// First  line
	if (!fgets(line, maxline, file))
		exitErrorf(EXIT_FAILURE, true, "Graph file incomplete");
	sscanf(line, "%ld\t%ld\t%i\n", &long_var, &long_var2,
	       &wordLength);
	nodeCounter = (IDnum) long_var;
	sequenceCount = (IDnum) long_var2;
	graph = emptyGraph(sequenceCount, wordLength);
	resetWordFilter(wordLength);
	allocateNodeSpace(graph, nodeCounter);
	printf("Graph has %ld nodes and %ld sequences\n", (long) nodeCounter,
	       (long) sequenceCount);

	// Read nodes
	if (!fgets(line, maxline, file))
		exitErrorf(EXIT_FAILURE, true, "Graph file incomplete");
	while (strncmp(line, "NODE", 4) == 0) {
		strtok(line, "\t\n");
		sscanf(strtok(NULL, "\t\n"), "%ld", &long_var);
		nodeID = (IDnum) long_var;
		sscanf(strtok(NULL, "\t\n"), "%lld", &longlong_var);

		if (longlong_var < 50) {
			if (fgets(line, maxline, file) == NULL)
				finished = true;
			if (fgets(line, maxline, file) == NULL)
				finished = true;
			if (fgets(line, maxline, file) == NULL)
				finished = true;
			continue;
		}
		
		node = addEmptyNodeToGraph(graph, nodeID);
		node->length = (Coordinate) longlong_var;
		for (cat = 0; cat < CATEGORIES; cat++) {
			sscanf(strtok(NULL, "\t\n"), "%lld", &longlong_var);
			coverage = (Coordinate) longlong_var;
			setVirtualCoverage(node, cat, coverage);
			sscanf(strtok(NULL, "\t\n"), "%lld",
			       &longlong_var);
			originalCoverage = (Coordinate) longlong_var;
			setOriginalVirtualCoverage(node, cat,
						   originalCoverage);
		}

		arrayLength = node->length / 4;
		if (node->length % 4 > 0)
			arrayLength++;
		node->descriptor =
		    callocOrExit(arrayLength, Descriptor);

		index = 0;
		while ((c = fgetc(file)) != '\n' && c != EOF) {
			if (c == 'A')
				writeNucleotideInDescriptor(ADENINE,
							    node->
							    descriptor,
							    index++);
			else if (c == 'C')
				writeNucleotideInDescriptor(CYTOSINE,
							    node->
							    descriptor,
							    index++);
			else if (c == 'G')
				writeNucleotideInDescriptor(GUANINE,
							    node->
							    descriptor,
							    index++);
			else if (c == 'T')
				writeNucleotideInDescriptor(THYMINE,
							    node->
							    descriptor,
							    index++);
		}

		twin = node->twinNode;
		twin->length = node->length;
		twin->descriptor =
		    callocOrExit(arrayLength, Descriptor);
		index = 0;
		while ((c = fgetc(file)) != '\n' && c != EOF) {
			if (c == 'A')
				writeNucleotideInDescriptor(ADENINE,
							    twin->
							    descriptor,
							    index++);
			else if (c == 'C')
				writeNucleotideInDescriptor(CYTOSINE,
							    twin->
							    descriptor,
							    index++);
			else if (c == 'G')
				writeNucleotideInDescriptor(GUANINE,
							    twin->
							    descriptor,
							    index++);
			else if (c == 'T')
				writeNucleotideInDescriptor(THYMINE,
							    twin->
							    descriptor,
							    index++);
		}

		if (fgets(line, maxline, file) == NULL)
			finished = true;
	}

	// Read arcs
	while (!finished && line[0] == 'A')
		if (fgets(line, maxline, file) == NULL)
			finished = true;

	// Read sequences
	while (!finished && line[0] != 'N') {
		sscanf(line, "SEQ\t%ld\n", &long_var);
		seqID = (IDnum) long_var;
		marker = NULL;
		if (!fgets(line, maxline, file))
			exitErrorf(EXIT_FAILURE, true, "Graph file incomplete");

		while (!finished && line[0] != 'N' && line[0] != 'S') {
			sCount =
			    sscanf(line, "%ld\t%lld\t%lld\t%lld\t%lld\n",
				   &long_var, &longlong_var, &longlong_var2, &longlong_var3,
				   &longlong_var4);
			nodeID = (IDnum) long_var;
			startOffset = (Coordinate) longlong_var;
			start = (Coordinate) longlong_var2;
			finish = (Coordinate) longlong_var3;
			finishOffset = (Coordinate) longlong_var4;
			if (sCount != 5) {
				printf
				    ("ERROR: reading in graph - only %d items read for line '%s'",
				     sCount, line);
				abort();
				exit(1);
			}
			if (getNodeInGraph(graph, nodeID)) {
				newMarker =
				    newPassageMarker(seqID, start, finish,
						     startOffset, finishOffset);
				transposePassageMarker(newMarker,
						       getNodeInGraph(graph,
								      nodeID));
				connectPassageMarkers(marker, newMarker, graph);
				marker = newMarker;
			}
			if (fgets(line, maxline, file) == NULL)
				finished = true;
		}
	}

	// Node reads
	while (!finished) {
		sscanf(line, "NR\t%ld\t%ld\n", &long_var, &long_var2);
		nodeID = (IDnum) long_var;
		readCount = (IDnum) long_var2;
		if (!readStartsAreActivated(graph))
			activateReadStarts(graph);

		if (getNodeInGraph(graph, nodeID)) {
			graph->nodeReadCounts[nodeID + graph->nodeCount] =
			    readCount;
			array = mallocOrExit(readCount, ShortReadMarker);
			graph->nodeReads[nodeID + graph->nodeCount] = array;
		}

		readCount = 0;
		if (!fgets(line, maxline, file))
			exitErrorf(EXIT_FAILURE, true, "Graph file incomplete");
		while (!finished && line[0] != 'N') {
			if (getNodeInGraph(graph, nodeID)) {
				sscanf(line, "%ld\t%lld\t%hd\n", &long_var,
				       &longlong_var, &short_var);
				seqID = (IDnum) long_var;
				startOffset = (Coordinate) longlong_var;
				length = (ShortLength) short_var;
				array[readCount].readID = seqID;
				array[readCount].position = startOffset;
				array[readCount].offset = length;
				readCount++;
			}
			if (fgets(line, maxline, file) == NULL)
				finished = true;
		}
	}

	//printf("New graph has %d nodes\n", graph->nodeCount);

	fclose(file);
	//puts("Done, exiting");
	renumberNodes(graph);
	return graph;
}

Graph *readPreGraphFile(char *preGraphFilename, boolean * double_strand)
{
	FILE *file = fopen(preGraphFilename, "r");
	const int maxline = MAXLINE;
	char line[MAXLINE];

	Graph *graph;
	IDnum nodeCounter, sequenceCount;

	Node *node, *twin;
	IDnum nodeID = 0;
	Coordinate index, nodeLength;
	char c;
	int wordLength, wordShift;
	size_t arrayLength;
	short short_var;
	long long_var, long_var2;
	long long longlong_var;

	if (file == NULL)
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", preGraphFilename);

	printf("Reading pre-graph file %s\n", preGraphFilename);

	// First  line
	if (!fgets(line, maxline, file))
		exitErrorf(EXIT_FAILURE, true, "PreGraph file incomplete");
	sscanf(line, "%ld\t%ld\t%i\t%hi\n", &long_var, &long_var2,
	       &wordLength, &short_var);
	nodeCounter = (IDnum) long_var;
	sequenceCount = (IDnum) long_var2;
	*double_strand = (boolean) short_var;
	wordShift = wordLength - 1;
	graph = emptyGraph(sequenceCount, wordLength);
	resetWordFilter(wordLength);
	allocateNodeSpace(graph, nodeCounter);
	printf("Graph has %ld nodes and %ld sequences\n", (long) nodeCounter,
	       (long) sequenceCount);

	// Read nodes
	if (!fgets(line, maxline, file))
		exitErrorf(EXIT_FAILURE, true, "PreGraph file incomplete");
	while (line[0] == 'N') {
		nodeID++;
		node = addEmptyNodeToGraph(graph, nodeID);

		sscanf(line, "%*s\t%*i\t%lli\n", &longlong_var);
		node->length = (Coordinate) longlong_var;
		nodeLength = node->length;
		arrayLength = node->length / 4;
		if (node->length % 4 > 0)
			arrayLength++;
		node->descriptor =
		    callocOrExit(arrayLength, Descriptor);

		twin = node->twinNode;
		twin->length = nodeLength;
		twin->descriptor =
		    callocOrExit(arrayLength, Descriptor);


		index = 0;
		while ((c = getc(file)) != '\n') {
			if (c == 'A') {
				if (index - wordShift >= 0)
					writeNucleotideInDescriptor(ADENINE,
								    node->
								    descriptor,
								    index - wordShift);
				if (nodeLength - index - 1 >= 0) {
#ifndef COLOR
					writeNucleotideInDescriptor(THYMINE,
								    twin->
								    descriptor,
								    nodeLength - index - 1);
#else
					writeNucleotideInDescriptor(ADENINE,
								    twin->
								    descriptor,
								    nodeLength - index - 1);
#endif
				}
			} else if (c == 'C') {
				if (index - wordShift >= 0)
					writeNucleotideInDescriptor(CYTOSINE,
								    node->
								    descriptor,
								    index - wordShift);
				if (nodeLength - index - 1 >= 0) {
#ifndef COLOR
					writeNucleotideInDescriptor(GUANINE,
								    twin->
								    descriptor,
								    nodeLength - index - 1);
#else
					writeNucleotideInDescriptor(CYTOSINE,
								    twin->
								    descriptor,
								    nodeLength - index - 1);
#endif
				}
			} else if (c == 'G') {
				if (index - wordShift >= 0)
					writeNucleotideInDescriptor(GUANINE,
								    node->
								    descriptor,
								    index - wordShift);
				if (nodeLength - index - 1 >= 0) {
#ifndef COLOR
					writeNucleotideInDescriptor(CYTOSINE,
								    twin->
								    descriptor,
								    nodeLength - index - 1);
#else
					writeNucleotideInDescriptor(GUANINE,
								    twin->
								    descriptor,
								    nodeLength - index - 1);
#endif
				}
			} else if (c == 'T') {
				if (index - wordShift >= 0)
					writeNucleotideInDescriptor(THYMINE,
								    node->
								    descriptor,
								    index - wordShift);
				if (nodeLength - index - 1 >= 0) {
#ifndef COLOR
					writeNucleotideInDescriptor(ADENINE,
								    twin->
								    descriptor,
								    nodeLength - index - 1);
#else
					writeNucleotideInDescriptor(THYMINE,
								    twin->
								    descriptor,
								    nodeLength - index - 1);
#endif
				}
			}
			
			index++;
		}

		if (fgets(line, maxline, file) == NULL) {
			fclose(file);
			return graph;
		}
	}

	fclose(file);
	return graph;
}

// Prints out the information relative to the topology of a node into a new file
// Internal to exportDOTGraph()
void DOTNode(Node * node, FILE * outfile)
{
	IDnum ID;
	Arc *arc;
	Node *otherNode;

	ID = node->ID;
	if (ID < 0)
		return;

	fprintf(outfile, "\t%d [label=\"<left>|%d|<right>\"]\n", ID, ID);

	for (arc = node->arc; arc != NULL; arc = arc->next) {
		otherNode = arc->destination;
		if (!(otherNode->ID >= ID || otherNode->ID <= -ID)) {
			continue;
		}

		if (otherNode->ID > 0)
			fprintf(outfile, "\t%d:right -> %d:left\n", ID,
				otherNode->ID);
		else
			fprintf(outfile, "\t%d:right -> %d:right\n", ID,
				-otherNode->ID);
	}

	for (arc = node->twinNode->arc; arc != NULL; arc = arc->next) {
		otherNode = arc->destination;
		if (!(otherNode->ID >= ID || otherNode->ID <= -ID)) {
			continue;
		}

		if (otherNode->ID > 0)
			fprintf(outfile, "\t%d:left -> %d:left\n", ID,
				otherNode->ID);
		else
			fprintf(outfile, "\t%d:left -> %d:right\n", ID,
				-otherNode->ID);
	}
}

// Exports the topology of a graph into a new file (designated by its filename)
void exportDOTGraph(char *filename, Graph * graph)
{
	IDnum nodeIndex;
	Node *currentNode;

	FILE *outfile = fopen(filename, "w");
	if (outfile == NULL) {
		puts("Couldn't open file, sorry");
		return;
	} else
		puts("Writing into file...");

	fprintf(outfile, "digraph G {\n");
	fprintf(outfile, "\tRANKDIR=LR\n");
	fprintf(outfile, "\tnode [shape=record]\n");

	for (nodeIndex = 1; nodeIndex <= graph->nodeCount; nodeIndex++) {
		currentNode = getNodeInGraph(graph, nodeIndex);
		DOTNode(currentNode, outfile);
	}

	fprintf(outfile, "}\n");
	fclose(outfile);
}

TightString *expandNode(Node * node, int WORDLENGTH)
{
	Nucleotide nucleotide;
	Coordinate index;
	TightString *tString =
	    newTightString(node->length + WORDLENGTH - 1);
	Node *twin = node->twinNode;
	Coordinate length = node->length;

	for (index = 0; index < WORDLENGTH; index++) {
		nucleotide =
		    getNucleotideInDescriptor(twin->descriptor,
					      length - index - 1);
#ifndef COLOR
		writeNucleotideAtPosition(3 - nucleotide, index, tString);
#else
		writeNucleotideAtPosition(nucleotide, index, tString);
#endif
	}

	for (index = 1; index < node->length; index++) {
		nucleotide =
		    getNucleotideInDescriptor(node->descriptor, index);
		writeNucleotideAtPosition(nucleotide,
					  index + WORDLENGTH - 1, tString);
	}

	return tString;
}

char *expandNodeFragment(Node * node, Coordinate contigStart,
			 Coordinate contigFinish, int wordLength)
{
	Nucleotide nucleotide;
	Coordinate index;
	Node *twin = node->twinNode;
	Coordinate length = contigFinish - contigStart;
	int wordShift = wordLength - 1;
	char *string;

	if (length >= wordShift) {
		string = callocOrExit(length + wordLength, char);

		for (index = 0; index < wordShift; index++) {
			nucleotide =
			    getNucleotideInDescriptor(twin->descriptor,
						      twin->length - contigStart -
						      index - 1);
	#ifndef COLOR
			nucleotide = 3 - nucleotide;
	#endif

			switch (nucleotide) {
			case ADENINE:
				string[index] = 'A';
				break;
			case CYTOSINE:
				string[index] = 'C';
				break;
			case GUANINE:
				string[index] = 'G';
				break;
			case THYMINE:
				string[index] = 'T';
				break;
			}

		}

		for (index = 0; index < length; index++) {
			nucleotide =
			    getNucleotideInDescriptor(node->descriptor,
						      contigStart + index);
			switch (nucleotide) {
			case ADENINE:
				string[index + wordShift] = 'A';
				break;
			case CYTOSINE:
				string[index + wordShift] = 'C';
				break;
			case GUANINE:
				string[index + wordShift] = 'G';
				break;
			case THYMINE:
				string[index + wordShift] = 'T';
				break;
			}
		}

		string[length + wordShift] = '\0';
	} else {
		string = callocOrExit(length + 1, char);

		for (index = 0; index < length; index++) {
			nucleotide =
			    getNucleotideInDescriptor(node->descriptor, contigStart + index);
			switch (nucleotide) {
			case ADENINE:
				string[index] = 'A';
				break;
			case CYTOSINE:
				string[index] = 'C';
				break;
			case GUANINE:
				string[index] = 'G';
				break;
			case THYMINE:
				string[index] = 'T';
				break;
			}
		}

		string[length] = '\0';
	}

	return string;
}

boolean readStartsAreActivated(Graph * graph)
{
	return graph->nodeReads != NULL;
}

void activateReadStarts(Graph * graph)
{
	graph->nodeReads =
	    callocOrExit(2 * graph->nodeCount + 1, ShortReadMarker *);
	graph->nodeReadCounts =
	    callocOrExit(2 * graph->nodeCount + 1, IDnum);
}

void deactivateReadStarts(Graph * graph)
{
	free(graph->nodeReads);
	free(graph->nodeReadCounts);

	graph->nodeReads = NULL;
	graph->nodeReadCounts = NULL;
}

boolean findIDnumInArray(IDnum query, IDnum * array, IDnum arrayLength)
{
	IDnum leftIndex = 0;
	IDnum rightIndex = arrayLength;
	IDnum middleIndex;

	if (arrayLength == 0)
		return false;

	while (true) {
		middleIndex = leftIndex + (rightIndex - leftIndex) / 2;

		if (array[middleIndex] == query)
			return true;
		else if (leftIndex >= rightIndex)
			return false;
		else if (array[middleIndex] > query)
			rightIndex = middleIndex;
		else if (leftIndex == middleIndex)
			leftIndex++;
		else
			leftIndex = middleIndex;
	}
}

static inline int compareShortReadMarkers(const void *A, const void *B)
{
	IDnum a = ((ShortReadMarker *) A)->readID;
	IDnum b = ((ShortReadMarker *) B)->readID;

	if (a > b)
		return 1;
	if (a == b)
		return 0;
	return -1;
}

static inline int compareIDnums(const void *A, const void *B)
{
	IDnum a = *((IDnum *) A);
	IDnum b = *((IDnum *) B);

	if (a > b)
		return 1;
	if (a == b)
		return 0;
	return -1;
}

void incrementReadStartCount(Node * node, Graph * graph)
{
	graph->nodeReadCounts[node->ID + graph->nodeCount]++;
}

void createNodeReadStartArrays(Graph * graph)
{
	IDnum index;

	if (graph->nodeReads == NULL)
		return;

	for (index = 0; index <= 2 * (graph->nodeCount); index++) {
		if (graph->nodeReadCounts[index] != 0) {
			graph->nodeReads[index] =
			    mallocOrExit(graph->nodeReadCounts[index],
				   ShortReadMarker);
			graph->nodeReadCounts[index] = 0;
		} else {
			graph->nodeReads[index] = NULL;
		}
	}
}

void orderNodeReadStartArrays(Graph * graph)
{
	IDnum index;

	if (graph->nodeReads == NULL)
		return;

	for (index = 0; index <= 2 * (graph->nodeCount); index++)
		if (graph->nodeReadCounts[index] != 0)
			qsort(graph->nodeReads[index],
			      graph->nodeReadCounts[index],
			      sizeof(ShortReadMarker),
			      compareShortReadMarkers);
}

void addReadStart(Node * node, IDnum seqID, Coordinate position,
		  Graph * graph, Coordinate offset)
{
	IDnum nodeIndex = getNodeID(node) + graph->nodeCount;

	ShortReadMarker *array = graph->nodeReads[nodeIndex];
	IDnum arrayLength = graph->nodeReadCounts[nodeIndex];

	if (node->status)
		return;
	node->status = true;

	array[arrayLength].readID = seqID;
	array[arrayLength].position = position;
	array[arrayLength].offset = (ShortLength) offset;
	graph->nodeReadCounts[nodeIndex]++;
}

void blurLastShortReadMarker(Node * node, Graph * graph)
{
	IDnum nodeIndex = getNodeID(node) + nodeCount(graph);
	IDnum index = graph->nodeReadCounts[nodeIndex] - 1;
	ShortReadMarker *marker;

	if (index >= 0)
		marker = &(graph->nodeReads[nodeIndex][index]);
	else
		abort();

	setShortReadMarkerPosition(marker, -1);
}

ShortReadMarker *commonNodeReads(Node * nodeA, Node * nodeB, Graph * graph,
				 IDnum * length)
{
	IDnum targetID, targetLength, targetIndex, targetVal;
	IDnum sourceID, sourceLength, sourceIndex, sourceVal;
	IDnum mergeLength;
	ShortReadMarker *mergeArray, *targetArray, *sourceArray;

	if (graph->nodeReads == NULL) {
		*length = 0;
		return NULL;
	}

	if (nodeA == NULL || nodeB == NULL) {
		*length = 0;
		return NULL;
	}

	targetID = getNodeID(nodeA) + graph->nodeCount;
	targetArray = graph->nodeReads[targetID];
	targetLength = graph->nodeReadCounts[targetID];

	sourceID = getNodeID(nodeB) + graph->nodeCount;
	sourceArray = graph->nodeReads[sourceID];
	sourceLength = graph->nodeReadCounts[sourceID];

	if (sourceArray == NULL || targetArray == NULL) {
		*length = 0;
		return NULL;
	}

	mergeArray =
	    mallocOrExit(sourceLength +
		    targetLength, ShortReadMarker);

	mergeLength = 0;
	sourceIndex = 0;
	targetIndex = 0;
	sourceVal = sourceArray[0].readID;
	targetVal = targetArray[0].readID;

	while (sourceIndex < sourceLength && targetIndex < targetLength) {
		switch (compareIDnums(&sourceVal, &targetVal)) {
		case -1:
			mergeArray[mergeLength].readID = sourceVal;
			mergeArray[mergeLength].position = -1;
			mergeArray[mergeLength].offset = -1;
			mergeLength++;
			sourceIndex++;
			if (sourceIndex < sourceLength)
				sourceVal =
				    sourceArray[sourceIndex].readID;
			break;
		case 0:
			mergeArray[mergeLength].readID = sourceVal;
			mergeArray[mergeLength].position = -1;
			mergeArray[mergeLength].offset = -1;
			mergeLength++;
			sourceIndex++;
			if (sourceIndex < sourceLength)
				sourceVal =
				    sourceArray[sourceIndex].readID;
			targetIndex++;
			if (targetIndex < targetLength)
				targetVal =
				    targetArray[targetIndex].readID;
			break;
		case 1:
			mergeArray[mergeLength].readID = targetVal;
			mergeArray[mergeLength].position = -1;
			mergeArray[mergeLength].offset = -1;
			mergeLength++;
			targetIndex++;
			if (targetIndex < targetLength)
				targetVal =
				    targetArray[targetIndex].readID;
		}
	}

	while (sourceIndex < sourceLength) {
		mergeArray[mergeLength].readID =
		    sourceArray[sourceIndex].readID;
		mergeArray[mergeLength].position = -1;
		mergeArray[mergeLength].offset = -1;
		mergeLength++;
		sourceIndex++;
	}

	while (targetIndex < targetLength) {
		mergeArray[mergeLength].readID =
		    targetArray[targetIndex].readID;
		mergeArray[mergeLength].position = -1;
		mergeArray[mergeLength].offset = -1;
		mergeLength++;
		targetIndex++;
	}

	*length = mergeLength;
	return mergeArray;
}

ShortReadMarker *extractFrontOfNodeReads(Node * node,
					 Coordinate breakpoint,
					 Graph * graph, IDnum * length,
					 PassageMarker * sourceMarker,
					 Coordinate * lengths)
{
	IDnum sourceID;
	IDnum mergeLength, newLength, sourceLength;
	IDnum sourceIndex;
	ShortReadMarker *mergeArray, *sourceArray, *newArray;
	ShortReadMarker *mergePtr, *sourcePtr, *newPtr;
	Coordinate finish;
	Coordinate revBreakpoint;

	if (graph->nodeReads == NULL) {
		*length = 0;
		return NULL;
	}

	if (node == NULL) {
		*length = 0;
		return NULL;
	}

	if (breakpoint == 0) {
		return commonNodeReads(node,
				       getTwinNode(getNode
						   (getPreviousInSequence
						    (sourceMarker))),
				       graph, length);
	}

	sourceID = getNodeID(node) + graph->nodeCount;
	sourceArray = graph->nodeReads[sourceID];
	sourceLength = graph->nodeReadCounts[sourceID];

	if (sourceArray == NULL) {
		*length = 0;
		return NULL;
	}

	revBreakpoint = node->length - breakpoint;

	mergeLength = 0;
	newLength = 0;
	sourcePtr = sourceArray;
	for (sourceIndex = 0; sourceIndex < sourceLength; sourceIndex++) {
		if (sourcePtr->position == -1) {
			newLength++;
			mergeLength++;
		} else {
			finish =
			    sourcePtr->position - sourcePtr->offset +
			    lengths[sourcePtr->readID - 1];
			if (sourcePtr->position < revBreakpoint)
				newLength++;
			if (finish > revBreakpoint)
				mergeLength++;
		}
		sourcePtr++;
	}

	newArray = mallocOrExit(newLength, ShortReadMarker);
	mergeArray = mallocOrExit(mergeLength, ShortReadMarker);

	mergePtr = mergeArray;
	newPtr = newArray;
	sourcePtr = sourceArray;
	mergeLength = 0;
	newLength = 0;
	for (sourceIndex = 0; sourceIndex < sourceLength; sourceIndex++) {
		if (sourcePtr->position == -1) {
			mergePtr->readID = sourcePtr->readID;
			setShortReadMarkerPosition(mergePtr, -1);
			setShortReadMarkerOffset(mergePtr, -1);
			mergePtr++;
			mergeLength++;
			newPtr->readID = sourcePtr->readID;
			setShortReadMarkerPosition(newPtr, -1);
			setShortReadMarkerOffset(newPtr, -1);
			newPtr++;
			newLength++;
		} else {
			finish =
			    sourcePtr->position - sourcePtr->offset +
			    lengths[sourcePtr->readID - 1];
			if (sourcePtr->position < revBreakpoint) {
				newPtr->readID = sourcePtr->readID;
				setShortReadMarkerPosition(newPtr,
							   sourcePtr->
							   position);
				setShortReadMarkerOffset(newPtr,
							 sourcePtr->
							 offset);
				newPtr++;
				newLength++;

				// Saddle back reads:
				if (finish > revBreakpoint) {
					mergePtr->readID =
					    sourcePtr->readID;
					setShortReadMarkerPosition
					    (mergePtr, 0);
					setShortReadMarkerOffset(mergePtr,
								 sourcePtr->
								 offset +
								 revBreakpoint
								 -
								 sourcePtr->
								 position);
					mergePtr++;
				}
			} else if (finish > revBreakpoint) {
				mergePtr->readID = sourcePtr->readID;
				setShortReadMarkerPosition(mergePtr,
							   sourcePtr->
							   position - revBreakpoint);
				setShortReadMarkerOffset(mergePtr,
							 sourcePtr->
							 offset);
				mergePtr++;
				mergeLength++;
			}
		}

		sourcePtr++;
	}

	free(sourceArray);
	graph->nodeReads[sourceID] = newArray;
	graph->nodeReadCounts[sourceID] = newLength;

	*length = mergeLength;
	return mergeArray;
}

ShortReadMarker *extractBackOfNodeReads(Node * node, Coordinate breakpoint,
					Graph * graph, IDnum * length,
					PassageMarker * sourceMarker,
					Coordinate * lengths)
{
	IDnum sourceID;
	IDnum mergeLength, newLength, sourceLength;
	IDnum sourceIndex;
	ShortReadMarker *mergeArray, *sourceArray, *newArray;
	ShortReadMarker *mergePtr, *sourcePtr, *newPtr;
	Coordinate finish;

	if (graph->nodeReads == NULL) {
		*length = 0;
		return NULL;
	}

	if (node == NULL) {
		*length = 0;
		return NULL;
	}

	if (breakpoint == 0) {
		return
		    commonNodeReads(getNode
				    (getPreviousInSequence(sourceMarker)),
				    node, graph, length);
	}

	sourceID = getNodeID(node) + graph->nodeCount;
	sourceArray = graph->nodeReads[sourceID];
	sourceLength = graph->nodeReadCounts[sourceID];

	if (sourceArray == NULL) {
		*length = 0;
		return NULL;
	}

	mergeLength = 0;
	newLength = 0;
	sourcePtr = sourceArray;
	for (sourceIndex = 0; sourceIndex < sourceLength; sourceIndex++) {
		if (sourcePtr->position == -1) {
			mergeLength++;
			newLength++;
		} else {
			finish =
			    sourcePtr->position - sourcePtr->offset +
			    lengths[sourcePtr->readID - 1];
			if (sourcePtr->position < breakpoint)
				mergeLength++;
			if (finish > breakpoint)
				newLength++;
		}
		sourcePtr++;
	}

	newArray = mallocOrExit(newLength, ShortReadMarker);
	mergeArray = mallocOrExit(mergeLength, ShortReadMarker);

	mergePtr = mergeArray;
	newPtr = newArray;
	sourcePtr = sourceArray;
	for (sourceIndex = 0; sourceIndex < sourceLength; sourceIndex++) {
		if (sourcePtr->position == -1) {
			mergePtr->readID = sourcePtr->readID;
			setShortReadMarkerPosition(mergePtr, -1);
			setShortReadMarkerOffset(mergePtr, -1);
			mergePtr++;

			newPtr->readID = sourcePtr->readID;
			setShortReadMarkerPosition(newPtr, -1);
			setShortReadMarkerOffset(newPtr, -1);
			newPtr++;

			sourcePtr++;
			continue;
		} else {
			finish =
			    sourcePtr->position - sourcePtr->offset +
			    lengths[sourcePtr->readID - 1];

			if (sourcePtr->position < breakpoint) {
				mergePtr->readID = sourcePtr->readID;
				setShortReadMarkerPosition(mergePtr,
							   sourcePtr->
							   position);
				setShortReadMarkerOffset(mergePtr,
							 sourcePtr->
							 offset);
				mergePtr++;

				// Saddle back reads:
				if (finish > breakpoint) {
					newPtr->readID = sourcePtr->readID;
					setShortReadMarkerPosition(newPtr,
								   0);
					setShortReadMarkerOffset(newPtr,
								 sourcePtr->
								 offset +
								 breakpoint
								 -
								 sourcePtr->
								 position);
					newPtr++;
				}
			} else if (finish > breakpoint) {
				newPtr->readID = sourcePtr->readID;
				setShortReadMarkerPosition(newPtr,
							   sourcePtr->
							   position -
							   breakpoint);
				setShortReadMarkerOffset(newPtr,
							 sourcePtr->
							 offset);
				newPtr++;
			}
		}

		sourcePtr++;
	}

	free(sourceArray);
	graph->nodeReads[sourceID] = newArray;
	graph->nodeReadCounts[sourceID] = newLength;

	*length = mergeLength;
	return mergeArray;
}

void spreadReadIDs(ShortReadMarker * reads, IDnum readCount, Node * node,
		   Graph * graph)
{
	IDnum targetID, targetLength, targetIndex, targetVal;
	IDnum sourceLength, sourceIndex, sourceVal;
	IDnum mergeLength;
	ShortReadMarker *sourceArray, *targetArray, *mergeArray;
	ShortReadMarker *sourcePtr, *targetPtr, *mergePtr;
	Coordinate targetPosition;
	//ShortLength nodeLength = (ShortLength) getNodeLength(node);
	ShortLength targetOffset;

	if (graph->nodeReads == NULL || reads == NULL || node == NULL)
		return;

	targetID = getNodeID(node) + graph->nodeCount;
	targetArray = graph->nodeReads[targetID];
	targetLength = graph->nodeReadCounts[targetID];
	targetPtr = targetArray;

	sourceArray = reads;
	sourceLength = readCount;
	sourcePtr = sourceArray;

	if (targetArray == NULL) {
		mergeArray =
		    mallocOrExit(sourceLength, ShortReadMarker);
		mergePtr = mergeArray;

		sourceIndex = 0;
		while (sourceIndex < sourceLength) {
			mergePtr->readID = sourcePtr->readID;
			setShortReadMarkerPosition(mergePtr, -1);
			setShortReadMarkerOffset(mergePtr, -1);
			mergePtr++;
			sourcePtr++;
			sourceIndex++;
		}

		graph->nodeReads[targetID] = mergeArray;
		graph->nodeReadCounts[targetID] = sourceLength;
		return;
	}

	mergeArray =
	    mallocOrExit(sourceLength +
		    targetLength, ShortReadMarker);
	mergePtr = mergeArray;

	mergeLength = 0;
	sourceIndex = 0;
	targetIndex = 0;
	sourceVal = sourcePtr->readID;
	targetVal = targetPtr->readID;
	targetPosition = targetPtr->position;
	targetOffset = targetPtr->offset;

	while (sourceIndex < sourceLength && targetIndex < targetLength) {
		if (sourceVal < targetVal) {
			mergePtr->readID = sourceVal;
			setShortReadMarkerPosition(mergePtr, -1);
			setShortReadMarkerOffset(mergePtr, -1);
			sourceIndex++;
			sourcePtr++;
			if (sourceIndex < sourceLength)
				sourceVal = sourcePtr->readID;
		} else if (sourceVal == targetVal) {
			mergePtr->readID = sourceVal;
			setShortReadMarkerPosition(mergePtr, -1);
			setShortReadMarkerOffset(mergePtr, -1);
			sourceIndex++;
			sourcePtr++;
			if (sourceIndex < sourceLength)
				sourceVal = sourcePtr->readID;
			targetIndex++;
			targetPtr++;
			if (targetIndex < targetLength) {
				targetVal = targetPtr->readID;
				targetPosition = targetPtr->position;
				targetOffset = targetPtr->offset;
			}
		} else {
			mergePtr->readID = targetVal;
			setShortReadMarkerPosition(mergePtr,
						   targetPosition);
			setShortReadMarkerOffset(mergePtr, targetOffset);
			targetIndex++;
			targetPtr++;
			if (targetIndex < targetLength) {
				targetVal = targetPtr->readID;
				targetPosition = targetPtr->position;
				targetOffset = targetPtr->offset;
			}
		}

		mergeLength++;
		mergePtr++;
	}

	while (sourceIndex < sourceLength) {
		mergePtr->readID = sourcePtr->readID;
		setShortReadMarkerPosition(mergePtr, -1);
		setShortReadMarkerOffset(mergePtr, -1);
		mergeLength++;
		mergePtr++;
		sourceIndex++;
		sourcePtr++;
	}

	while (targetIndex < targetLength) {
		mergePtr->readID = targetPtr->readID;
		setShortReadMarkerPosition(mergePtr, targetPtr->position);
		setShortReadMarkerOffset(mergePtr, targetPtr->offset);
		mergeLength++;
		mergePtr++;
		targetIndex++;
		targetPtr++;
	}

	free(targetArray);
	graph->nodeReads[targetID] = mergeArray;
	graph->nodeReadCounts[targetID] = mergeLength;
}

static inline Coordinate min(Coordinate A, Coordinate B)
{
	return A < B ? A : B;
}

static inline ShortLength min_short(ShortLength A, ShortLength B)
{
	return A < B ? A : B;
}

void injectShortReads(ShortReadMarker * sourceArray, IDnum sourceLength,
		      Node * target, Graph * graph)
{
	IDnum targetID = getNodeID(target) + graph->nodeCount;
	ShortReadMarker *targetArray = graph->nodeReads[targetID];
	IDnum targetLength = graph->nodeReadCounts[targetID];
	ShortReadMarker *targetPtr = targetArray;
	ShortReadMarker *sourcePtr = sourceArray;
	ShortReadMarker *mergeArray, *mergePtr;
	IDnum mergeLength;
	Coordinate targetPosition, sourcePosition;
	ShortLength targetOffset, sourceOffset;
	IDnum targetIndex, targetVal, sourceIndex, sourceVal;

	if (sourceLength == 0) {
		free(sourceArray);
		return;
	}

	if (targetLength == 0) {
		free(targetArray);
		graph->nodeReads[targetID] = sourceArray;
		graph->nodeReadCounts[targetID] = sourceLength;
		return;
	}

	mergeArray =
	    mallocOrExit(sourceLength +
		    targetLength, ShortReadMarker);
	mergePtr = mergeArray;

	mergeLength = 0;
	sourceIndex = 0;
	targetIndex = 0;
	targetVal = targetPtr->readID;
	targetPosition = targetPtr->position;
	targetOffset = targetPtr->offset;
	sourceVal = sourcePtr->readID;
	sourcePosition = sourcePtr->position;
	sourceOffset = sourcePtr->offset;

	while (sourceIndex < sourceLength && targetIndex < targetLength) {
		if (sourceVal < targetVal) {
			mergePtr->readID = sourceVal;
			setShortReadMarkerPosition(mergePtr,
						   sourcePosition);
			setShortReadMarkerOffset(mergePtr, sourceOffset);
			sourceIndex++;
			if (sourceIndex < sourceLength) {
				sourcePtr++;
				sourceVal = sourcePtr->readID;
				sourcePosition = sourcePtr->position;
				sourceOffset = sourcePtr->offset;
			}
		} else if (sourceVal == targetVal) {
			mergePtr->readID = sourceVal;
			if (sourcePosition == -1 && targetPosition == -1) {
				setShortReadMarkerPosition(mergePtr, -1);
				setShortReadMarkerOffset(mergePtr, -1);
			} else if (sourcePosition == -1) {
				setShortReadMarkerPosition(mergePtr,
							   targetPosition);
				setShortReadMarkerOffset(mergePtr,
							 targetOffset);
			} else if (targetPosition == -1) {
				setShortReadMarkerPosition(mergePtr,
							   sourcePosition);
				setShortReadMarkerOffset(mergePtr,
							 sourceOffset);
			} else {
				setShortReadMarkerPosition(mergePtr,
							   min
							   (sourcePosition,
							    targetPosition));
				setShortReadMarkerOffset(mergePtr,
							 min_short
							 (sourceOffset,
							  targetOffset));
			}
			sourceIndex++;
			if (sourceIndex < sourceLength) {
				sourcePtr++;
				sourceVal = sourcePtr->readID;
				sourcePosition = sourcePtr->position;
				sourceOffset = sourcePtr->offset;
			}
			targetIndex++;
			if (targetIndex < targetLength) {
				targetPtr++;
				targetVal = targetPtr->readID;
				targetPosition = targetPtr->position;
				targetOffset = targetPtr->offset;
			}
		} else {
			mergePtr->readID = targetVal;
			setShortReadMarkerPosition(mergePtr,
						   targetPosition);
			setShortReadMarkerOffset(mergePtr, targetOffset);
			targetIndex++;
			if (targetIndex < targetLength) {
				targetPtr++;
				targetVal = targetPtr->readID;
				targetPosition = targetPtr->position;
				targetOffset = targetPtr->offset;
			}
		}

		mergeLength++;
		mergePtr++;
	}

	while (sourceIndex < sourceLength) {
		mergePtr->readID = sourcePtr->readID;
		setShortReadMarkerPosition(mergePtr, sourcePtr->position);
		setShortReadMarkerOffset(mergePtr, sourcePtr->offset);
		mergeLength++;
		mergePtr++;
		sourceIndex++;
		sourcePtr++;
	}

	while (targetIndex < targetLength) {
		mergePtr->readID = targetPtr->readID;
		setShortReadMarkerPosition(mergePtr, targetPtr->position);
		setShortReadMarkerOffset(mergePtr, targetPtr->offset);
		mergeLength++;
		mergePtr++;
		targetIndex++;
		targetPtr++;
	}

	free(targetArray);
	graph->nodeReads[targetID] = mergeArray;
	graph->nodeReadCounts[targetID] = mergeLength;

	free(sourceArray);
}

void mergeNodeReads(Node * target, Node * source, Graph * graph)
{
	IDnum sourceID, sourceLength;
	ShortReadMarker *sourceArray;

	if (graph->nodeReads == NULL || source == NULL || target == NULL)
		return;

	sourceID = getNodeID(source) + graph->nodeCount;
	sourceArray = graph->nodeReads[sourceID];
	sourceLength = graph->nodeReadCounts[sourceID];

	if (sourceArray == NULL)
		return;

	graph->nodeReads[sourceID] = NULL;
	graph->nodeReadCounts[sourceID] = 0;

	injectShortReads(sourceArray, sourceLength, target, graph);
}

void foldSymmetricalNodeReads(Node * node, Graph * graph)
{
	IDnum targetID, targetLength, targetIndex;
	IDnum sourceID, sourceLength, sourceIndex;
	IDnum targetVal = 0;
	IDnum sourceVal = 0;
	IDnum mergeLength;
	ShortReadMarker *sourceArray, *targetArray, *mergeArray,
	    *mergeArray2;
	ShortReadMarker *sourcePtr, *targetPtr, *mergePtr, *mergePtr2;

	if (graph->nodeReads == NULL || node == NULL)
		return;

	sourceID = getNodeID(node) + graph->nodeCount;
	sourceArray = graph->nodeReads[sourceID];
	sourceLength = graph->nodeReadCounts[sourceID];
	sourcePtr = sourceArray;

	targetID = -getNodeID(node) + graph->nodeCount;
	targetArray = graph->nodeReads[targetID];
	targetLength = graph->nodeReadCounts[targetID];
	targetPtr = targetArray;

	if (sourceArray == NULL && targetArray == NULL)
		return;

	mergeArray =
	    mallocOrExit(sourceLength +
		    targetLength, ShortReadMarker);
	mergeArray2 =
	    mallocOrExit(sourceLength +
		    targetLength, ShortReadMarker);
	mergePtr = mergeArray;
	mergePtr2 = mergeArray2;

	mergeLength = 0;
	sourceIndex = 0;
	targetIndex = 0;
	if (targetIndex < targetLength)
		targetVal = targetPtr->readID;
	if (sourceIndex < sourceLength)
		sourceVal = sourcePtr->readID;

	while (sourceIndex < sourceLength && targetIndex < targetLength) {
		if (sourceVal < targetVal) {
			mergePtr->readID = sourceVal;
			setShortReadMarkerPosition(mergePtr, -1);
			setShortReadMarkerOffset(mergePtr, -1);
			mergePtr2->readID = sourceVal;
			setShortReadMarkerPosition(mergePtr2, -1);
			setShortReadMarkerOffset(mergePtr2, -1);
			sourceIndex++;
			sourcePtr++;
			if (sourceIndex < sourceLength)
				sourceVal = sourcePtr->readID;
		} else if (sourceVal == targetVal) {
			mergePtr->readID = sourceVal;
			setShortReadMarkerPosition(mergePtr, -1);
			setShortReadMarkerOffset(mergePtr, -1);
			mergePtr2->readID = sourceVal;
			setShortReadMarkerPosition(mergePtr2, -1);
			setShortReadMarkerOffset(mergePtr2, -1);
			sourceIndex++;
			sourcePtr++;
			if (sourceIndex < sourceLength)
				sourceVal = sourcePtr->readID;
			targetIndex++;
			targetPtr++;
			if (targetIndex < targetLength)
				targetVal = targetPtr->readID;
		} else {
			mergePtr->readID = targetVal;
			setShortReadMarkerPosition(mergePtr, -1);
			setShortReadMarkerOffset(mergePtr, -1);
			mergePtr2->readID = targetVal;
			setShortReadMarkerPosition(mergePtr2, -1);
			setShortReadMarkerOffset(mergePtr2, -1);
			targetIndex++;
			targetPtr++;
			if (targetIndex < targetLength)
				targetVal = targetPtr->readID;
		}

		mergeLength++;
		mergePtr++;
		mergePtr2++;
	}

	while (sourceIndex < sourceLength) {
		mergePtr->readID = sourcePtr->readID;
		setShortReadMarkerPosition(mergePtr, -1);
		setShortReadMarkerOffset(mergePtr, -1);
		mergePtr2->readID = sourcePtr->readID;
		setShortReadMarkerPosition(mergePtr2, -1);
		setShortReadMarkerOffset(mergePtr2, -1);
		mergeLength++;
		mergePtr++;
		mergePtr2++;
		sourceIndex++;
		sourcePtr++;
	}

	while (targetIndex < targetLength) {
		mergePtr->readID = targetPtr->readID;
		setShortReadMarkerPosition(mergePtr, -1);
		setShortReadMarkerOffset(mergePtr, -1);
		mergePtr2->readID = targetPtr->readID;
		setShortReadMarkerPosition(mergePtr2, -1);
		setShortReadMarkerOffset(mergePtr2, -1);
		mergeLength++;
		mergePtr++;
		mergePtr2++;
		targetIndex++;
		targetPtr++;
	}

	free(targetArray);
	graph->nodeReads[targetID] = mergeArray;
	graph->nodeReadCounts[targetID] = mergeLength;

	free(sourceArray);
	graph->nodeReads[sourceID] = mergeArray2;
	graph->nodeReadCounts[sourceID] = mergeLength;
}

void shareReadStarts(Node * target, Node * source, Graph * graph)
{
	ShortReadMarker *sourceArray;
	IDnum sourceLength, sourceID;

	if (graph->nodeReads == NULL)
		return;

	if (target == NULL || source == NULL)
		return;

	sourceID = source->ID + graph->nodeCount;
	sourceArray = graph->nodeReads[sourceID];
	sourceLength = graph->nodeReadCounts[sourceID];

	if (sourceArray == NULL)
		return;

	spreadReadIDs(sourceArray, sourceLength, target, graph);
}

ShortReadMarker **getNodeToReadMappings(Graph * graph)
{
	return graph->nodeReads;
}

IDnum getShortReadMarkerID(ShortReadMarker * marker)
{
	return marker->readID;
}

inline ShortLength getShortReadMarkerOffset(ShortReadMarker * marker)
{
	return marker->offset;
}

inline void setShortReadMarkerOffset(ShortReadMarker * marker,
				     ShortLength offset)
{
	marker->offset = offset;
}

IDnum *getNodeReadCounts(Graph * graph)
{
	return graph->nodeReadCounts;
}

int getWordLength(Graph * graph)
{
	return graph->wordLength;
}

void displayArcMemory()
{
	printf("ARC MEMORY %lld allocated %lld free\n",
	       (long long) RecycleBin_memory_usage(arcMemory),
	       (long long) recycleBinFreeSpace(arcMemory));
}

void displayNodeMemory()
{
	printf("NODE MEMORY %lld allocated %lld free\n",
	       (long long) RecycleBin_memory_usage(nodeMemory),
	       (long long) recycleBinFreeSpace(nodeMemory));
}

ShortReadMarker *getNodeReads(Node * node, Graph * graph)
{
	IDnum id = node->ID + graph->nodeCount;
	return graph->nodeReads[id];
}

IDnum getNodeReadCount(Node * node, Graph * graph)
{
	IDnum id = node->ID + graph->nodeCount;
	return graph->nodeReadCounts[id];
}

inline Coordinate getShortReadMarkerPosition(ShortReadMarker * marker)
{
	return marker->position;
}

inline void setShortReadMarkerPosition(ShortReadMarker * marker,
				       Coordinate position)
{
	if (position < -100)
		return;

	marker->position = position;
}

ShortReadMarker *getShortReadMarkerAtIndex(ShortReadMarker * array,
					   IDnum index)
{
	return &(array[index]);
}

void destroyGraph(Graph * graph)
{
	IDnum index;
	Node *node;
	for (index = 1; index <= graph->nodeCount; index++) {
		node = getNodeInGraph(graph, index);
		if (node != NULL)
			destroyNode(node, graph);
	}

	if (graph->gapMarkers)
		deactivateGapMarkers(graph);

	free(graph->nodes);
	destroyRecycleBin(nodeMemory);
	destroyRecycleBin(arcMemory);
	destroyAllPassageMarkers();
	free(graph->arcLookupTable);
	free(graph->nodeReads);
	free(graph->nodeReadCounts);
	free(graph);
}

void checkNodeReads(IDnum index, Graph * graph)
{
	IDnum ref = index + graph->nodeCount;
	IDnum arrayLength = graph->nodeReadCounts[ref];
	ShortReadMarker *array = graph->nodeReads[ref];
	IDnum i;

	//return;

	if (arrayLength > graph->sequenceCount)
		abort();

	//if (arrayLength > 10000)
	//      printf("Array length %d %d\n", arrayLength, index);

	for (i = 1; i < arrayLength; i++) {
		if (array[i].readID <= array[i - 1].readID)
			abort();
		if (array[i].position >= 0 && array[i].offset < 0)
			abort();
		if (array[i - 1].position >= 0 && array[i - 1].offset < 0)
			abort();
	}
}

void setInsertLengths(Graph * graph, Category cat, Coordinate insertLength,
		      Coordinate insertLength_std_dev)
{
	graph->insertLengths[cat] = insertLength;
	graph->insertLengths_var[cat] =
	    insertLength_std_dev * insertLength_std_dev;
}

Coordinate getInsertLength(Graph * graph, Category cat)
{
	return graph->insertLengths[cat / 2];
}

double getInsertLength_var(Graph * graph, Category cat)
{
	return graph->insertLengths_var[cat / 2];
}

void activateGapMarkers(Graph * graph)
{
	graph->gapMarkers =
	    callocOrExit(graph->nodeCount + 1, GapMarker *);
	gapMarkerMemory = newRecycleBin(sizeof(GapMarker), GAPBLOCKSIZE);
}

void deactivateGapMarkers(Graph * graph)
{
	free(graph->gapMarkers);
	graph->gapMarkers = NULL;
	destroyRecycleBin(gapMarkerMemory);
	gapMarkerMemory = NULL;
}

static GapMarker *allocateGapMarker()
{
	return (GapMarker *) allocatePointer(gapMarkerMemory);
}

void appendGap(Node * node, Coordinate length, Graph * graph)
{
	IDnum nodeID = getNodeID(node);
	GapMarker *marker = allocateGapMarker();
	GapMarker *tmp;

	marker->length = length;

	if (nodeID > 0) {
		marker->position = node->length;
		marker->next = graph->gapMarkers[nodeID];
		graph->gapMarkers[nodeID] = marker;
	} else {
		for (tmp = graph->gapMarkers[-nodeID]; tmp != NULL;
		     tmp = tmp->next)
			tmp->position += length;

		marker->position = 0;
		marker->next = graph->gapMarkers[-nodeID];
		graph->gapMarkers[-nodeID] = marker;
	}

	addBufferToDescriptor(node, length);
}

void appendNodeGaps(Node * destination, Node * source, Graph * graph)
{
	IDnum destinationID = getNodeID(destination);
	IDnum sourceID = getNodeID(source);
	GapMarker *marker;

	if (graph->gapMarkers == NULL)
		return;

	if (destinationID > 0 && sourceID > 0) {
		for (marker = graph->gapMarkers[sourceID]; marker != NULL;
		     marker = marker->next)
			marker->position += destination->length;
	} else if (destinationID > 0 && sourceID < 0) {
		sourceID = -sourceID;
		for (marker = graph->gapMarkers[sourceID]; marker != NULL;
		     marker = marker->next)
			marker->position =
			    source->length + destination->length -
			    marker->position - marker->length;
	} else if (destinationID < 0 && sourceID > 0) {
		destinationID = -destinationID;
		for (marker = graph->gapMarkers[destinationID];
		     marker != NULL; marker = marker->next)
			marker->position += source->length;

		for (marker = graph->gapMarkers[sourceID]; marker != NULL;
		     marker = marker->next)
			marker->position =
			    source->length - marker->position -
			    marker->length;
	} else {
		destinationID = -destinationID;
		sourceID = -sourceID;
		for (marker = graph->gapMarkers[destinationID];
		     marker != NULL; marker = marker->next)
			marker->position += source->length;
	}

	if (graph->gapMarkers[destinationID] == NULL)
		graph->gapMarkers[destinationID] =
		    graph->gapMarkers[sourceID];
	else {
		marker = graph->gapMarkers[destinationID];
		while (marker->next != NULL)
			marker = marker->next;
		marker->next = graph->gapMarkers[sourceID];
	}

	graph->gapMarkers[sourceID] = NULL;
}

GapMarker *getGap(Node * node, Graph * graph)
{
	IDnum nodeID = getNodeID(node);

	if (graph->gapMarkers == NULL)
		return NULL;

	if (nodeID < 0)
		nodeID = -nodeID;

	return graph->gapMarkers[nodeID];
}

GapMarker *getNextGap(GapMarker * marker)
{
	return marker->next;
}

Coordinate getGapStart(GapMarker * marker)
{
	return marker->position;
}

Coordinate getGapFinish(GapMarker * marker)
{
	return marker->position + marker->length;
}

void reallocateNodeDescriptor(Node * node, Coordinate length) {
	Coordinate arrayLength, index, shift;
	Node * twin = node->twinNode;
	Descriptor * array;
	Nucleotide nucleotide;

	if (length < node->length)
		exitErrorf(EXIT_FAILURE, true, "Sum of node lengths smaller than first!");

	shift = length - node->length;

	arrayLength = length / 4;
	if (length % 4)
		arrayLength++;

	node->descriptor = reallocOrExit(node->descriptor, arrayLength, Descriptor);

	array = callocOrExit(arrayLength, Descriptor);
	for (index = node->length - 1; index >= 0; index--) {
		nucleotide = getNucleotideInDescriptor(twin->descriptor, index);
		writeNucleotideInDescriptor(nucleotide, array, index + shift);
	}
	
	free(twin->descriptor);
	twin->descriptor = array;
}
