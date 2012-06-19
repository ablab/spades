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
#include <ctype.h>

#include "globals.h"
#include "preGraph.h"
#include "recycleBin.h"
#include "tightString.h"
#include "run.h"
#include "utility.h"

#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

struct preArc_st {
	PreArc *nextLeft;
	PreArc *nextRight;
	IDnum multiplicity;
	IDnum preNodeIDLeft;
	IDnum preNodeIDRight;
};

struct preNode_st {
	PreArc *preArcLeft;
	PreArc *preArcRight;
	Descriptor *descriptor;
	Coordinate length;
};

struct preGraph_st {
	PreNode *preNodes;
	IDnum sequenceCount;
	IDnum preNodeCount;
	int wordLength;
	boolean double_strand;
};

static RecycleBin *preArcMemory = NULL;

#define BLOCKSIZE 10000

PreArc *allocatePreArc_pg()
{
	if (preArcMemory == NULL)
		preArcMemory = newRecycleBin(sizeof(PreArc), BLOCKSIZE);

	return allocatePointer(preArcMemory);
}

void deallocatePreArc_pg(PreArc * preArc)
{
	deallocatePointer(preArcMemory, preArc);
}

// Returns the length of the preNode's descriptor list
Coordinate getPreNodeLength_pg(IDnum preNodeID, PreGraph * preGraph)
{
	IDnum ID = preNodeID;

	if (ID < 0)
		ID = -ID;

	return (preGraph->preNodes[ID]).length;
}

// Returns the number of preNodes in the preGraph
IDnum preNodeCount_pg(PreGraph * preGraph)
{
	return preGraph->preNodeCount;
}

// returns the number of sequences used to buid the preGraph
IDnum sequenceCount_pg(PreGraph * preGraph)
{
	return preGraph->sequenceCount;
}

PreArc *getPreArcBetweenPreNodes_pg(IDnum originPreNodeID,
				    IDnum destinationPreNodeID,
				    PreGraph * preGraph)
{
	PreArc *preArc;

	if (originPreNodeID == 0 || destinationPreNodeID == 0) {
		return NULL;
	}

	for (preArc = getPreArc_pg(originPreNodeID, preGraph);
	     preArc != NULL;
	     preArc = getNextPreArc_pg(preArc, originPreNodeID)) {
		if (getDestination_pg(preArc, originPreNodeID) ==
		    destinationPreNodeID) {
			return preArc;
		}
	}

	return NULL;
}

static void addPreArcToPreNode_pg(PreArc * preArc, IDnum preNodeID,
				  PreGraph * preGraph)
{
	IDnum ID = preNodeID;
	PreNode *preNode;
	PreArc **preArcPtr;

	if (ID < 0)
		ID = -ID;

	preNode = &(preGraph->preNodes[ID]);

	if (preNodeID > 0)
		preArcPtr = &(preNode->preArcRight);
	else
		preArcPtr = &(preNode->preArcLeft);

	if (preNodeID == preArc->preNodeIDLeft) {
		preArc->nextLeft = *preArcPtr;
		*preArcPtr = preArc;
	}

	if (preNodeID == preArc->preNodeIDRight) {
		preArc->nextRight = *preArcPtr;
		*preArcPtr = preArc;
	}
}

// Creates an preArc from preNode origin to preNode destination.
// If this preArc already exists, increments its multiplicity by 1.
PreArc *createPreArc_pg(IDnum originPreNodeID, IDnum destinationPreNodeID,
			PreGraph * preGraph)
{
	PreArc *preArc;


	if (originPreNodeID == 0 || destinationPreNodeID == 0)
		return NULL;

	preArc =
	    getPreArcBetweenPreNodes_pg(originPreNodeID,
					destinationPreNodeID, preGraph);

	if (preArc != NULL) {
		preArc->multiplicity++;
		return preArc;
	}
	// If not found
	preArc = allocatePreArc_pg();
	preArc->preNodeIDLeft = originPreNodeID;
	preArc->preNodeIDRight = -destinationPreNodeID;
	preArc->multiplicity = 1;

	addPreArcToPreNode_pg(preArc, originPreNodeID, preGraph);

	// Hairpin case
	if (destinationPreNodeID == -originPreNodeID) {
		preArc->multiplicity++;
		return preArc;
	}

	addPreArcToPreNode_pg(preArc, -destinationPreNodeID, preGraph);

	return preArc;
}

void createAnalogousPreArc_pg(IDnum originPreNodeID,
			      IDnum destinationPreNodeID,
			      PreArc * refPreArc, PreGraph * preGraph)
{
	PreArc *preArc;

	if (originPreNodeID == 0 || destinationPreNodeID == 0)
		return;

	preArc =
	    getPreArcBetweenPreNodes_pg(originPreNodeID,
					destinationPreNodeID, preGraph);

	if (preArc != NULL) {
		preArc->multiplicity += refPreArc->multiplicity;
		return;
	}
	// If not found
	preArc = allocatePreArc_pg();
	preArc->preNodeIDLeft = originPreNodeID;
	preArc->preNodeIDRight = -destinationPreNodeID;
	preArc->multiplicity = refPreArc->multiplicity;

	addPreArcToPreNode_pg(preArc, originPreNodeID, preGraph);

	// Hairpin case
	if (destinationPreNodeID == -originPreNodeID) {
		preArc->multiplicity++;
		return;
	}

	addPreArcToPreNode_pg(preArc, -destinationPreNodeID, preGraph);
}

void changeMultiplicity_pg(PreArc * preArc, IDnum variation)
{
	if (preArc == NULL)
		return;
	preArc->multiplicity += variation;
}

static void setNextPreArc_pg(PreArc * preArc, IDnum preNodeID,
			     PreArc * nextPreArc)
{
	if (preNodeID == preArc->preNodeIDLeft)
		preArc->nextLeft = nextPreArc;
	if (preNodeID == preArc->preNodeIDRight)
		preArc->nextRight = nextPreArc;
}

void removePreArcFromList_pg(PreArc * preArc, IDnum preNodeID,
			     PreGraph * preGraph)
{
	IDnum ID = preNodeID;
	PreNode *preNode;
	PreArc **preArcPtr;
	PreArc *tempPreArc;

	if (ID < 0)
		ID = -ID;

	preNode = &(preGraph->preNodes[ID]);

	if (preNodeID > 0)
		preArcPtr = &(preNode->preArcRight);
	else
		preArcPtr = &(preNode->preArcLeft);

	if (*preArcPtr == preArc) {
		*preArcPtr = getNextPreArc_pg(preArc, preNodeID);
		return;
	}

	for (tempPreArc = *preArcPtr; tempPreArc != NULL;
	     tempPreArc = getNextPreArc_pg(tempPreArc, preNodeID))
		if (getNextPreArc_pg(tempPreArc, preNodeID) == preArc)
			setNextPreArc_pg(tempPreArc, preNodeID,
					 getNextPreArc_pg(preArc,
							  preNodeID));
}

void destroyPreArc_pg(PreArc * preArc, PreGraph * preGraph)
{
	IDnum leftID, rightID;

	if (preArc == NULL)
		return;

	leftID = preArc->preNodeIDLeft;
	rightID = preArc->preNodeIDRight;

	// Removing preArc from list
	removePreArcFromList_pg(preArc, leftID, preGraph);

	// Removing preArc's twin from list
	if (rightID != leftID)
		removePreArcFromList_pg(preArc, rightID, preGraph);

	deallocatePreArc_pg(preArc);
}

void destroyPreNode_pg(IDnum preNodeID, PreGraph * preGraph)
{
	PreNode *preNode;
	IDnum ID = preNodeID;

	//printf("Destroying %ld\n and twin %ld\n", getPreNodeID(preNode), getPreNodeID(twin));

	if (ID < 0)
		ID = -ID;

	preNode = &(preGraph->preNodes[ID]);

	// PreNode preArcs:
	while (preNode->preArcLeft != NULL)
		destroyPreArc_pg(preNode->preArcLeft, preGraph);
	while (preNode->preArcRight != NULL)
		destroyPreArc_pg(preNode->preArcRight, preGraph);

	// Descriptors
	free(preNode->descriptor);

	// Flag as destroyed
	preNode->descriptor = NULL;
}

void destroyPreGraph_pg(PreGraph * preGraph)
{
	IDnum index;
	PreNode *preNode = &(preGraph->preNodes[1]);

	// Descriptors
	for (index = 1; index <= preGraph->preNodeCount; index++) {
		free(preNode->descriptor);
		preNode++;
	}

	// Arcs
	destroyRecycleBin(preArcMemory);

	// Nodes
	free(preGraph->preNodes);

	// Graph
	free(preGraph);

}

static Nucleotide getNucleotideInDescriptor_pg(Descriptor * descriptor,
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

PreNode *getPreNodeInPreGraph_pg(PreGraph * preGraph, IDnum preNodeID)
{
	PreNode *preNode;
	if (preNodeID <= 0)
		abort();
	else {
		preNode = &(preGraph->preNodes[preNodeID]);
		if (preNode->descriptor != NULL)
			return preNode;
		else
			return NULL;
	}
	return NULL;
}

PreArc *getPreArc_pg(IDnum preNodeID, PreGraph * preGraph)
{
	IDnum ID = preNodeID;
	PreNode *preNode;

	if (ID < 0)
		ID = -ID;

	preNode = &(preGraph->preNodes[ID]);

	if (preNodeID > 0)
		return preNode->preArcRight;
	else
		return preNode->preArcLeft;
}

PreArc *getNextPreArc_pg(PreArc * preArc, IDnum preNodeID)
{
	if (preNodeID == preArc->preNodeIDLeft)
		return preArc->nextLeft;
	else
		return preArc->nextRight;
}

IDnum getMultiplicity_pg(PreArc * preArc)
{
	if (preArc == NULL)
		return 0;

	return preArc->multiplicity;
}

IDnum getOtherEnd_pg(PreArc * preArc, IDnum preNodeID)
{
	if (preNodeID == preArc->preNodeIDLeft)
		return preArc->preNodeIDRight;
	else
		return preArc->preNodeIDLeft;
}

IDnum getDestination_pg(PreArc * preArc, IDnum preNodeID)
{
	if (preArc == NULL)
		return 0;

	if (preNodeID == preArc->preNodeIDLeft)
		return -preArc->preNodeIDRight;
	else
		return -preArc->preNodeIDLeft;
}

static void writeNucleotideInDescriptor_pg(Nucleotide nucleotide,
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

static inline Descriptor *mergeDescriptors_pg(Descriptor * descr,
					      Coordinate destinationLength,
					      Descriptor * copy,
					      Coordinate sourceLength,
					      int wordLength)
{
	Descriptor *readPtr, *writePtr;
	Descriptor readCopy = 0;
	int readOffset, writeOffset;
	size_t arrayLength;
	Coordinate newLength =
	    destinationLength + sourceLength + wordLength - 1;
	Descriptor *new;
	Coordinate index;

	// Specify new array
	arrayLength = newLength / 4;
	if (newLength % 4)
		arrayLength++;
	new = callocOrExit(arrayLength, Descriptor);
	for (index = 0; index < arrayLength; index++)
		new[index] = 0;

	// Copying first descriptor
	readPtr = descr;
	writePtr = new;
	writeOffset = 0;
	for (index = 0; index < destinationLength + wordLength - 1;
	     index++) {
		(*writePtr) >>= 2;
		if (writeOffset == 0)
			readCopy = *readPtr;
		(*writePtr) += (readCopy & 3) << 6;

		/*switch ((readCopy & 3)) {
		   case ADENINE:
		   printf("A%ld", index);
		   break;
		   case CYTOSINE:
		   printf("C%ld", index);
		   break;
		   case GUANINE:
		   printf("G%ld", index);
		   break;
		   case THYMINE:
		   printf("T%ld", index);
		   break;
		   } */
		readCopy >>= 2;

		writeOffset++;
		if (writeOffset == 4) {
			writePtr++;
			readPtr++;
			writeOffset = 0;
		}
	}

	//puts("");

	// Skipping initial k-1 letters in second descriptor
	readPtr = &(copy[(wordLength - 1) / 4]);
	readCopy = *readPtr;
	readOffset = (wordLength - 1) % 4;
	readCopy >>= (readOffset * 2);

	// Going on copying second descriptor
	for (index = 0; index < sourceLength; index++) {
		(*writePtr) >>= 2;
		if (readOffset == 0)
			readCopy = *readPtr;
		(*writePtr) += (readCopy & 3) << 6;
		/*switch ((readCopy & 3)) {
		   case ADENINE:
		   printf("A%ld", index);
		   break;
		   case CYTOSINE:
		   printf("C%ld", index);
		   break;
		   case GUANINE:
		   printf("G%ld", index);
		   break;
		   case THYMINE:
		   printf("T%ld", index);
		   break;
		   default:
		   printf("?%ld;", index);
		   } */
		readCopy >>= 2;

		writeOffset++;
		if (writeOffset == 4) {
			writePtr++;
			writeOffset = 0;
		}

		readOffset++;
		if (readOffset == 4) {
			readPtr++;
			readOffset = 0;
		}
	}

	//puts("");

	if (writeOffset != 0) {
		while (writeOffset != 4) {
			(*writePtr) >>= 2;
			writeOffset++;
		}
	}

	return new;
}

static inline Descriptor *mergeDescriptorsH2H_pg(Descriptor * descr,
						 Coordinate
						 destinationLength,
						 Descriptor * copy,
						 Coordinate sourceLength,
						 int wordLength)
{
	Descriptor *readPtr, *writePtr;
	Descriptor readCopy;
	int readOffset, writeOffset;
	size_t arrayLength;
	Coordinate newLength =
	    destinationLength + sourceLength + wordLength - 1;
	Descriptor *new;
	Coordinate index;

	// Specify new array
	arrayLength = newLength / 4;
	if (newLength % 4)
		arrayLength++;
	new = callocOrExit(arrayLength, Descriptor);
	for (index = 0; index < arrayLength; index++)
		new[index] = 0;

	// Copying first descriptor (including final (k-1)-mer)
	readPtr = descr;
	readCopy = *readPtr;
	writePtr = new;
	writeOffset = 0;
	readOffset = 0;
	for (index = 0; index < destinationLength + wordLength - 1;
	     index++) {
		(*writePtr) >>= 2;
		if (writeOffset == 0)
			readCopy = *readPtr;
		(*writePtr) += (readCopy & 3) << 6;
		/*switch ((readCopy & 3)) {
		   case ADENINE:
		   printf("A(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case CYTOSINE:
		   printf("C(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case GUANINE:
		   printf("G(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case THYMINE:
		   printf("T(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   default:
		   printf("?(%ld %i %i);", index, writeOffset, readOffset);
		   } */
		readCopy >>= 2;

		writeOffset++;
		if (writeOffset == 4) {
			writePtr++;
			readPtr++;
			writeOffset = 0;
		}
	}

	//puts("");

	// Going to end of second descriptor 
	readPtr = &(copy[(sourceLength - 1) / 4]);
	readCopy = *readPtr;
	readOffset = (sourceLength - 1) % 4;
	readCopy <<= ((3 - readOffset) * 2);

	//printf("Read copy %x\n", readCopy);

	// Going on copying reverse complement of second descriptor
	for (index = 0; index < sourceLength; index++) {
		(*writePtr) >>= 2;
		if (readOffset == 3)
			readCopy = *readPtr;
#ifndef COLOR
		(*writePtr) += 192 - (readCopy & 192);
#else
		(*writePtr) += (readCopy & 192);
#endif
		/*switch (3 - ((readCopy & 192) >> 6)) {
		   case ADENINE:
		   printf("A(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case CYTOSINE:
		   printf("C(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case GUANINE:
		   printf("G(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case THYMINE:
		   printf("T(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   default:
		   printf("?(%ld %i %i);", index, writeOffset, readOffset);
		   } */
		readCopy <<= 2;

		writeOffset++;
		if (writeOffset == 4) {
			writePtr++;
			writeOffset = 0;
		}

		readOffset--;
		if (readOffset == -1) {
			readPtr--;
			readOffset = 3;
		}
	}

	//puts("");

	if (writeOffset != 0) {
		while (writeOffset != 4) {
			(*writePtr) >>= 2;
			writeOffset++;
		}
	}

	return new;
}

static inline Descriptor *mergeDescriptorsF2F_pg(Descriptor * descr,
						 Coordinate
						 destinationLength,
						 Descriptor * copy,
						 Coordinate sourceLength,
						 int wordLength)
{
	Descriptor *readPtr, *writePtr;
	Descriptor readCopy;
	int readOffset, writeOffset;
	size_t arrayLength;
	Coordinate newLength =
	    destinationLength + sourceLength + wordLength - 1;
	Descriptor *new;
	Coordinate index;

	// Specify new array
	arrayLength = newLength / 4;
	if (newLength % 4)
		arrayLength++;
	new = callocOrExit(arrayLength, Descriptor);
	for (index = 0; index < arrayLength; index++)
		new[index] = 0;

	writePtr = new;
	writeOffset = 0;

	// Going to end of first descriptor 
	readPtr = &(copy[(sourceLength + wordLength - 2) / 4]);
	readCopy = *readPtr;
	readOffset = (sourceLength + wordLength - 2) % 4;
	readCopy <<= ((3 - readOffset) * 2);

	// Copying reverse complement of first descriptor (minus final (k-1)-mer)
	for (index = 0; index < sourceLength; index++) {
		(*writePtr) >>= 2;
		if (readOffset == 3)
			readCopy = *readPtr;
#ifndef COLOR
		(*writePtr) += 192 - (readCopy & 192);
#else
		(*writePtr) += (readCopy & 192);
#endif
		/*switch (3 - ((readCopy & 192) >> 6)) {
		   case ADENINE:
		   printf("A(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case CYTOSINE:
		   printf("C(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case GUANINE:
		   printf("G(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case THYMINE:
		   printf("T(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   default:
		   printf("?(%ld %i %i);", index, writeOffset, readOffset);
		   } */
		readCopy <<= 2;

		writeOffset++;
		if (writeOffset == 4) {
			writePtr++;
			writeOffset = 0;
		}

		readOffset--;
		if (readOffset == -1) {
			readPtr--;
			readOffset = 3;
		}
	}

	//puts("");

	// Going on copying second descriptor
	readPtr = descr;
	readCopy = *readPtr;
	readOffset = 0;

	for (index = 0; index < destinationLength + wordLength - 1;
	     index++) {
		(*writePtr) >>= 2;
		if (readOffset == 0)
			readCopy = *readPtr;
		(*writePtr) += (readCopy & 3) << 6;
		/*switch ((readCopy & 3)) {
		   case ADENINE:
		   printf("A(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case CYTOSINE:
		   printf("C(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case GUANINE:
		   printf("G(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case THYMINE:
		   printf("T(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   default:
		   printf("?(%ld %i %i);", index, writeOffset, readOffset);
		   } */
		readCopy >>= 2;

		writeOffset++;
		if (writeOffset == 4) {
			writePtr++;
			writeOffset = 0;
		}

		readOffset++;
		if (readOffset == 4) {
			readPtr++;
			readOffset = 0;
		}
	}

	//puts("");

	if (writeOffset != 0) {
		while (writeOffset != 4) {
			(*writePtr) >>= 2;
			writeOffset++;
		}
	}

	return new;
}

void setMultiplicity_pg(PreArc * preArc, IDnum mult)
{
	preArc->multiplicity = mult;
}

static void updatePreArcData_pg(PreArc * preArc, IDnum oldPreNodeID,
				IDnum newPreNodeID)
{
	if (preArc->preNodeIDLeft == oldPreNodeID)
		preArc->preNodeIDLeft = newPreNodeID;
	if (preArc->preNodeIDRight == oldPreNodeID)
		preArc->preNodeIDRight = newPreNodeID;
}

// Reshuffles the preGraph->preNodes array to remove NULL pointers
// Beware that preNode IDs are accordingly reshuffled (all pointers remain valid though)
void renumberPreNodes_pg(PreGraph * preGraph)
{
	IDnum preNodeIndex;
	PreNode *currentPreNode, *destinationPreNode;
	IDnum counter = 0;
	IDnum preNodes = preGraph->preNodeCount;
	IDnum newIndex;
	PreArc *preArc;

	puts("Renumbering preNodes");
	printf("Initial preNode count %d\n", preGraph->preNodeCount);

	for (preNodeIndex = 1; preNodeIndex <= preNodes; preNodeIndex++) {
		currentPreNode = &(preGraph->preNodes[preNodeIndex]);

		if (currentPreNode->descriptor == NULL)
			counter++;
		else if (counter != 0) {
			newIndex = preNodeIndex - counter;
			destinationPreNode =
			    &(preGraph->preNodes[newIndex]);

			destinationPreNode->preArcLeft =
			    currentPreNode->preArcLeft;
			destinationPreNode->preArcRight =
			    currentPreNode->preArcRight;
			destinationPreNode->descriptor =
			    currentPreNode->descriptor;
			destinationPreNode->length =
			    currentPreNode->length;

			for (preArc = getPreArc_pg(newIndex, preGraph);
			     preArc != NULL;
			     preArc = getNextPreArc_pg(preArc, newIndex))
				updatePreArcData_pg(preArc, preNodeIndex,
						    newIndex);
			for (preArc = getPreArc_pg(-newIndex, preGraph);
			     preArc != NULL;
			     preArc = getNextPreArc_pg(preArc, -newIndex))
				updatePreArcData_pg(preArc, -preNodeIndex,
						    -newIndex);
		}
	}

	preGraph->preNodeCount -= counter;
	preGraph->preNodes = reallocOrExit(preGraph->preNodes,
				     preGraph->preNodeCount +
				      1, PreNode);

	printf("Destroyed %d preNodes\n", counter);
}

// Allocate memory for an empty preGraph created with sequenceCount different sequences
PreGraph *emptyPreGraph_pg(IDnum sequenceCount, int wordLength, boolean double_strand)
{
	PreGraph *newPreGraph = mallocOrExit(1, PreGraph);
	newPreGraph->sequenceCount = sequenceCount;
	newPreGraph->wordLength = wordLength;
	newPreGraph->preNodeCount = 0;
	newPreGraph->double_strand = double_strand;
	return newPreGraph;
}

static Descriptor *newDescriptor_pg(Coordinate length, FILE * file,
				    Kmer * initialKmer, int wordLength)
{
	char letter;
	Nucleotide nucleotide;
	Coordinate totalLength = length + wordLength - 1;
	size_t arrayLength = totalLength / 4;
	Descriptor *res;
	Coordinate index;
	Kmer kmerCopy;

	if (totalLength % 4 > 0)
		arrayLength++;

	res = callocOrExit(arrayLength, Descriptor);

	copyKmers(&kmerCopy, initialKmer);
	for (index = wordLength - 2; index >= 0; index--)
		writeNucleotideInDescriptor_pg(popNucleotide(&kmerCopy), res,
					       index);

	for (index = wordLength - 1; index < totalLength; index++) {
		letter = getc(file);
		while (!isalpha(letter))
			letter = getc(file);

		//printf("%c", letter);
		switch (letter) {
		case 'A':
			nucleotide = ADENINE;
			break;
		case 'C':
			nucleotide = CYTOSINE;
			break;
		case 'G':
			nucleotide = GUANINE;
			break;
		case 'T':
			nucleotide = THYMINE;
			break;
		default:
			fflush(stdout);
			abort();
		}

		writeNucleotideInDescriptor_pg(nucleotide, res, index);
		pushNucleotide(initialKmer, nucleotide);
	}

	//printf(" ");

	return res;
}

void allocatePreNodeSpace_pg(PreGraph * preGraph, IDnum preNodeCount)
{
	preGraph->preNodes = callocOrExit(preNodeCount + 1, PreNode);
	preGraph->preNodeCount = preNodeCount;
}

void addPreNodeToPreGraph_pg(PreGraph * preGraph, Coordinate start,
			     Coordinate finish, FILE * file,
			     Kmer * initialKmer, IDnum ID)
{
	PreNode *newnd = &(preGraph->preNodes[ID]);

	newnd->preArcLeft = NULL;
	newnd->preArcRight = NULL;

	newnd->length = finish - start;

	newnd->descriptor =
	    newDescriptor_pg(newnd->length, file, initialKmer,
			     preGraph->wordLength);
}

static void exportPreNode_pg(FILE * outfile, PreNode * preNode, IDnum ID,
			     int wordLength)
{
	Coordinate index;
	Nucleotide nucleotide;

	if (preNode == NULL)
		return;

	fprintf(outfile, "NODE\t%ld\t%lld\n", (long) ID, (long long) preNode->length);

	if (preNode->length == 0) {
		fprintf(outfile, "\n");
		return;
	}

	for (index = 0; index < preNode->length + wordLength - 1; index++) {
		nucleotide =
		    getNucleotideInDescriptor_pg(preNode->descriptor,
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

void exportPreGraph_pg(char *filename, PreGraph * preGraph)
{
	IDnum index;
	FILE *outfile;
	PreNode *preNode;
	int wordLength = getWordLength_pg(preGraph);

	if (preGraph == NULL) {
		return;
	}

	outfile = fopen(filename, "w");
	if (outfile == NULL) {
		puts("Couldn't open file, sorry");
		return;
	} else
		printf("Writing into pregraph file %s...\n", filename);

	// General data
	fprintf(outfile, "%ld\t%ld\t%i\t%hi\n", (long) preGraph->preNodeCount,
		(long) preGraph->sequenceCount, preGraph->wordLength, (short) preGraph->double_strand);

	// PreNode info
	for (index = 1; index <= preGraph->preNodeCount; index++) {
		preNode = getPreNodeInPreGraph_pg(preGraph, index);
		exportPreNode_pg(outfile, preNode, index, wordLength);
	}

	fclose(outfile);
}

int getWordLength_pg(PreGraph * preGraph)
{
	return preGraph->wordLength;
}

void displayPreArcMemory_pg()
{
	if (preArcMemory == NULL)
		return;
	printf("ARC MEMORY %lld allocated %lld free\n",
	       (long long) RecycleBin_memory_usage(preArcMemory),
	       (long long) recycleBinFreeSpace(preArcMemory));
}

boolean hasSinglePreArc_pg(IDnum preNodeID, PreGraph * preGraph)
{
	IDnum ID = preNodeID;
	PreNode *preNode;
	PreArc *preArc;

	if (ID < 0)
		ID = -ID;

	preNode = &(preGraph->preNodes[ID]);

	if (preNodeID > 0)
		preArc = preNode->preArcRight;
	else
		preArc = preNode->preArcLeft;

	return (preArc != NULL
		&& getNextPreArc_pg(preArc, preNodeID) == NULL);
}

char simplePreArcCount_pg(IDnum preNodeID, PreGraph * preGraph)
{
	PreNode *preNode;
	PreArc *preArc;
	char count = 0;
	IDnum ID = preNodeID;

	if (ID < 0)
		ID = -ID;

	preNode = &(preGraph->preNodes[ID]);

	if (preNodeID > 0)
		preArc = preNode->preArcRight;
	else
		preArc = preNode->preArcLeft;

	for (; preArc != NULL;
	     preArc = getNextPreArc_pg(preArc, preNodeID))
		count++;

	return count;
}

boolean isLoop_pg(PreArc * preArc)
{
	return (preArc->preNodeIDLeft == preArc->preNodeIDRight
		|| preArc->preNodeIDLeft == -preArc->preNodeIDRight);
}

void setPreNodeDescriptor_pg(Descriptor * descr, Coordinate length, IDnum preNodeID, PreGraph * preGraph) {
	PreNode * preNode;

	if (preNodeID < 0) 
		preNodeID = -preNodeID;

	preNode = getPreNodeInPreGraph_pg(preGraph, preNodeID);
	free(preNode->descriptor);
	preNode->descriptor = descr;
	preNode->length = length;	
}

static void appendPositiveDescriptor_pg(Descriptor ** writePtr, int * writeOffset, IDnum preNodeID, PreGraph * preGraph, boolean initial) {
	PreNode * preNode = getPreNodeInPreGraph_pg(preGraph, preNodeID);
	Descriptor * readPtr = preNode->descriptor;
	Descriptor readCopy;
	int wordLength = getWordLength_pg(preGraph);
	Coordinate length = preNode->length;
	Coordinate index;
	int readOffset = 0;

	if (initial) {
		index = 0;
		readPtr = preNode->descriptor;
		readCopy = *readPtr;
		readOffset = 0;
	} else {
		index = wordLength - 1;
		readPtr = &(preNode->descriptor[(wordLength - 1) / 4]);
		readCopy = *readPtr;
		readOffset = (wordLength - 1) % 4;
		readCopy >>= (readOffset * 2);
	}

	for (; index < length + wordLength - 1; index++) {
		(**writePtr) >>= 2;
		if (readOffset == 0)
			readCopy = *readPtr;
		(**writePtr) += (readCopy & 3) << 6;
		readCopy >>= 2;

		if (++(*writeOffset) == 4) {
			(*writePtr)++;
			*writeOffset = 0;
		}

		if (++readOffset == 4) {
			readPtr++;
			readOffset = 0;
		}
	}
}

static void appendNegativeDescriptor_pg(Descriptor ** writePtr, int * writeOffset, IDnum preNodeID, PreGraph * preGraph, boolean initial) {
	PreNode * preNode = getPreNodeInPreGraph_pg(preGraph, preNodeID);
	Descriptor * readPtr = preNode->descriptor;
	Descriptor readCopy;
	int wordLength = getWordLength_pg(preGraph);
	Coordinate length = preNode->length;
	Coordinate index;
	int readOffset;

	if (initial) 
		length += wordLength - 1;

	readPtr = &(preNode->descriptor[(length - 1) / 4]);
	readCopy = *readPtr;
	readOffset = (length - 1) % 4;
	readCopy <<= ((3 - readOffset) * 2);

	for (index = 0; index < length; index++) {
		(**writePtr) >>= 2;
		if (readOffset == 3)
			readCopy = *readPtr;
#ifndef COLOR
		(**writePtr) += 192 - (readCopy & 192);
#else
		(**writePtr) += (readCopy & 192);
#endif
		readCopy <<= 2;

		(*writeOffset)++;
		if (*writeOffset == 4) {
			(*writePtr)++;
			*writeOffset = 0;
		}

		readOffset--;
		if (readOffset == -1) {
			readPtr--;
			readOffset = 3;
		}
	}
}

void appendDescriptors_pg(Descriptor ** start, int * writeOffset, IDnum preNodeID, PreGraph* preGraph, boolean initial) {
	if (preNodeID > 0)
		appendPositiveDescriptor_pg(start, writeOffset, preNodeID, preGraph, initial);
	else
		appendNegativeDescriptor_pg(start, writeOffset, -preNodeID, preGraph, initial);
}
