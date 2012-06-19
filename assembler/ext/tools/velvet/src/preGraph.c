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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "globals.h"
#include "allocArray.h"
#include "preGraph.h"
#include "recycleBin.h"
#include "tightString.h"
#include "run.h"
#include "utility.h"

#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

struct preMarker_st {
	PreMarker * previous;
	PreMarker * next;
	IDnum referenceStart;
	IDnum preNodeStart;
	IDnum length;
	IDnum referenceID;
	IDnum preNodeID; /* SF TODO only the sign seems to matter. Could replace with char or bit field */
} ATTRIBUTE_PACKED;

typedef struct preArc_st PreArc;

struct preArc_st {
	PreArcI nextLeft; /* Index of the previous PreArc */
	PreArcI nextRight; /* Index of the next PreArc */
	IDnum multiplicity;
	IDnum preNodeIDLeft;
	IDnum preNodeIDRight;
} ATTRIBUTE_PACKED;

struct preNode_st {
	PreArcI preArcLeft;
	PreArcI preArcRight;
	Descriptor *descriptor;
	IDnum length;
}  ATTRIBUTE_PACKED;

struct preGraph_st {
	PreNode *preNodes;
	IDnum * nodeReferenceMarkerCounts;
	PreMarker ** nodeReferenceMarkers;
	IDnum sequenceCount;
	IDnum referenceCount;
	IDnum preNodeCount;
	int wordLength;
	boolean double_strand;
};

static AllocArray *preArcMemory = NULL;

DECLARE_FAST_ACCESSORS(PREARC, PreArc, preArcMemory)

PreArcI allocatePreArc_pg()
{
#ifdef _OPENMP
	return allocArrayArrayAllocate (preArcMemory); 
#else
	if (preArcMemory == NULL)
		preArcMemory = newAllocArray(sizeof(PreArc), "PreArc");
	return allocArrayAllocate (preArcMemory);
#endif

}

void deallocatePreArc_pg(PreArcI preArc)
{
#ifdef _OPENMP
	allocArrayArrayFree (preArcMemory, preArc);
#else
	allocArrayFree (preArcMemory, preArc);
#endif
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

PreArcI getPreArcBetweenPreNodes_pg(IDnum originPreNodeID,
				    IDnum destinationPreNodeID,
				    PreGraph * preGraph)
{
	PreArcI preArc;

	if (originPreNodeID == 0 || destinationPreNodeID == 0) {
		return NULL_IDX;
	}

	for (preArc = getPreArc_pg(originPreNodeID, preGraph);
	     preArc != NULL_IDX;
	     preArc = getNextPreArc_pg(preArc, originPreNodeID)) {
		if (getDestination_pg(preArc, originPreNodeID) ==
		    destinationPreNodeID) {
			return preArc;
		}
	}

	return NULL_IDX;
}

static void addPreArcToPreNode_pg(PreArcI preArc, IDnum preNodeID,
				  PreGraph * preGraph)
{
	IDnum ID = preNodeID;
	PreNode *preNode;
	PreArcI *preArcPtr;
	PreArc *preArcVal;

	if (ID < 0)
		ID = -ID;

	preNode = &(preGraph->preNodes[ID]);

	if (preNodeID > 0)
		preArcPtr = &(preNode->preArcRight);
	else
		preArcPtr = &(preNode->preArcLeft);

	preArcVal = PREARC_I2P (preArc);
	preArcVal = PREARC_I2P (preArc);

	if (preNodeID == preArcVal->preNodeIDLeft) {
		preArcVal->nextLeft = *preArcPtr;
		*preArcPtr = preArc;
	}

	if (preNodeID == preArcVal->preNodeIDRight) {
		preArcVal->nextRight = *preArcPtr;
		*preArcPtr = preArc;
	}
}

// Creates an preArc from preNode origin to preNode destination.
// If this preArc already exists, increments its multiplicity by 1.
PreArcI createPreArc_pg(IDnum originPreNodeID, IDnum destinationPreNodeID,
			PreGraph * preGraph)
{
	PreArcI preArc;
	PreArc *preArcVal;


	if (originPreNodeID == 0 || destinationPreNodeID == 0)
		return NULL_IDX;

	preArc =
	    getPreArcBetweenPreNodes_pg(originPreNodeID,
					destinationPreNodeID, preGraph);

	if (preArc != NULL_IDX) {
		PREARC_FI2P (preArc)->multiplicity++;
		if (destinationPreNodeID == -originPreNodeID) 
			PREARC_FI2P (preArc)->multiplicity++;
		return preArc;
	}
	// If not found
	preArc = allocatePreArc_pg();
	preArcVal = PREARC_FI2P (preArc);
	preArcVal->preNodeIDLeft = originPreNodeID;
	preArcVal->preNodeIDRight = -destinationPreNodeID;
	preArcVal->multiplicity = 1;

	addPreArcToPreNode_pg(preArc, originPreNodeID, preGraph);

	// Hairpin case
	if (destinationPreNodeID == -originPreNodeID) {
		preArcVal->multiplicity++;
		return preArc;
	}

	addPreArcToPreNode_pg(preArc, -destinationPreNodeID, preGraph);

	return preArc;
}

void createAnalogousPreArc_pg(IDnum originPreNodeID,
			      IDnum destinationPreNodeID,
			      PreArcI refPreArc, PreGraph * preGraph)
{
	PreArcI preArc;
	PreArc *preArcVal;

	if (originPreNodeID == 0 || destinationPreNodeID == 0)
		return;

	preArc =
	    getPreArcBetweenPreNodes_pg(originPreNodeID,
					destinationPreNodeID, preGraph);

	if (preArc != NULL_IDX) {
		PREARC_FI2P (preArc)->multiplicity += PREARC_FI2P (refPreArc)->multiplicity;
		return;
	}
	// If not found
	preArc = allocatePreArc_pg();
	preArcVal = PREARC_FI2P (preArc);
	preArcVal->preNodeIDLeft = originPreNodeID;
	preArcVal->preNodeIDRight = -destinationPreNodeID;
	preArcVal->multiplicity = PREARC_FI2P (refPreArc)->multiplicity;

	addPreArcToPreNode_pg(preArc, originPreNodeID, preGraph);

	// Hairpin case
	if (destinationPreNodeID == -originPreNodeID)
		return;

	addPreArcToPreNode_pg(preArc, -destinationPreNodeID, preGraph);
}

static void setNextPreArc_pg(PreArcI preArc, IDnum preNodeID,
			     PreArcI nextPreArc)
{
	PreArc *preArcVal;

	preArcVal = PREARC_FI2P (preArc);
	if (preNodeID == preArcVal->preNodeIDLeft)
		preArcVal->nextLeft = nextPreArc;
	if (preNodeID == preArcVal->preNodeIDRight)
		preArcVal->nextRight = nextPreArc;
}

void removePreArcFromList_pg(PreArcI preArc, IDnum preNodeID,
			     PreGraph * preGraph)
{
	IDnum ID = preNodeID;
	PreNode *preNode;
	PreArcI *preArcPtr;
	PreArcI tempPreArc;

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

	for (tempPreArc = *preArcPtr; tempPreArc != NULL_IDX;
	     tempPreArc = getNextPreArc_pg(tempPreArc, preNodeID))
		if (getNextPreArc_pg(tempPreArc, preNodeID) == preArc)
			setNextPreArc_pg(tempPreArc, preNodeID,
					 getNextPreArc_pg(preArc,
							  preNodeID));
}

void destroyPreArc_pg(PreArcI preArc, PreGraph * preGraph)
{
	IDnum leftID, rightID;
	PreArc *preArcVal;

	if (preArc == NULL_IDX)
		return;

	preArcVal = PREARC_FI2P (preArc);
	leftID = preArcVal->preNodeIDLeft;
	rightID = preArcVal->preNodeIDRight;

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
	IDnum index;
	PreMarker * preMarker;

	//velvetLog("Destroying %ld\n", (long) preNodeID);

	if (ID < 0)
		ID = -ID;

	preNode = &(preGraph->preNodes[ID]);

	// PreNode preArcs:
	while (preNode->preArcLeft != NULL_IDX)
		destroyPreArc_pg(preNode->preArcLeft, preGraph);
	while (preNode->preArcRight != NULL_IDX)
		destroyPreArc_pg(preNode->preArcRight, preGraph);

	// PreMarkers
	if (preGraph->nodeReferenceMarkers) {
		for (index = 0; index < preGraph->nodeReferenceMarkerCounts[ID]; index++) {
			preMarker = &(preGraph->nodeReferenceMarkers[ID][index]);
			if (preMarker->previous != NULL)
				preMarker->previous->next = NULL;
			if (preMarker->next != NULL)
				preMarker->next->previous = NULL;
			preMarker->preNodeID = 0;
			preMarker->referenceID = 0;
		}
		if (preGraph->nodeReferenceMarkers[ID]) 
			free(preGraph->nodeReferenceMarkers[ID]);
		preGraph->nodeReferenceMarkers[ID] = NULL;
		preGraph->nodeReferenceMarkerCounts[ID] = 0;
	}

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
#ifdef _OPENMP
	destroyAllocArrayArray(preArcMemory);
#else
	destroyAllocArray(preArcMemory);
#endif

	// Nodes
	free(preGraph->preNodes);

        // PreMarkers
        if (preGraph->nodeReferenceMarkerCounts) {
                free(preGraph->nodeReferenceMarkerCounts);
                free(preGraph->nodeReferenceMarkers);
        }

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

PreArcI getPreArc_pg(IDnum preNodeID, PreGraph * preGraph)
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

PreArcI getNextPreArc_pg(PreArcI preArc, IDnum preNodeID)
{
	PreArc *preArcVal;

	preArcVal = PREARC_FI2P (preArc);

	if (preNodeID == preArcVal->preNodeIDLeft) {
		return preArcVal->nextLeft;
	} else {
		return preArcVal->nextRight;
	}
}

IDnum getMultiplicity_pg(PreArcI preArc)
{
	if (preArc == NULL_IDX)
		return 0;

	return PREARC_FI2P (preArc)->multiplicity;
}

IDnum getOtherEnd_pg(PreArcI preArc, IDnum preNodeID)
{
	PreArc *preArcVal;

	preArcVal = PREARC_FI2P (preArc);
	if (preNodeID == preArcVal->preNodeIDLeft)
		return preArcVal->preNodeIDRight;
	else
		return preArcVal->preNodeIDLeft;
}

IDnum getDestination_pg(PreArcI preArc, IDnum preNodeID)
{
	PreArc *preArcVal;

	if (preArc == NULL_IDX)
		return 0;

	preArcVal = PREARC_FI2P (preArc);

	if (preNodeID == preArcVal->preNodeIDLeft)
		return -preArcVal->preNodeIDRight;
	else
		return -preArcVal->preNodeIDLeft;
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
		   velvetLog("A%ld", index);
		   break;
		   case CYTOSINE:
		   velvetLog("C%ld", index);
		   break;
		   case GUANINE:
		   velvetLog("G%ld", index);
		   break;
		   case THYMINE:
		   velvetLog("T%ld", index);
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

	//velvetLog("\n");

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
		   velvetLog("A%ld", index);
		   break;
		   case CYTOSINE:
		   velvetLog("C%ld", index);
		   break;
		   case GUANINE:
		   velvetLog("G%ld", index);
		   break;
		   case THYMINE:
		   velvetLog("T%ld", index);
		   break;
		   default:
		   velvetLog("?%ld;", index);
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

	//velvetLog("\n");

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
		   velvetLog("A(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case CYTOSINE:
		   velvetLog("C(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case GUANINE:
		   velvetLog("G(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case THYMINE:
		   velvetLog("T(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   default:
		   velvetLog("?(%ld %i %i);", index, writeOffset, readOffset);
		   } */
		readCopy >>= 2;

		writeOffset++;
		if (writeOffset == 4) {
			writePtr++;
			readPtr++;
			writeOffset = 0;
		}
	}

	//velvetLog("\n");

	// Going to end of second descriptor 
	readPtr = &(copy[(sourceLength - 1) / 4]);
	readCopy = *readPtr;
	readOffset = (sourceLength - 1) % 4;
	readCopy <<= ((3 - readOffset) * 2);

	//velvetLog("Read copy %x\n", readCopy);

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
		   velvetLog("A(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case CYTOSINE:
		   velvetLog("C(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case GUANINE:
		   velvetLog("G(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case THYMINE:
		   velvetLog("T(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   default:
		   velvetLog("?(%ld %i %i);", index, writeOffset, readOffset);
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

	//velvetLog("\n");

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
		   velvetLog("A(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case CYTOSINE:
		   velvetLog("C(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case GUANINE:
		   velvetLog("G(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case THYMINE:
		   velvetLog("T(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   default:
		   velvetLog("?(%ld %i %i);", index, writeOffset, readOffset);
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

	//velvetLog("\n");

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
		   velvetLog("A(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case CYTOSINE:
		   velvetLog("C(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case GUANINE:
		   velvetLog("G(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   case THYMINE:
		   velvetLog("T(%ld %i %i) ", index, writeOffset, readOffset);
		   break;
		   default:
		   velvetLog("?(%ld %i %i);", index, writeOffset, readOffset);
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

	//velvetLog("\n");

	if (writeOffset != 0) {
		while (writeOffset != 4) {
			(*writePtr) >>= 2;
			writeOffset++;
		}
	}

	return new;
}

void setMultiplicity_pg(PreArcI preArc, IDnum mult)
{
	PREARC_FI2P (preArc)->multiplicity = mult;
}

static void updatePreArcData_pg(PreArcI preArc, IDnum oldPreNodeID,
				IDnum newPreNodeID)
{
	PreArc *preArcVal;

	preArcVal = PREARC_FI2P (preArc);
	if (preArcVal->preNodeIDLeft == oldPreNodeID)
		preArcVal->preNodeIDLeft = newPreNodeID;
	if (preArcVal->preNodeIDRight == oldPreNodeID)
		preArcVal->preNodeIDRight = newPreNodeID;
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
	IDnum preMarkerIndex;
	PreMarker * preMarker;
	PreArcI preArc;

	velvetLog("Renumbering preNodes\n");
	velvetLog("Initial preNode count %li\n", (long) preGraph->preNodeCount);

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
			     preArc != NULL_IDX;
			     preArc = getNextPreArc_pg(preArc, newIndex))
				updatePreArcData_pg(preArc, preNodeIndex,
						    newIndex);
			for (preArc = getPreArc_pg(-newIndex, preGraph);
			     preArc != NULL_IDX;
			     preArc = getNextPreArc_pg(preArc, -newIndex))
				updatePreArcData_pg(preArc, -preNodeIndex,
						    -newIndex);
			
			if (preGraph->nodeReferenceMarkers) {
				preGraph->nodeReferenceMarkerCounts[newIndex] = preGraph->nodeReferenceMarkerCounts[preNodeIndex];
				preGraph->nodeReferenceMarkers[newIndex] = preGraph->nodeReferenceMarkers[preNodeIndex];

				for (preMarkerIndex = 0; preMarkerIndex < preGraph->nodeReferenceMarkerCounts[newIndex]; preMarkerIndex++) {
					preMarker = &(preGraph->nodeReferenceMarkers[newIndex][preMarkerIndex]);
					if (preMarker->preNodeID == preNodeIndex)
						preMarker->preNodeID = newIndex;
					else if (preMarker->preNodeID == -preNodeIndex)
						preMarker->preNodeID = -newIndex;
					else 
						abort();
				}
			}
		}
	}

	preGraph->preNodeCount -= counter;
	preGraph->preNodes = reallocOrExit(preGraph->preNodes,
				     preGraph->preNodeCount +
				      1, PreNode);

	velvetLog("Destroyed %li preNodes\n", (long) counter);
}

// Allocate memory for an empty preGraph created with sequenceCount different sequences
PreGraph *emptyPreGraph_pg(IDnum sequenceCount, IDnum referenceCount, int wordLength, boolean double_strand)
{
	PreGraph *newPreGraph = mallocOrExit(1, PreGraph);
	newPreGraph->sequenceCount = sequenceCount;
	newPreGraph->wordLength = wordLength;
	newPreGraph->preNodeCount = 0;
	newPreGraph->double_strand = double_strand;
	newPreGraph->referenceCount = referenceCount;
	newPreGraph->preNodes = NULL;
	newPreGraph->nodeReferenceMarkerCounts = NULL;
	newPreGraph->nodeReferenceMarkers = NULL;

#ifdef _OPENMP
	preArcMemory = newAllocArrayArray(omp_get_max_threads(), sizeof(PreArc), "PreArc");
#endif

	return newPreGraph;
}

static Descriptor *newDescriptor_pg(Coordinate length, SequencesReader *seqReadInfo,
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
		if (seqReadInfo->m_bIsBinary) {
			letter = **seqReadInfo->m_ppCurrString;
			*seqReadInfo->m_ppCurrString += 1;   // increment the pointer
		} else {
			letter = getc(seqReadInfo->m_pFile);
			while (!isalpha(letter))
				letter = getc(seqReadInfo->m_pFile);
		}
		//velvetLog("%c", letter);
		switch (letter) {
		case 'N':
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

	//velvetLog(" ");

	return res;
}

void allocatePreNodeSpace_pg(PreGraph * preGraph, IDnum preNodeCount)
{
	preGraph->preNodes = callocOrExit(preNodeCount + 1, PreNode);
	preGraph->preNodeCount = preNodeCount;
}

void allocatePreMarkerCountSpace_pg(PreGraph * preGraph)
{
	preGraph->nodeReferenceMarkerCounts = callocOrExit(preGraph->preNodeCount + 1, IDnum);
	preGraph->nodeReferenceMarkers = callocOrExit(preGraph->preNodeCount + 1, PreMarker *);
}

void incrementNodeReferenceMarkerCount_pg(PreGraph * preGraph, IDnum preNodeID) {
	if (preNodeID < 0)
		preNodeID = -preNodeID;

	preGraph->nodeReferenceMarkerCounts[preNodeID]++;
}

void allocatePreMarkerSpace_pg(PreGraph * preGraph) {
	IDnum index;
	
	if (!preGraph->nodeReferenceMarkers)
		return;

	for (index = 1; index <= preGraph->preNodeCount; index++) {
		if (preGraph->nodeReferenceMarkerCounts[index])
			preGraph->nodeReferenceMarkers[index] = callocOrExit(preGraph->nodeReferenceMarkerCounts[index], PreMarker);
		else 
			preGraph->nodeReferenceMarkers[index] = NULL;
		preGraph->nodeReferenceMarkerCounts[index] = 0;
	}
}

PreMarker * addPreMarker_pg(PreGraph * preGraph, IDnum nodeID, IDnum seqID, Coordinate * start, PreMarker * previous) {
	PreMarker * preMarker;
	IDnum positive_nodeID;

	if (nodeID < 0) 
		abort();
	else
		positive_nodeID = nodeID;

	//printf("Adding preMarker %li\n", (long) *start);

	preMarker = &(preGraph->nodeReferenceMarkers[positive_nodeID][(preGraph->nodeReferenceMarkerCounts[positive_nodeID])++]);
	preMarker->previous = previous;
	if (previous) 
		previous->next = preMarker;
	preMarker->next = NULL;
	preMarker->referenceStart = *start;
	preMarker->length = preGraph->preNodes[positive_nodeID].length;
	preMarker->preNodeStart = 0;
	preMarker->preNodeID = nodeID;
	preMarker->referenceID = seqID;

	*start += preMarker->length;

	return preMarker;
}
void addPreNodeToPreGraph_pg(PreGraph * preGraph, Coordinate start,
			     Coordinate finish, SequencesReader *seqReadInfo,
			     Kmer * initialKmer, IDnum ID)
{
	PreNode *newnd = &(preGraph->preNodes[ID]);

	newnd->preArcLeft = NULL_IDX;
	newnd->preArcRight = NULL_IDX;

	newnd->length = finish - start;

	newnd->descriptor =
	    newDescriptor_pg(newnd->length, seqReadInfo, initialKmer,
			     preGraph->wordLength);
}

static void exportPreNode_pg(FILE * outfile, PreNode * preNode, IDnum ID,
			     int wordLength)
{
	Coordinate index;
	Nucleotide nucleotide;

	if (preNode == NULL)
		return;

	velvetFprintf(outfile, "NODE\t%ld\t%lld\n", (long) ID, (long long) preNode->length);

	if (preNode->length == 0) {
		velvetFprintf(outfile, "\n");
		return;
	}

	for (index = 0; index < preNode->length + wordLength - 1; index++) {
		nucleotide =
		    getNucleotideInDescriptor_pg(preNode->descriptor,
						 index);
		switch (nucleotide) {
		case ADENINE:
			velvetFprintf(outfile, "A");
			break;
		case CYTOSINE:
			velvetFprintf(outfile, "C");
			break;
		case GUANINE:
			velvetFprintf(outfile, "G");
			break;
		case THYMINE:
			velvetFprintf(outfile, "T");
			break;
		}
	}

	velvetFprintf(outfile, "\n");
}

static void exportPreMarker(FILE * outfile, PreMarker* preMarker) {
	velvetFprintf(outfile, "%li\t%lli\t%lli\t%lli\n", (long) preMarker->preNodeID, (long long) preMarker->preNodeStart, (long long) preMarker->referenceStart, (long long) preMarker->length);
}

static void exportPreReference_pg(FILE * outfile, IDnum refIndex, PreGraph * preGraph) {
	PreMarker * preMarker;
	IDnum nodeID, index;

	velvetFprintf(outfile, "SEQ\t%li\n", (long) refIndex);

	for (nodeID = 1; nodeID <= preGraph->preNodeCount; nodeID++) {
		for (index = 0; index < preGraph->nodeReferenceMarkerCounts[nodeID]; index++) {
			preMarker = &(preGraph->nodeReferenceMarkers[nodeID][index]);
			if (preMarker->referenceID == refIndex && !preMarker->previous) {
				for (;preMarker;preMarker = preMarker->next) {
					exportPreMarker(outfile, preMarker);
				}
			}
		}
	}
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
		velvetLog("Couldn't open file, sorry\n");
		return;
	} else
		velvetLog("Writing into pregraph file %s...\n", filename);

	// General data
	velvetFprintf(outfile, "%ld\t%ld\t%i\t%hi\n", (long) preGraph->preNodeCount,
		(long) preGraph->sequenceCount, preGraph->wordLength, (short) preGraph->double_strand);

	// PreNode info
	for (index = 1; index <= preGraph->preNodeCount; index++) {
		preNode = getPreNodeInPreGraph_pg(preGraph, index);
		exportPreNode_pg(outfile, preNode, index, wordLength);
	}

	// Reference sequence info
	for (index = 1; index <= preGraph->referenceCount; index++)
		exportPreReference_pg(outfile, index, preGraph);


	fclose(outfile);
}

int getWordLength_pg(PreGraph * preGraph)
{
	return preGraph->wordLength;
}

boolean hasSinglePreArc_pg(IDnum preNodeID, PreGraph * preGraph)
{
	IDnum ID = preNodeID;
	PreNode *preNode;
	PreArcI preArc;

	if (ID < 0)
		ID = -ID;

	preNode = &(preGraph->preNodes[ID]);

	if (preNodeID > 0)
		preArc = preNode->preArcRight;
	else
		preArc = preNode->preArcLeft;

	return (preArc != NULL_IDX
		&& getNextPreArc_pg(preArc, preNodeID) == NULL_IDX);
}

char simplePreArcCount_pg(IDnum preNodeID, PreGraph * preGraph)
{
	PreNode *preNode;
	PreArcI preArc;
	char count = 0;
	IDnum ID = preNodeID;

	if (ID < 0)
		ID = -ID;

	preNode = &(preGraph->preNodes[ID]);

	if (preNodeID > 0)
		preArc = preNode->preArcRight;
	else
		preArc = preNode->preArcLeft;

	for (; preArc != NULL_IDX;
	     preArc = getNextPreArc_pg(preArc, preNodeID))
		count++;

	return count;
}

boolean isLoop_pg(PreArcI preArc)
{
	PreArc *preArcVal = PREARC_FI2P (preArc);

	return (preArcVal->preNodeIDLeft == preArcVal->preNodeIDRight
		|| preArcVal->preNodeIDLeft == -preArcVal->preNodeIDRight);
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

boolean referenceMarkersAreActivated_pg(PreGraph * preGraph) {
	return preGraph->nodeReferenceMarkers != NULL;
}

static void copyPreMarker(PreMarker * dest, PreMarker * source, IDnum preNodeAID, PreGraph * preGraph) {
	dest->previous = source->previous;
	dest->next = source->next;

	dest->preNodeStart = source->preNodeStart;
	dest->length = source->length;
	dest->referenceID = source->referenceID;
	dest->referenceStart = source->referenceStart;

	if (source->preNodeID > 0)
		dest->preNodeID = preNodeAID;
	else 
		dest->preNodeID = -preNodeAID;

	if (source->previous)
		source->previous->next = dest;
	if (source->next)
		source->next->previous = dest;

	source->referenceID = 0;
	source->preNodeID = 0;
	source->previous = NULL;
	source->next = NULL;
}

static PreMarker * reallocOrExitReferenceMarkers(PreGraph * preGraph, IDnum preNodeID, IDnum length) {
	PreMarker * array = callocOrExit(length, PreMarker);
	PreMarker * writer = array;
	PreMarker * reader = preGraph->nodeReferenceMarkers[preNodeID];
	IDnum index;

	for (index = 0; index < preGraph->nodeReferenceMarkerCounts[preNodeID]; index++) {
		copyPreMarker(writer, reader, preNodeID, preGraph);	
		writer++;
		reader++;
	} 

	free(preGraph->nodeReferenceMarkers[preNodeID]);

	return array;
}

static void concatenateReferenceMarkers_H2T_pg(IDnum preNodeAID, IDnum preNodeBID, PreGraph * preGraph, Coordinate totalOffset) {
	IDnum index;
	IDnum countA = preGraph->nodeReferenceMarkerCounts[preNodeAID];
	IDnum countB = preGraph->nodeReferenceMarkerCounts[preNodeBID];
	Coordinate lengthA = preGraph->preNodes[preNodeAID].length + totalOffset;
	PreMarker * markerA, *next, *markerB;
	IDnum counter = 0;
	
	for (index = 0 ; index < countA; index++) {
		markerA = &(preGraph->nodeReferenceMarkers[preNodeAID][index]);

		if (markerA->preNodeID > 0)
			next = markerA->next;
		else 
			next = markerA->previous;

		if (!next)
			continue;

		if (markerA->preNodeID == preNodeAID && next->preNodeID != preNodeBID)
			continue;
		if (markerA->preNodeID == -preNodeAID && next->preNodeID != -preNodeBID)
			continue;

		next->referenceID = 0;
		next->preNodeID = 0;

		markerA->length += next->length;
		if (markerA->preNodeID > 0) {
			markerA->next = next->next;
			if (next->next) 
				next->next->previous = markerA;
		} else {
			markerA->previous = next->previous;
			if (next->previous)
				next->previous->next = markerA;
			markerA->referenceStart = next->referenceStart;
		}
		next->next = NULL;
		next->previous = NULL;
	}

	for (index = 0; index < countB; index++) 
		if (preGraph->nodeReferenceMarkers[preNodeBID][index].referenceID)
			counter++;

	if (counter == 0)
		return;

	if (countA)
		preGraph->nodeReferenceMarkers[preNodeAID] = reallocOrExitReferenceMarkers(preGraph, preNodeAID, countA + counter);
	else 
		preGraph->nodeReferenceMarkers[preNodeAID] = callocOrExit(counter, PreMarker);

	for (index = 0; index < countB; index++) {
		markerB = &(preGraph->nodeReferenceMarkers[preNodeBID][index]);
		if (markerB->referenceID) {
			markerA = &(preGraph->nodeReferenceMarkers[preNodeAID][countA++]);
			copyPreMarker(markerA, markerB, preNodeAID, preGraph);
			markerA->preNodeStart += lengthA;
		}
	}

	preGraph->nodeReferenceMarkerCounts[preNodeAID] = countA;
}

static void concatenateReferenceMarkers_H2H_pg(IDnum preNodeAID, IDnum preNodeBID, PreGraph * preGraph, Coordinate totalOffset) {
	IDnum index;
	IDnum countA = preGraph->nodeReferenceMarkerCounts[preNodeAID];
	IDnum countB = preGraph->nodeReferenceMarkerCounts[preNodeBID];
	Coordinate lengthA = preGraph->preNodes[preNodeAID].length + totalOffset;
	Coordinate lengthB = preGraph->preNodes[preNodeBID].length;
	PreMarker * markerA, *next, *markerB;
	IDnum counter = 0;

	for (index = 0 ; index < countA; index++) {
		markerA = &(preGraph->nodeReferenceMarkers[preNodeAID][index]);

		if (markerA->preNodeID > 0)
			next = markerA->next;
		else 
			next = markerA->previous;

		
		if ((!next) 
		    || (markerA->preNodeID == preNodeAID && next->preNodeID != -preNodeBID)
		    || (markerA->preNodeID == -preNodeAID && next->preNodeID != preNodeBID))
			continue;
	
		next->referenceID = 0;
		next->preNodeID = 0;

		markerA->length += next->length;
		if (markerA->preNodeID > 0) {
			markerA->next = next->next;
			if (next->next)
				next->next->previous = markerA;
		} else {
			markerA->previous = next->previous;
			if (next->previous)
				next->previous->next = markerA;
			markerA->referenceStart = next->referenceStart;
		}
		next->next = NULL;
		next->previous = NULL;
	}

	for (index = 0; index < countB; index++) 
		if (preGraph->nodeReferenceMarkers[preNodeBID][index].referenceID)
			counter++;

	if (counter == 0)
		return;

	if (countA)
		preGraph->nodeReferenceMarkers[preNodeAID] = reallocOrExitReferenceMarkers(preGraph, preNodeAID, countA + counter);
	else
		preGraph->nodeReferenceMarkers[preNodeAID] = callocOrExit(counter, PreMarker);

	for (index = 0; index < countB; index++) {
		markerB = &(preGraph->nodeReferenceMarkers[preNodeBID][index]);
		if (markerB->referenceID) {
			markerA = &(preGraph->nodeReferenceMarkers[preNodeAID][countA++]);
			copyPreMarker(markerA, markerB, preNodeAID, preGraph);
			markerA->preNodeID *= -1;
			markerA->preNodeStart = lengthA + lengthB - markerA->preNodeStart - markerA->length;
		}
	}

	preGraph->nodeReferenceMarkerCounts[preNodeAID] = countA;
}

static void concatenateReferenceMarkers_T2T_pg(IDnum preNodeAID, IDnum preNodeBID, PreGraph * preGraph, Coordinate totalOffset) {
	IDnum index;
	IDnum countA = preGraph->nodeReferenceMarkerCounts[preNodeAID];
	IDnum countB = preGraph->nodeReferenceMarkerCounts[preNodeBID];
	Coordinate lengthB = preGraph->preNodes[preNodeBID].length;
	PreMarker * markerA, *next, *markerB;
	IDnum counter = 0;

	for (index = 0 ; index < countA; index++) {
		markerA = &(preGraph->nodeReferenceMarkers[preNodeAID][index]);

		if (markerA->preNodeID < 0)
			next = markerA->next;
		else 
			next = markerA->previous;

		if (!next 
		    || (markerA->preNodeID == preNodeAID && next->preNodeID != -preNodeBID)
		    || (markerA->preNodeID == -preNodeAID && next->preNodeID != preNodeBID)) {
			markerA->preNodeStart += lengthB;
			continue;
		}

		next->referenceID = 0;
		next->preNodeID = 0;

		markerA->length += next->length;
		markerA->preNodeStart = lengthB - next->preNodeStart - next->length;
		if (markerA->preNodeID < 0) {
			markerA->next = next->next;
			if (next->next)
				next->next->previous = markerA;
		} else {
			markerA->previous = next->previous;
			if (next->previous) 
				next->previous->next = markerA;
			markerA->referenceStart = next->referenceStart;
		}
		next->next = NULL;
		next->previous = NULL;
	}

	for (index = 0; index < countB; index++) 
		if (preGraph->nodeReferenceMarkers[preNodeBID][index].referenceID)
			counter++;

	if (counter == 0)
		return;

	if (countA)
		preGraph->nodeReferenceMarkers[preNodeAID] = reallocOrExitReferenceMarkers(preGraph, preNodeAID, countA + counter);
	else
		preGraph->nodeReferenceMarkers[preNodeAID] = callocOrExit(counter, PreMarker);

	for (index = 0; index < countB; index++) {
		markerB = &(preGraph->nodeReferenceMarkers[preNodeBID][index]);
		if (markerB->referenceID) {
			markerA = &(preGraph->nodeReferenceMarkers[preNodeAID][countA++]);
			copyPreMarker(markerA, markerB, preNodeAID, preGraph);
			markerA->preNodeID *= -1;
			markerA->preNodeStart = lengthB - markerA->preNodeStart - markerA->length;
		}
	}

	preGraph->nodeReferenceMarkerCounts[preNodeAID] = countA;
}

static void concatenateReferenceMarkers_T2H_pg(IDnum preNodeAID, IDnum preNodeBID, PreGraph * preGraph, Coordinate totalOffset) {
	IDnum index;
	IDnum countA = preGraph->nodeReferenceMarkerCounts[preNodeAID];
	IDnum countB = preGraph->nodeReferenceMarkerCounts[preNodeBID];
	PreMarker * markerA, *next, *markerB;
	Coordinate lengthB = preGraph->preNodes[preNodeBID].length;
	IDnum counter = 0;
	
	for (index = 0 ; index < countA; index++) {
		markerA = &(preGraph->nodeReferenceMarkers[preNodeAID][index]);

		if (markerA->preNodeID < 0)
			next = markerA->next;
		else 
			next = markerA->previous;

		if (!next
		    || (markerA->preNodeID == preNodeAID && next->preNodeID != preNodeBID)
		    || (markerA->preNodeID == -preNodeAID && next->preNodeID != -preNodeBID)) {
			markerA->preNodeStart += lengthB;
			continue;
		}

		next->referenceID = 0;
		next->preNodeID = 0;

		markerA->length += next->length;
		markerA->preNodeStart = next->preNodeStart;
		if (markerA->preNodeID < 0) {
			markerA->next = next->next;
			if (next->next)
				next->next->previous = markerA;
		} else {
			markerA->previous = next->previous;
			if (next->previous)
				next->previous->next = markerA;
			markerA->referenceStart = next->referenceStart;
		}
		next->next = NULL;
		next->previous = NULL;
	}

	for (index = 0; index < countB; index++) 
		if (preGraph->nodeReferenceMarkers[preNodeBID][index].referenceID)
			counter++;

	if (counter == 0)
		return;

	if (countA)
		preGraph->nodeReferenceMarkers[preNodeAID] = reallocOrExitReferenceMarkers(preGraph, preNodeAID, countA + counter);
	else
		preGraph->nodeReferenceMarkers[preNodeAID] = callocOrExit(counter, PreMarker);

	for (index = 0; index < countB; index++) {
		markerB = &(preGraph->nodeReferenceMarkers[preNodeBID][index]);
		if (markerB->referenceID) {
			markerA = &(preGraph->nodeReferenceMarkers[preNodeAID][countA++]);
			copyPreMarker(markerA, markerB, preNodeAID, preGraph);
		}
	}

	preGraph->nodeReferenceMarkerCounts[preNodeAID] = countA;
}

void concatenateReferenceMarkers_pg(IDnum preNodeAID, IDnum preNodeBID, PreGraph * preGraph, Coordinate totalOffset) {
	if (!referenceMarkersAreActivated_pg(preGraph))
		return;

	if (preNodeAID > 0 && preNodeBID > 0)
		concatenateReferenceMarkers_H2T_pg(preNodeAID, preNodeBID, preGraph, totalOffset);
	else if (preNodeAID > 0) 
		concatenateReferenceMarkers_H2H_pg(preNodeAID, -preNodeBID, preGraph, totalOffset);
	else if (preNodeBID > 0)
		concatenateReferenceMarkers_T2T_pg(-preNodeAID, preNodeBID, preGraph, totalOffset);
	else
		concatenateReferenceMarkers_T2H_pg(-preNodeAID, -preNodeBID, preGraph, totalOffset);
}

boolean hasPreMarkers(IDnum nodeID, PreGraph * preGraph) {
	if (nodeID < 0)
		nodeID = -nodeID;
	return preGraph->nodeReferenceMarkerCounts[nodeID] > 0;
}
