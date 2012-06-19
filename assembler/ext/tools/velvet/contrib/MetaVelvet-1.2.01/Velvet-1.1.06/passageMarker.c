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
#ifdef _OPENMP
#include <omp.h>
#endif

#include "globals.h"
#include "allocArray.h"
#include "graph.h"
#include "recycleBin.h"
#include "passageMarker.h"
#include "tightString.h"
#include "utility.h"

typedef struct passage_st PassageMarker;

struct passage_st {
	struct node_st *node;
	PassageMarkerI nextInNode;
	PassageMarkerI previousInNode;
	PassageMarkerI twinMarker;
	PassageMarkerI nextInSequence;
	IDnum start;
	IDnum finishOffset;
	IDnum sequenceID;
	boolean status;
} ATTRIBUTE_PACKED;

static AllocArray *markerMemory = NULL;
DECLARE_FAST_ACCESSORS (PM, PassageMarker, markerMemory)

static RecycleBin *listMemory = NULL;
static const int LISTBLOCKSIZE = 10000;

PassageMarkerI allocatePassageMarker()
{
  if (markerMemory == NULL){
    char passagemarker[50] = "PassageMarker";
    markerMemory =
      newAllocArray (sizeof(PassageMarker), passagemarker);
  }

  return allocArrayAllocate (markerMemory);
}

static void deallocatePassageMarker(PassageMarkerI marker)
{
	allocArrayFree(markerMemory, marker);
}

PassageMarkerList *allocatePassageMarkerList()
{
	if (listMemory == NULL)
		listMemory =
		    newRecycleBin(sizeof(PassageMarkerList),
				  LISTBLOCKSIZE);

	return (PassageMarkerList *) allocatePointer(listMemory);
}

void deallocatePassageMarkerList(PassageMarkerList * marker)
{
	deallocatePointer(listMemory, marker);
}

void setNextInSequence(PassageMarkerI previous, PassageMarkerI next)
{
	if (previous == NULL_IDX)
		return;

	PM_FI2P (previous)->nextInSequence = next;
}

void extractPassageMarker(PassageMarkerI marker)
{
	PassageMarker *twin;
	PassageMarker *markerVal;

	if (marker == NULL_IDX)
		return;

	markerVal = PM_FI2P (marker);
	if (markerVal->node == NULL_IDX)
		return;

	if (markerVal->previousInNode == marker)
		setMarker(markerVal->node, markerVal->nextInNode);
	else
		setNextInNode(markerVal->previousInNode, markerVal->nextInNode);

	markerVal->previousInNode = NULL_IDX;
	markerVal->nextInNode = NULL_IDX;
	markerVal->node = NULL_IDX;

	twin = PM_FI2P (markerVal->twinMarker);
	twin->nextInNode = NULL_IDX;
	twin->previousInNode = NULL_IDX;
	twin->node = NULL_IDX;
}

void destroyPassageMarker(PassageMarkerI marker)
{
	PassageMarker *markerVal;
	PassageMarker *twinVal;
	PassageMarkerI twin;

	if (marker == NULL_IDX)
		return;

	markerVal = PM_FI2P (marker);
	twin = markerVal->twinMarker;
	extractPassageMarker(marker);

	if (markerVal->nextInSequence != NULL_IDX
	    && PM_FI2P (PM_FI2P (markerVal->nextInSequence)->twinMarker)->nextInSequence == twin)
		PM_FI2P (PM_FI2P (markerVal->nextInSequence)->twinMarker)->nextInSequence = NULL_IDX;

	twinVal = PM_FI2P (twin);
	if (twinVal->nextInSequence != NULL_IDX
	    && PM_FI2P (PM_FI2P (twinVal->nextInSequence)->twinMarker)->nextInSequence == marker)
		PM_FI2P (PM_FI2P (twinVal->nextInSequence)->twinMarker)->nextInSequence = NULL_IDX;

	deallocatePassageMarker(twin);
	deallocatePassageMarker(marker);

	//velvetLog("Done destroying passage marker\n");
}

void destroyAllPassageMarkers()
{
	if (markerMemory != NULL)
		destroyAllocArray(markerMemory);
	if (listMemory != NULL)
		destroyRecycleBin(listMemory);
}


void setPreviousInSequence(PassageMarkerI previous, PassageMarkerI next)
{
	if (next == NULL_IDX)
		return;
	else if (previous == NULL_IDX)
		PM_FI2P (PM_FI2P (next)->twinMarker)->nextInSequence = NULL_IDX;
	else
		PM_FI2P (PM_FI2P (next)->twinMarker)->nextInSequence = PM_FI2P (previous)->twinMarker;
}

void disconnectNextPassageMarker(PassageMarkerI marker, Graph * graph)
{
	PassageMarkerI middle = getNextInSequence(marker);
	PassageMarkerI next = getNextInSequence(middle);

	setPreviousInSequence(marker, next);
	concatenatePassageMarkers(marker, middle);
	setNextInSequence(middle, NULL_IDX);
	setPreviousInSequence(NULL_IDX, middle);
}

void deleteNextPassageMarker(PassageMarkerI marker, Graph * graph)
{
	PassageMarkerI middle = getNextInSequence(marker);
	PassageMarkerI next = getNextInSequence(middle);

	setPreviousInSequence(marker, next);
	setNextInSequence(marker, next);
	setNextInSequence(middle, NULL_IDX);
	setPreviousInSequence(NULL_IDX, middle);
}

PassageMarkerI getNextInNode(PassageMarkerI marker)
{
	if (marker == NULL_IDX)
		return NULL_IDX;

	return PM_FI2P (marker)->nextInNode;
}

void setNextInNode(PassageMarkerI marker, PassageMarkerI next)
{
	PassageMarker *markerVal;

	// DEBUG
	if (next == marker || next == getTwinMarker(marker))
		abort();

	if (marker == NULL_IDX)
		return;

	markerVal = PM_FI2P (marker);
	if (next == NULL_IDX) {
		markerVal->nextInNode = NULL_IDX;
		PM_FI2P (markerVal->twinMarker)->nextInNode = NULL_IDX;
	} else {
		PassageMarker *nextVal;

		if (markerVal->twinMarker == NULL_IDX) {
			velvetLog("Dead marker in node %li %li\n",
			       (long) getNodeID(getNode(marker)),
			       (long) getPassageMarkerSequenceID(marker));
			abort();
		}
		nextVal = PM_FI2P (next);
		markerVal->nextInNode = next;
		PM_FI2P (markerVal->twinMarker)->nextInNode = nextVal->twinMarker;
		nextVal->previousInNode = marker;
		PM_FI2P (nextVal->twinMarker)->previousInNode = markerVal->twinMarker;
	}
}

void setTopOfTheNode(PassageMarkerI marker)
{
	if (marker == NULL_IDX)
		return;

	PM_FI2P (marker)->previousInNode = marker;
}

PassageMarkerI getNextInSequence(PassageMarkerI marker)
{
	if (marker != NULL_IDX)
	{
		PassageMarker *markerVal;

		markerVal = PM_FI2P (marker);
		if (markerVal->nextInSequence == NULL_IDX)
			return NULL_IDX;
		return markerVal->nextInSequence;
	}
	return NULL_IDX;
}

PassageMarkerI getPreviousInSequence(PassageMarkerI marker)
{
	PassageMarker *twinVal;

	if (marker == NULL_IDX)
		return NULL_IDX;

	twinVal = PM_FI2P (PM_FI2P (marker)->twinMarker);
	if (twinVal->nextInSequence == NULL_IDX)
		return NULL_IDX;

	return PM_FI2P (twinVal->nextInSequence)->twinMarker;
}

void
connectPassageMarkers(PassageMarkerI previous, PassageMarkerI next,
		      Graph * graph)
{
	if (previous != NULL_IDX)
		setNextInSequence(previous, next);

	if (next != NULL_IDX)
		setPreviousInSequence(previous, next);
}

char *readPassageMarker(PassageMarkerI marker)
{
	PassageMarker *markerVal;
	char *s = mallocOrExit(100, char);

	if (marker == NULL_IDX)
		return s;

	markerVal = PM_FI2P (marker);
	sprintf(s, "MARKER %ld (%lld -> %lld):", (long) markerVal->sequenceID,
		(long long) markerVal->start, (long long) getPassageMarkerFinish(marker));

	if (getPreviousInSequence(marker) == NULL_IDX)
		sprintf(s, "%s START -> %ld", s,
			(long) getNodeID(getNode(marker)));
	else
		sprintf(s, "%s %ld -> %ld", s,
			(long) getNodeID(getNode(getPreviousInSequence(marker))),
			(long) getNodeID(getNode(marker)));

	if (getNextInSequence(marker) == NULL_IDX)
		sprintf(s, "%s -> FINISH", s);
	else
		sprintf(s, "%s -> %ld ", s,
			(long) getNodeID(getNode(getNextInSequence(marker))));

	return s;
}

PassageMarkerI addPassageMarker(IDnum sequenceID, Coordinate start,
				Node * node)
{
	PassageMarkerI marker = allocatePassageMarker();
	PassageMarkerI twinMarker = allocatePassageMarker();
	PassageMarker *markerVal;
	PassageMarker *twinVal;

	markerVal = PM_FI2P (marker);
	twinVal = PM_FI2P (twinMarker);

	markerVal->sequenceID = sequenceID;
	markerVal->start = start;
	markerVal->node = node;
	markerVal->nextInSequence = NULL_IDX;
	markerVal->finishOffset = 0;
	markerVal->twinMarker = twinMarker;
	markerVal->status = false;

	twinVal->sequenceID = -sequenceID;
	twinVal->start = start + getNodeLength(node);
	twinVal->node = getTwinNode(node);
	twinVal->nextInSequence = NULL_IDX;
	twinVal->finishOffset = 0;
	twinVal->twinMarker = marker;
	twinVal->status = false;

	setNextInNode(marker, getMarker(node));
	setMarker(node, marker);

	return marker;
}

PassageMarkerList *copyPassageMarkerList(PassageMarkerList * list)
{
	PassageMarkerList *copy;
	PassageMarkerList *result = NULL;
	PassageMarkerList *pointer;

	if (list == NULL)
		return NULL;

	for (pointer = list; pointer != NULL; pointer = pointer->next) {
		copy = allocatePassageMarkerList();
		copy->marker = pointer->marker;
		copy->next = result;
		result = copy;
	}

	return result;
}

void incrementFinishOffset(PassageMarkerI marker, Coordinate offset)
{
	PM_FI2P (marker)->finishOffset += offset;
}

void incrementStartOffset(PassageMarkerI marker, Coordinate offset)
{
	PM_FI2P (PM_FI2P (marker)->twinMarker)->finishOffset += offset;
}

Coordinate getFinishOffset(PassageMarkerI marker)
{
	return PM_FI2P (marker)->finishOffset;
}

void setFinishOffset(PassageMarkerI marker, Coordinate offset)
{
	PM_FI2P (marker)->finishOffset = offset;
}

Coordinate getStartOffset(PassageMarkerI marker)
{
	return PM_FI2P (PM_FI2P (marker)->twinMarker)->finishOffset;
}

void setStartOffset(PassageMarkerI marker, Coordinate offset)
{
	PM_FI2P (PM_FI2P (marker)->twinMarker)->finishOffset = offset;
}

void transposePassageMarker(PassageMarkerI marker, Node * node)
{
	PassageMarker *markerVal;
	PassageMarker *twinMarkerVal;

	markerVal = PM_FI2P (marker);
	twinMarkerVal = PM_FI2P (markerVal->twinMarker);
	insertPassageMarker(marker, node);
	markerVal->node = node;
	insertPassageMarker(markerVal->twinMarker, getTwinNode(node));
	twinMarkerVal->node = getTwinNode(node);
}

PassageMarkerI getTwinMarker(PassageMarkerI marker)
{
	return PM_FI2P (marker)->twinMarker;
}

IDnum getPassageMarkerSequenceID(PassageMarkerI marker)
{
	return PM_FI2P (marker)->sequenceID;
}

IDnum getAbsolutePassMarkerSeqID(PassageMarkerI marker)
{
	IDnum ID = PM_FI2P (marker)->sequenceID;

	if (ID > 0)
		return ID;
	else
		return -ID;
}

Node *getNode(PassageMarkerI marker)
{
	if (marker == NULL_IDX)
		return NULL;

	return PM_FI2P (marker)->node;
}

void concatenatePassageMarkers(PassageMarkerI marker,
			       PassageMarkerI next)
{
	PassageMarker *markerVal;
	PassageMarker *nextVal;

	if (marker == NULL_IDX || next == NULL_IDX)
		return;

	markerVal = PM_FI2P (marker);
	nextVal = PM_FI2P (next);

	markerVal->finishOffset = nextVal->finishOffset;
	PM_FI2P (markerVal->twinMarker)->start = PM_FI2P (nextVal->twinMarker)->start;
	markerVal->nextInSequence = nextVal->nextInSequence;
}

boolean getPassageMarkerStatus(PassageMarkerI marker)
{
	return PM_FI2P (marker)->status;
}

void setPassageMarkerStatus(PassageMarkerI marker, boolean status)
{
	PassageMarker *markerVal = PM_FI2P (marker);

	markerVal->status = status;
	PM_FI2P (markerVal->twinMarker)->status = status;
}

boolean isDestinationToMarker(PassageMarkerI marker, Node * node)
{
	PassageMarker *markerVal = PM_FI2P (marker);

	if (markerVal->nextInSequence == NULL_IDX)
		return false;

	return PM_FI2P (markerVal->nextInSequence)->node == node;
}

boolean isTerminal(PassageMarkerI marker)
{
	PassageMarker *markerVal;

	if (marker == NULL_IDX)
		return false;

	markerVal = PM_FI2P (marker);
	return markerVal->nextInSequence == NULL_IDX;
}

boolean isInitial(PassageMarkerI marker)
{
	PassageMarker *markerVal;

	if (marker == NULL_IDX)
		return false;

	markerVal = PM_FI2P (marker);
	if (markerVal->twinMarker == NULL_IDX) {
		velvetLog("Unpaired marker seq %ld start %lld node %ld\n",
		       (long) markerVal->sequenceID, (long long) markerVal->start,
		       (long) getNodeID(markerVal->node));
		velvetLog("SNAFU\n");
		abort();
	}

	return PM_FI2P (markerVal->twinMarker)->nextInSequence == NULL_IDX;
}

Coordinate getPassageMarkerStart(PassageMarkerI marker)
{
	return PM_FI2P (marker)->start;
}

void setPassageMarkerStart(PassageMarkerI marker, Coordinate start)
{
	PM_FI2P (marker)->start = start;
}

Coordinate getPassageMarkerFinish(PassageMarkerI marker)
{
	PassageMarker *twinMarkerVal;

	twinMarkerVal = PM_FI2P (PM_FI2P (marker)->twinMarker);
	if (twinMarkerVal->start == -10)
		return -10;

	return twinMarkerVal->start;
}

void setPassageMarkerFinish(PassageMarkerI marker, Coordinate finish)
{
	PassageMarker *twinMarkerVal;

	twinMarkerVal = PM_FI2P (PM_FI2P (marker)->twinMarker);
	if (finish == -10)
		twinMarkerVal->start = -10;

	twinMarkerVal->start = finish;
}

Coordinate getPassageMarkerLength(PassageMarkerI marker)
{
	PassageMarker *markerVal;
	PassageMarker *twinMarkerVal;

	markerVal = PM_FI2P (marker);
	twinMarkerVal = PM_FI2P (markerVal->twinMarker);
	if (markerVal->start == -10 || twinMarkerVal->start == -10)
		return 0;

	else if (markerVal->sequenceID > 0)
		return twinMarkerVal->start - markerVal->start;
	else
		return markerVal->start - twinMarkerVal->start;
}

int passageMarkerDirection(PassageMarkerI marker)
{
	if (PM_FI2P (marker)->sequenceID > 0)
		return 1;
	else
		return -1;
}

PassageMarkerI addUncertainPassageMarker(IDnum sequenceID, Node * node)
{
	PassageMarkerI marker = allocatePassageMarker();
	PassageMarkerI twinMarker = allocatePassageMarker();
	PassageMarker *markerVal = PM_FI2P (marker);
	PassageMarker *twinMarkerVal = PM_FI2P (twinMarker);

	markerVal->sequenceID = sequenceID;
	markerVal->start = -10;
	markerVal->node = node;
	markerVal->nextInSequence = NULL_IDX;
	markerVal->finishOffset = 0;
	markerVal->twinMarker = twinMarker;
	markerVal->status = false;

	twinMarkerVal->sequenceID = -sequenceID;
	twinMarkerVal->start = -10;
	if (node == NULL)
		twinMarkerVal->node = NULL;
	else
		twinMarkerVal->node = getTwinNode(node);
	twinMarkerVal->nextInSequence = NULL_IDX;
	twinMarkerVal->finishOffset = 0;
	twinMarkerVal->twinMarker = marker;
	twinMarkerVal->status = false;

	if (node != NULL) {
		setNextInNode(marker, getMarker(node));
		setMarker(node, marker);
	}

	return marker;
}

PassageMarkerList *newPassageMarkerList(PassageMarkerI marker,
					PassageMarkerList * next)
{
	PassageMarkerList *list = allocatePassageMarkerList();
	list->marker = marker;
	list->next = next;
	return list;
}

PassageMarkerI newPassageMarker(IDnum seqID, Coordinate start,
				Coordinate finish, Coordinate startOffset,
				Coordinate finishOffset)
{
	PassageMarkerI marker;
	PassageMarkerI twinMarker;
	PassageMarker *markerVal;
	PassageMarker *twinMarkerVal;

#ifdef _OPENMP
	#pragma omp critical
	{
#endif
		marker = allocatePassageMarker();
		twinMarker = allocatePassageMarker();
#ifdef _OPENMP
	}
#endif
	markerVal = PM_FI2P (marker);
	twinMarkerVal = PM_FI2P (twinMarker);

//      velvetLog("Values %d\t%d\t%d\t%d\t%d\n", seqID, start, finish, startOffset, finishOffset);

	markerVal->sequenceID = seqID;
	markerVal->node = NULL;
	markerVal->nextInSequence = NULL_IDX;
	markerVal->twinMarker = twinMarker;
	markerVal->nextInNode = NULL_IDX;
	markerVal->status = false;

	twinMarkerVal->sequenceID = -seqID;
	twinMarkerVal->node = NULL;
	twinMarkerVal->nextInSequence = NULL_IDX;
	twinMarkerVal->twinMarker = marker;
	twinMarkerVal->nextInNode = NULL_IDX;
	twinMarkerVal->status = false;

	setPassageMarkerStart(marker, start);
	setPassageMarkerFinish(marker, finish);
	setStartOffset(marker, startOffset);
	setFinishOffset(marker, finishOffset);

	if (getPassageMarkerLength(marker) < 0) {
		velvetLog("Negative marker %ld %lld %lld %lld\n",
		       (long) getPassageMarkerSequenceID(marker),
		       (long long) getPassageMarkerStart(marker),
		       (long long) getPassageMarkerFinish(marker),
		       (long long) getPassageMarkerLength(marker));
		abort();
	}

	return marker;
}

void exportMarker(FILE * outfile, PassageMarkerI marker,
		  TightString * sequences, int WORDLENGTH)
{
	PassageMarker *markerVal = PM_FI2P (marker);
	PassageMarkerI current;

	if (markerVal->sequenceID > 0) {
		if (!isInitial(marker)) {
			return;
		}
		current = marker;
	} else {
		if (!isTerminal(marker)) {
			return;
		}
		current = markerVal->twinMarker;
	}

	velvetFprintf(outfile, "SEQ\t%li\n", (long) PM_FI2P (current)->sequenceID);
	for (; current != NULL_IDX; current = PM_FI2P (current)->nextInSequence) {
		velvetFprintf(outfile, "%ld\t%lld\t%lld\t%lld\t%lld",
			(long) getNodeID(PM_FI2P (current)->node), (long long) getStartOffset(current),
			(long long) getPassageMarkerStart(current),
			(long long) getPassageMarkerFinish(current),
			(long long) getFinishOffset(current));
		velvetFprintf(outfile, "\n");
	}
}
