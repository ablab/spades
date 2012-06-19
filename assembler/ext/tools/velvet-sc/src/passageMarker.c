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
#include "tightString.h"
#include "utility.h"

struct passage_st {
	struct node_st *node;
	PassageMarker *nextInNode;
	PassageMarker *previousInNode;
	PassageMarker *twinMarker;
	PassageMarker *nextInSequence;
	Coordinate start;
	Coordinate finishOffset;
	IDnum sequenceID;
	boolean status;
};

static RecycleBin *markerMemory = NULL;
static RecycleBin *listMemory = NULL;
static const int MARKERBLOCKSIZE = 1000000;
static const int LISTBLOCKSIZE = 10000;

PassageMarker *allocatePassageMarker()
{
	if (markerMemory == NULL)
		markerMemory =
		    newRecycleBin(sizeof(PassageMarker), MARKERBLOCKSIZE);

	return (PassageMarker *) allocatePointer(markerMemory);
}

static void deallocatePassageMarker(PassageMarker * marker)
{
	deallocatePointer(markerMemory, marker);
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

void setNextInSequence(PassageMarker * previous, PassageMarker * next)
{
	if (previous == NULL)
		return;

	previous->nextInSequence = next;
}

void extractPassageMarker(PassageMarker * marker)
{
	PassageMarker *twin;

	if (marker == NULL)
		return;

	if (marker->node == NULL)
		return;

	if (marker->previousInNode == marker)
		setMarker(marker->node, marker->nextInNode);
	else
		setNextInNode(marker->previousInNode, marker->nextInNode);

	marker->previousInNode = NULL;
	marker->nextInNode = NULL;
	marker->node = NULL;

	twin = marker->twinMarker;
	twin->nextInNode = NULL;
	twin->previousInNode = NULL;
	twin->node = NULL;
}

void destroyPassageMarker(PassageMarker * marker)
{
	PassageMarker *twin = marker->twinMarker;

	if (marker == NULL)
		return;

	extractPassageMarker(marker);

	if (marker->nextInSequence != NULL
	    && marker->nextInSequence->twinMarker->nextInSequence == twin)
		marker->nextInSequence->twinMarker->nextInSequence = NULL;

	if (twin->nextInSequence != NULL
	    && twin->nextInSequence->twinMarker->nextInSequence == marker)
		twin->nextInSequence->twinMarker->nextInSequence = NULL;

	deallocatePassageMarker(twin);
	deallocatePassageMarker(marker);

	//puts("Done destroying passage marker");
}

void destroyAllPassageMarkers()
{
	if (markerMemory != NULL)
		destroyRecycleBin(markerMemory);
	if (listMemory != NULL)
		destroyRecycleBin(listMemory);
}


void setPreviousInSequence(PassageMarker * previous, PassageMarker * next)
{
	if (next == NULL)
		return;
	else if (previous == NULL)
		next->twinMarker->nextInSequence = NULL;
	else
		next->twinMarker->nextInSequence = previous->twinMarker;
}

void disconnectNextPassageMarker(PassageMarker * marker, Graph * graph)
{
	PassageMarker *middle = getNextInSequence(marker);
	PassageMarker *next = getNextInSequence(middle);

	setPreviousInSequence(marker, next);
	concatenatePassageMarkers(marker, middle);
	setNextInSequence(middle, NULL);
	setPreviousInSequence(NULL, middle);
}

PassageMarker *getNextInNode(PassageMarker * marker)
{
	if (marker == NULL)
		return NULL;

	return marker->nextInNode;
}

void setNextInNode(PassageMarker * marker, PassageMarker * next)
{
	if (marker == NULL)
		return;

	if (next == NULL) {
		marker->nextInNode = NULL;
		marker->twinMarker->nextInNode = NULL;
	} else {
		if (marker->twinMarker == NULL) {
			printf("Dead marker in node %d %d\n",
			       getNodeID(getNode(marker)),
			       getPassageMarkerSequenceID(marker));
			abort();
		}
		marker->nextInNode = next;
		marker->twinMarker->nextInNode = next->twinMarker;
		next->previousInNode = marker;
		next->twinMarker->previousInNode = marker->twinMarker;
	}
}

void setTopOfTheNode(PassageMarker * marker)
{
	if (marker == NULL)
		return;

	marker->previousInNode = marker;
}

PassageMarker *getNextInSequence(PassageMarker * marker)
{
	if (marker == NULL || marker->nextInSequence == NULL)
		return NULL;

	return marker->nextInSequence;
}

PassageMarker *getPreviousInSequence(PassageMarker * marker)
{
	if (marker == NULL)
		return NULL;

	if (marker->twinMarker->nextInSequence == NULL)
		return NULL;

	return marker->twinMarker->nextInSequence->twinMarker;
}

void
connectPassageMarkers(PassageMarker * previous, PassageMarker * next,
		      Graph * graph)
{
	if (previous != NULL)
		setNextInSequence(previous, next);

	if (next != NULL)
		setPreviousInSequence(previous, next);
}

char *readPassageMarker(PassageMarker * marker)
{
	char *s = mallocOrExit(100, char);

	if (marker == NULL)
		return s;

	sprintf(s, "MARKER %ld (%lld -> %lld):", (long) marker->sequenceID,
		(long long) marker->start, (long long) getPassageMarkerFinish(marker));

	if (getPreviousInSequence(marker) == NULL)
		sprintf(s, "%s START -> %ld", s,
			(long) getNodeID(getNode(marker)));
	else
		sprintf(s, "%s %ld -> %ld", s,
			(long) getNodeID(getNode(getPreviousInSequence(marker))),
			(long) getNodeID(getNode(marker)));

	if (getNextInSequence(marker) == NULL)
		sprintf(s, "%s -> FINISH", s);
	else
		sprintf(s, "%s -> %ld ", s,
			(long) getNodeID(getNode(getNextInSequence(marker))));

	return s;
}

char *readPassageMarkerSequence(PassageMarker * marker,
				TightString ** sequences, int WORDLENGTH)
{
	TightString *sequence =
	    sequences[getAbsolutePassMarkerSeqID(marker) - 1];
	int i;
	char *s = NULL;

	if (marker == NULL)
		return s;

	s = mallocOrExit(getPassageMarkerLength(marker) + 1, char);

	if (getPassageMarkerSequenceID(marker) > 0)
		for (i = 0; i < getPassageMarkerLength(marker); i++)
			s[i] =
			    getNucleotideChar(getPassageMarkerStart(marker)
					      + i + WORDLENGTH - 1,
					      sequence);
	else
		for (i = 0; i < getPassageMarkerLength(marker); i++)
			s[i] =
			    getInverseNucleotideChar(getPassageMarkerStart
						     (marker) - i - 1,
						     sequence);

	s[getPassageMarkerLength(marker)] = '\0';

	return s;
}

PassageMarker *addPassageMarker(IDnum sequenceID, Coordinate start,
				Node * node)
{
	PassageMarker *marker = allocatePassageMarker();
	PassageMarker *twinMarker = allocatePassageMarker();

	marker->sequenceID = sequenceID;
	marker->start = start;
	marker->node = node;
	marker->nextInSequence = NULL;
	marker->finishOffset = 0;
	marker->twinMarker = twinMarker;
	marker->status = false;

	twinMarker->sequenceID = -sequenceID;
	twinMarker->start = start + getNodeLength(node);
	twinMarker->node = getTwinNode(node);
	twinMarker->nextInSequence = NULL;
	twinMarker->finishOffset = 0;
	twinMarker->twinMarker = marker;
	twinMarker->status = false;

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

PassageMarker *copyPassageMarker(PassageMarker * marker)
{
	PassageMarker *twin = marker->twinMarker;
	PassageMarker *copy = allocatePassageMarker();
	PassageMarker *twinCopy = allocatePassageMarker();

	copy->sequenceID = marker->sequenceID;
	copy->start = marker->start;
	copy->nextInNode = NULL;
	copy->previousInNode = NULL;
	copy->node = NULL;
	copy->nextInSequence = marker->nextInSequence;
	copy->finishOffset = marker->finishOffset;
	copy->status = false;

	twinCopy->sequenceID = twin->sequenceID;
	twinCopy->start = twin->start;
	twinCopy->nextInNode = NULL;
	twinCopy->previousInNode = NULL;
	twinCopy->node = NULL;
	twinCopy->nextInSequence = twin->nextInSequence;
	twinCopy->finishOffset = twin->finishOffset;
	twinCopy->status = false;

	copy->twinMarker = twinCopy;
	twinCopy->twinMarker = copy;

	return copy;
}

void incrementFinishOffset(PassageMarker * marker, Coordinate offset)
{
	marker->finishOffset += offset;
}

void incrementStartOffset(PassageMarker * marker, Coordinate offset)
{
	marker->twinMarker->finishOffset += offset;
}

Coordinate getFinishOffset(PassageMarker * marker)
{
	return marker->finishOffset;
}

void setFinishOffset(PassageMarker * marker, Coordinate offset)
{
	marker->finishOffset = offset;
}

Coordinate getStartOffset(PassageMarker * marker)
{
	return marker->twinMarker->finishOffset;
}

void setStartOffset(PassageMarker * marker, Coordinate offset)
{
	marker->twinMarker->finishOffset = offset;
}

void transposePassageMarker(PassageMarker * marker, Node * node)
{
	marker->node = node;
	marker->twinMarker->node = getTwinNode(node);
	insertPassageMarker(marker, node);
	insertPassageMarker(marker->twinMarker, getTwinNode(node));
}

PassageMarker *getTwinMarker(PassageMarker * marker)
{
	return marker->twinMarker;
}

IDnum getPassageMarkerSequenceID(PassageMarker * marker)
{
	return marker->sequenceID;
}

IDnum getAbsolutePassMarkerSeqID(PassageMarker * marker)
{
	IDnum ID = marker->sequenceID;

	if (ID > 0)
		return ID;
	else
		return -ID;
}

Node *getNode(PassageMarker * marker)
{
	if (marker == NULL)
		return NULL;

	return marker->node;
}

void concatenatePassageMarkers(PassageMarker * marker,
			       PassageMarker * next)
{

	if (marker == NULL || next == NULL)
		return;

	marker->finishOffset = next->finishOffset;
	marker->twinMarker->start = next->twinMarker->start;
	marker->nextInSequence = next->nextInSequence;
}

boolean getPassageMarkerStatus(PassageMarker * marker)
{
	return marker->status;
}

void setPassageMarkerStatus(PassageMarker * marker, boolean status)
{
	marker->status = status;
	marker->twinMarker->status = status;
}

boolean isDestinationToMarker(PassageMarker * marker, Node * node)
{
	if (marker->nextInSequence == NULL)
		return false;

	return marker->nextInSequence->node == node;
}

boolean isTerminal(PassageMarker * marker)
{
	if (marker == NULL)
		return false;

	return marker->nextInSequence == NULL;
}

boolean isInitial(PassageMarker * marker)
{
	if (marker == NULL)
		return false;

	if (marker->twinMarker == NULL) {
		printf("Unpaired marker seq %ld start %lld node %ld\n",
		       (long) marker->sequenceID, (long long) marker->start,
		       (long) getNodeID(marker->node));
		puts("SNAFU");
		abort();
	}

	return marker->twinMarker->nextInSequence == NULL;
}

Coordinate getPassageMarkerStart(PassageMarker * marker)
{
	return marker->start;
}

void setPassageMarkerStart(PassageMarker * marker, Coordinate start)
{
	marker->start = start;
}

Coordinate getPassageMarkerFinish(PassageMarker * marker)
{
	if (marker->twinMarker->start == -10)
		return -10;

	return marker->twinMarker->start;
}

void setPassageMarkerFinish(PassageMarker * marker, Coordinate finish)
{
	if (finish == -10)
		marker->twinMarker->start = -10;

	marker->twinMarker->start = finish;
}

Coordinate getPassageMarkerLength(PassageMarker * marker)
{
	if (marker->start == -10 || marker->twinMarker->start == -10)
		return 0;

	else if (marker->sequenceID > 0)
		return marker->twinMarker->start - marker->start;
	else
		return marker->start - marker->twinMarker->start;
}

int passageMarkerDirection(PassageMarker * marker)
{
	if (marker->sequenceID > 0)
		return 1;
	else
		return -1;
}

PassageMarker *addUncertainPassageMarker(IDnum sequenceID, Node * node)
{
	PassageMarker *marker = allocatePassageMarker();
	PassageMarker *twinMarker = allocatePassageMarker();

	marker->sequenceID = sequenceID;
	marker->start = -10;
	marker->node = node;
	marker->nextInSequence = NULL;
	marker->finishOffset = 0;
	marker->twinMarker = twinMarker;
	marker->status = false;

	twinMarker->sequenceID = -sequenceID;
	twinMarker->start = -10;
	if (node == NULL)
		twinMarker->node = NULL;
	else
		twinMarker->node = getTwinNode(node);
	twinMarker->nextInSequence = NULL;
	twinMarker->finishOffset = 0;
	twinMarker->twinMarker = marker;
	twinMarker->status = false;

	if (node != NULL) {
		setNextInNode(marker, getMarker(node));
		setMarker(node, marker);
	}

	return marker;
}

PassageMarkerList *newPassageMarkerList(PassageMarker * marker,
					PassageMarkerList * next)
{
	PassageMarkerList *list = allocatePassageMarkerList();
	list->marker = marker;
	list->next = next;
	return list;
}

PassageMarker *newPassageMarker(IDnum seqID, Coordinate start,
				Coordinate finish, Coordinate startOffset,
				Coordinate finishOffset)
{
	PassageMarker *marker = allocatePassageMarker();
	PassageMarker *twinMarker = allocatePassageMarker();

//      printf("Values %d\t%d\t%d\t%d\t%d\n", seqID, start, finish, startOffset, finishOffset);

	marker->sequenceID = seqID;
	marker->node = NULL;
	marker->nextInSequence = NULL;
	marker->twinMarker = twinMarker;
	marker->nextInNode = NULL;
	marker->status = false;

	twinMarker->sequenceID = -seqID;
	twinMarker->node = NULL;
	twinMarker->nextInSequence = NULL;
	twinMarker->twinMarker = marker;
	twinMarker->nextInNode = NULL;
	twinMarker->status = false;

	setPassageMarkerStart(marker, start);
	setPassageMarkerFinish(marker, finish);
	setStartOffset(marker, startOffset);
	setFinishOffset(marker, finishOffset);

	if (getPassageMarkerLength(marker) < 0) {
		printf("Negative marker %ld %lld %lld %lld\n",
		       (long) getPassageMarkerSequenceID(marker),
		       (long long) getPassageMarkerStart(marker),
		       (long long) getPassageMarkerFinish(marker),
		       (long long) getPassageMarkerLength(marker));
		abort();
	}

	return marker;
}

void exportMarker(FILE * outfile, PassageMarker * marker,
		  TightString ** sequences, int WORDLENGTH)
{
	PassageMarker *current;

	if (marker->sequenceID > 0) {
		if (!isInitial(marker)) {
			return;
		}
		current = marker;
	} else {
		if (!isTerminal(marker)) {
			return;
		}
		current = marker->twinMarker;
	}

	fprintf(outfile, "SEQ\t%d\n", current->sequenceID);
	for (; current != NULL; current = current->nextInSequence) {
		fprintf(outfile, "%ld\t%lld\t%lld\t%lld\t%lld",
			(long) getNodeID(current->node), (long long) getStartOffset(current),
			(long long) getPassageMarkerStart(current),
			(long long) getPassageMarkerFinish(current),
			(long long) getFinishOffset(current));
		fprintf(outfile, "\n");
	}
}
