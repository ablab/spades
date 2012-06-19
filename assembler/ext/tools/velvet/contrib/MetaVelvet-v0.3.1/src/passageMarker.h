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
#ifndef _PASSAGEMARKER_H_
#define _PASSAGEMARKER_H_

struct passageList_st {
	PassageMarker *marker;
	PassageMarkerList *next;
};

///////////////////////////////////////////////////////////////////
// PassageMarker lists
///////////////////////////////////////////////////////////////////
// You can always malloc a PassaegMarkerList but these routines manage the 
// memory for you, thus avoiding fragmentation
PassageMarkerList *newPassageMarkerList(PassageMarker * marker,
					PassageMarkerList * next);

void deallocatePassageMarkerList(PassageMarkerList * list);

///////////////////////////////////////////////////////////////////
// Creators/Destructors
///////////////////////////////////////////////////////////////////
PassageMarker *addPassageMarker(IDnum sequenceID, Coordinate start,
				Node * node);

PassageMarker *copyPassageMarker(PassageMarker * marker);

PassageMarker *addUncertainPassageMarker(IDnum sequenceID, Node * node);

PassageMarker *newPassageMarker(IDnum seqID, Coordinate start,
				Coordinate finish, Coordinate startOffset,
				Coordinate finishOffset);

// Deallocates but also removes all pointers towards that structure
void destroyPassageMarker(PassageMarker * marker);
void destroyAllPassageMarkers();

///////////////////////////////////////////////////////////////////
// Node 
///////////////////////////////////////////////////////////////////

// Current node
Node *getNode(PassageMarker * marker);

// Yank out of current node
void extractPassageMarker(PassageMarker * marker);

// Insert into a node
void transposePassageMarker(PassageMarker * marker, Node * destination);

///////////////////////////////////////////////////////////////////
// General Info
///////////////////////////////////////////////////////////////////
// Export into file
void exportMarker(FILE * outfile, PassageMarker * marker,
		  TightString ** sequences, int wordLength);

// General info for debugging
char *readPassageMarker(PassageMarker * marker);

// String description
char *readPassageMarkerSequence(PassageMarker * marker,
				TightString ** sequences, int WORDLENGTH);

// Sequence ID associated to the passage marker
IDnum getPassageMarkerSequenceID(PassageMarker * marker);
IDnum getAbsolutePassMarkerSeqID(PassageMarker * marker);
int passageMarkerDirection(PassageMarker * marker);

// Coordinates
Coordinate getPassageMarkerStart(PassageMarker * marker);
void setPassageMarkerStart(PassageMarker * marker, Coordinate start);
Coordinate getPassageMarkerFinish(PassageMarker * marker);
void setPassageMarkerFinish(PassageMarker * marker, Coordinate finish);
Coordinate getPassageMarkerLength(PassageMarker * marker);

// Offsets
Coordinate getStartOffset(PassageMarker * marker);
void setStartOffset(PassageMarker * marker, Coordinate offset);
void incrementStartOffset(PassageMarker * marker, Coordinate offset);
Coordinate getFinishOffset(PassageMarker * marker);
void setFinishOffset(PassageMarker * marker, Coordinate offset);
void incrementFinishOffset(PassageMarker * marker, Coordinate offset);

// Status
void setPassageMarkerStatus(PassageMarker * marker, boolean status);
boolean getPassageMarkerStatus(PassageMarker * marker);

///////////////////////////////////////////////////////////////////
// Marker Sequences
///////////////////////////////////////////////////////////////////

// Corresponding marker of reverse complement sequence
PassageMarker *getTwinMarker(PassageMarker * marker);

// Within a node
PassageMarker *getNextInNode(PassageMarker * marker);
void setNextInNode(PassageMarker * marker, PassageMarker * next);
void setTopOfTheNode(PassageMarker * marker);

// Within a sequence
PassageMarker *getNextInSequence(PassageMarker * marker);
void setNextInSequence(PassageMarker * previous, PassageMarker * next);
PassageMarker *getPreviousInSequence(PassageMarker * marker);
void setPreviousInSequence(PassageMarker * previous, PassageMarker * next);
void connectPassageMarkers(PassageMarker * previous, PassageMarker * next,
			   Graph * graph);

// End of read chains
boolean isTerminal(PassageMarker * marker);
boolean isInitial(PassageMarker * marker);

// Checks whether the node of the next marker is the one given in parameter
boolean isDestinationToMarker(PassageMarker * marker, Node * node);

// Bypasses the middle marker
void disconnectNextPassageMarker(PassageMarker * marker, Graph * graph);

// Merge two markers (cf concatenateGraph())
void concatenatePassageMarkers(PassageMarker * marker,
			       PassageMarker * nextMarker);

#endif
