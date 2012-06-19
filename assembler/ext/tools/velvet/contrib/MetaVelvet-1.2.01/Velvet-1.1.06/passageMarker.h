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
	PassageMarkerI marker;
	PassageMarkerList *next;
}  ATTRIBUTE_PACKED;

///////////////////////////////////////////////////////////////////
// PassageMarker lists
///////////////////////////////////////////////////////////////////
// You can always malloc a PassaegMarkerList but these routines manage the 
// memory for you, thus avoiding fragmentation
PassageMarkerList *newPassageMarkerList(PassageMarkerI marker,
					PassageMarkerList * next);

void deallocatePassageMarkerList(PassageMarkerList * list);

///////////////////////////////////////////////////////////////////
// Creators/Destructors
///////////////////////////////////////////////////////////////////
PassageMarkerI addPassageMarker(IDnum sequenceID, Coordinate start,
				Node * node);

PassageMarkerI addUncertainPassageMarker(IDnum sequenceID, Node * node);

PassageMarkerI newPassageMarker(IDnum seqID, Coordinate start,
				Coordinate finish, Coordinate startOffset,
				Coordinate finishOffset);

// Deallocates but also removes all pointers towards that structure
void destroyPassageMarker(PassageMarkerI marker);
void destroyAllPassageMarkers();

///////////////////////////////////////////////////////////////////
// Node 
///////////////////////////////////////////////////////////////////

// Current node
Node *getNode(PassageMarkerI marker);

// Yank out of current node
void extractPassageMarker(PassageMarkerI marker);

// Insert into a node
void transposePassageMarker(PassageMarkerI marker, Node * destination);

///////////////////////////////////////////////////////////////////
// General Info
///////////////////////////////////////////////////////////////////
// Export into file
void exportMarker(FILE * outfile, PassageMarkerI marker,
		  TightString * sequences, int wordLength);

// General info for debugging
char *readPassageMarker(PassageMarkerI marker);

// Sequence ID associated to the passage marker
IDnum getPassageMarkerSequenceID(PassageMarkerI marker);
IDnum getAbsolutePassMarkerSeqID(PassageMarkerI marker);
int passageMarkerDirection(PassageMarkerI marker);

// Coordinates
Coordinate getPassageMarkerStart(PassageMarkerI marker);
void setPassageMarkerStart(PassageMarkerI marker, Coordinate start);
Coordinate getPassageMarkerFinish(PassageMarkerI marker);
void setPassageMarkerFinish(PassageMarkerI marker, Coordinate finish);
Coordinate getPassageMarkerLength(PassageMarkerI marker);

// Offsets
Coordinate getStartOffset(PassageMarkerI marker);
void setStartOffset(PassageMarkerI marker, Coordinate offset);
void incrementStartOffset(PassageMarkerI marker, Coordinate offset);
Coordinate getFinishOffset(PassageMarkerI marker);
void setFinishOffset(PassageMarkerI marker, Coordinate offset);
void incrementFinishOffset(PassageMarkerI marker, Coordinate offset);

// Status
void setPassageMarkerStatus(PassageMarkerI marker, boolean status);
boolean getPassageMarkerStatus(PassageMarkerI marker);

///////////////////////////////////////////////////////////////////
// Marker Sequences
///////////////////////////////////////////////////////////////////

// Corresponding marker of reverse complement sequence
PassageMarkerI getTwinMarker(PassageMarkerI marker);

// Within a node
PassageMarkerI getNextInNode(PassageMarkerI marker);
void setNextInNode(PassageMarkerI marker, PassageMarkerI next);
void setTopOfTheNode(PassageMarkerI marker);

// Within a sequence
PassageMarkerI getNextInSequence(PassageMarkerI marker);
void setNextInSequence(PassageMarkerI previous, PassageMarkerI next);
PassageMarkerI getPreviousInSequence(PassageMarkerI marker);
void setPreviousInSequence(PassageMarkerI previous, PassageMarkerI next);
void connectPassageMarkers(PassageMarkerI previous, PassageMarkerI next,
			   Graph * graph);

// End of read chains
boolean isTerminal(PassageMarkerI marker);
boolean isInitial(PassageMarkerI marker);

// Checks whether the node of the next marker is the one given in parameter
boolean isDestinationToMarker(PassageMarkerI marker, Node * node);

// Bypasses the middle marker
void disconnectNextPassageMarker(PassageMarkerI marker, Graph * graph);
void deleteNextPassageMarker(PassageMarkerI marker, Graph * graph);

// Merge two markers (cf concatenateGraph())
void concatenatePassageMarkers(PassageMarkerI marker,
			       PassageMarkerI nextMarker);

#endif
