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
#ifndef _PREGRAPH_H_
#define _PREGRAPH_H_

////////////////////////////////////////////////////////////
// PreNode functions
////////////////////////////////////////////////////////////

void destroyPreNode_pg(IDnum preNode, PreGraph * preGraph);

// Locator
PreNode *getPreNodeInPreGraph_pg(PreGraph * preGraph, IDnum preNodeID);

// PreArc info 
PreArcI getPreArc_pg(IDnum preNodeID, PreGraph * preGraph);
boolean hasSinglePreArc_pg(IDnum preNodeID, PreGraph * graph);
char simplePreArcCount_pg(IDnum preNodeID, PreGraph * preGraph);

// Descriptor
Coordinate getPreNodeLength_pg(IDnum preNodeID, PreGraph * preGraph);
void setPreNodeDescriptor_pg(Descriptor * descr, Coordinate length, IDnum preNodeID, PreGraph * preGraph);
void appendDescriptors_pg(Descriptor ** start, int * writeOffset, IDnum preNodeID, PreGraph* preGraph, boolean initial);

////////////////////////////////////////////////////////////
// PreMarker functions
////////////////////////////////////////////////////////////

boolean referenceMarkersAreActivated_pg(PreGraph * preGraph);
void allocatePreMarkerCountSpace_pg(PreGraph * preGraph);
void incrementNodeReferenceMarkerCount_pg(PreGraph * preGraph, IDnum preNodeID);
void allocatePreMarkerSpace_pg(PreGraph * preGraph);
PreMarker * addPreMarker_pg(PreGraph * preGraph, IDnum nodeID, IDnum seqID, Coordinate * start, PreMarker * previous);
void concatenateReferenceMarkers_pg(IDnum preNodeAID, IDnum preNodeBID, PreGraph * preGraph, Coordinate totalOffset);
boolean hasPreMarkers(IDnum nodeID, PreGraph * preGraph);

////////////////////////////////////////////////////////////
// PreArc functions
////////////////////////////////////////////////////////////

// Creators/destructor
PreArcI createPreArc_pg(IDnum originID, IDnum destinationID,
			PreGraph * preGraph);
void createAnalogousPreArc_pg(IDnum originID, IDnum destinationID,
			      PreArcI refPreArc, PreGraph * preGraph);
void destroyPreArc_pg(PreArcI preArc, PreGraph * preGraph);

// Multiplicity
void setMultiplicity_pg(PreArcI preArc, IDnum mult);
IDnum getMultiplicity_pg(PreArcI preArc);

// Extremities
IDnum getDestination_pg(PreArcI preArc, IDnum nodeID);
IDnum getOtherEnd_pg(PreArcI preArc, IDnum preNodeID);

// Finding preArcs
PreArcI getPreArcBetweenPreNodes_pg(IDnum originID, IDnum destinationID,
				    PreGraph * preGraph);
PreArcI getNextPreArc_pg(PreArcI preArc, IDnum originPreNodeID);

// Misc
boolean isLoop_pg(PreArcI preArc);

////////////////////////////////////////////////////////////
// PreGraph functions
////////////////////////////////////////////////////////////

// Memory allocation
PreGraph *emptyPreGraph_pg(IDnum sequenceCount, IDnum referenceCount, int wordLength, boolean double_strand);
void allocatePreNodeSpace_pg(PreGraph * preGraph, IDnum preNodeCount);
void addPreNodeToPreGraph_pg(PreGraph * preGraph, Coordinate start,
			     Coordinate stop, SequencesReader *seqReadInfo,
			     Kmer * initialKmer, IDnum ID);

// Deallocation
void destroyPreGraph_pg(PreGraph * preGraph);

// Dimensions
IDnum preNodeCount_pg(PreGraph * preGraph);
IDnum sequenceCount_pg(PreGraph * preGraph);
void renumberPreNodes_pg(PreGraph * preGraph);

// File IO
void exportPreGraph_pg(char *filename, PreGraph * preGraph);

int getWordLength_pg(PreGraph * preGraph);

#endif
