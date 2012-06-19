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

//Creators/destructor
PreNode *newPreNode_pg(Coordinate start,
		       Coordinate finish,
		       FILE * file, Kmer * initialKmer, int wordLength);
void destroyPreNode_pg(IDnum preNode, PreGraph * preGraph);

// Locator
PreNode *getPreNodeInPreGraph_pg(PreGraph * preGraph, IDnum preNodeID);

// PreArc info 
PreArc *getPreArc_pg(IDnum preNodeID, PreGraph * preGraph);
boolean hasSinglePreArc_pg(IDnum preNodeID, PreGraph * graph);
char simplePreArcCount_pg(IDnum preNodeID, PreGraph * preGraph);

// Descriptor
Coordinate getPreNodeLength_pg(IDnum preNodeID, PreGraph * preGraph);
void setPreNodeDescriptor_pg(Descriptor * descr, Coordinate length, IDnum preNodeID, PreGraph * preGraph);
void appendDescriptors_pg(Descriptor ** start, int * writeOffset, IDnum preNodeID, PreGraph* preGraph, boolean initial);

////////////////////////////////////////////////////////////
// PreArc functions
////////////////////////////////////////////////////////////

// Creators/destructor
PreArc *createPreArc_pg(IDnum originID, IDnum destinationID,
			PreGraph * preGraph);
void createAnalogousPreArc_pg(IDnum originID, IDnum destinationID,
			      PreArc * refPreArc, PreGraph * preGraph);
void destroyPreArc_pg(PreArc * preArc, PreGraph * preGraph);

// Multiplicity
void setMultiplicity_pg(PreArc * preArc, IDnum mult);
IDnum getMultiplicity_pg(PreArc * preArc);
void changeMultiplicity_pg(PreArc * preArc, IDnum variation);

// Extremities
IDnum getDestination_pg(PreArc * preArc, IDnum nodeID);
IDnum getOtherEnd_pg(PreArc * preArc, IDnum preNodeID);

// Finding preArcs
PreArc *getPreArcBetweenPreNodes_pg(IDnum originID, IDnum destinationID,
				    PreGraph * preGraph);
PreArc *getNextPreArc_pg(PreArc * preArc, IDnum originPreNodeID);

// Misc
boolean isLoop_pg(PreArc * preArc);

////////////////////////////////////////////////////////////
// PreGraph functions
////////////////////////////////////////////////////////////

// Memory allocation
PreGraph *emptyPreGraph_pg(IDnum sequenceCount, int wordLength, boolean double_strand);
void allocatePreNodeSpace_pg(PreGraph * preGraph, IDnum preNodeCount);
void addPreNodeToPreGraph_pg(PreGraph * preGraph, Coordinate start,
			     Coordinate stop, FILE * file,
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

void displayPreArcMemory_pg();

int test_preGraph(int argc, char **argv);
#endif
