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
#ifndef _GRAPH_H_
#define _GRAPH_H_

////////////////////////////////////////////////////////////
// Node functions
////////////////////////////////////////////////////////////

//Creators/destructor
Node *newNode(IDnum sequenceID, Coordinate start, Coordinate finish,
	      Coordinate offset, IDnum ID, TightString ** sequences,
	      int WORDLENGTH);
Node *emptyNode();
void destroyNode(Node * node, Graph * graph);

// Locator
IDnum getNodeID(Node * node);
Node *getNodeInGraph(Graph * graph, IDnum nodeID);
Node *getTwinNode(Node * node);

// Arc info 
int arcCount(Node * node);
int simpleArcCount(Node * node);
Arc *getArc(Node * node);
boolean hasSingleArc(Node * node);

// Descriptor
Coordinate getNodeLength(Node * node);
void appendDescriptors(Node * target, Node * source);
void directlyAppendDescriptors(Node * target, Node * sourcei, Coordinate totalLength);
void appendSequence(Node * node, TightString ** reads,
		    PassageMarker * guide, Graph * graph);
void clipNodeLength(Node * node, Coordinate startClip,
		    Coordinate finishClip);
void splitNodeDescriptor(Node * source, Node * target, Coordinate offset);
void reduceNode(Node * node);
void reallocateNodeDescriptor(Node * node, Coordinate length);
Nucleotide getNucleotideInNode(Node * node, Coordinate index);

// Passage markers
void insertPassageMarker(PassageMarker * marker, Node * destination);
PassageMarker *getMarker(Node * node);
void setMarker(Node * node, PassageMarker * marker);
IDnum markerCount(Node * node);

// Short read marker creation
void incrementReadStartCount(Node * node, Graph * graph);
void addReadStart(Node * node, IDnum seqID, Coordinate position,
		  Graph * graph, Coordinate offset);
void blurLastShortReadMarker(Node * node, Graph * graph);

// Short read marker handling
ShortReadMarker *getNodeReads(Node * node, Graph * graph);
IDnum getNodeReadCount(Node * node, Graph * graph);
ShortReadMarker *commonNodeReads(Node * nodeA, Node * nodeB, Graph * graph,
				 IDnum * length);
ShortReadMarker *extractBackOfNodeReads(Node * node, Coordinate breakpoint,
					Graph * graph, IDnum * length,
					PassageMarker * sourceMarker,
					Coordinate * sequenceLengths);
ShortReadMarker *extractFrontOfNodeReads(Node * node,
					 Coordinate breakpoint,
					 Graph * graph, IDnum * length,
					 PassageMarker * sourceMarker,
					 Coordinate * sequenceLengths);

// Short read marker moving around
void foldSymmetricalNodeReads(Node * node, Graph * graph);
void spreadReadIDs(ShortReadMarker * reads, IDnum readCount, Node * node,
		   Graph * graph);
void injectShortReads(ShortReadMarker * sourceArray, IDnum sourceLength,
		      Node * target, Graph * graph);
void mergeNodeReads(Node * target, Node * source, Graph * graph);
void checkNodeReads(IDnum nodeIndex, Graph * graph);

// Virtual coverage
void setVirtualCoverage(Node * node, Category category,
			Coordinate coverage);
void incrementVirtualCoverage(Node * node, Category category,
			      Coordinate coverage);
Coordinate getVirtualCoverage(Node * node, Category category);

// Original virtual coverage
void setOriginalVirtualCoverage(Node * node, Category category,
				Coordinate coverage);
void incrementOriginalVirtualCoverage(Node * node, Category category,
				      Coordinate coverage);
Coordinate getOriginalVirtualCoverage(Node * node, Category category);

// Status
void setNodeStatus(Node * node, boolean status);
void setSingleNodeStatus(Node * node, boolean status);
boolean getNodeStatus(Node * node);

// Uniqueness
void setUniqueness(Node * node, boolean value);
boolean getUniqueness(Node * node);

// Gap markers
void appendGap(Node * node, Coordinate length, Graph * graph);
void appendNodeGaps(Node * destination, Node * source, Graph * graph);

// IO
char *readNode(Node * node);
TightString *expandNode(Node * node, int WORDLENGTH);
void appendNodeSequence(Node * node, TightString * sequence,
			Coordinate writeIndex);
char *expandNodeFragment(Node * node, Coordinate contigStart,
			 Coordinate contigFinish, int WORDLENGTH);

////////////////////////////////////////////////////////////
// Arc functions
////////////////////////////////////////////////////////////

// Creators/destructor
Arc *createArc(Node * origin, Node * destination, Graph * graph);
void createAnalogousArc(Node * origin, Node * destination, Arc * refArc,
			Graph * graph);
void destroyArc(Arc * arc, Graph * graph);

// Multiplicity
void setMultiplicity(Arc * arc, IDnum mult);
IDnum getMultiplicity(Arc * arc);
void changeMultiplicity(Arc * arc, IDnum variation);

// Extremities
Node *getOrigin(Arc * arc);
Node *getDestination(Arc * arc);

// Finding arcs
Arc *getArcBetweenNodes(Node * origin, Node * destination, Graph * graph);
Arc *getNextArc(Arc * arc);

// Lookup table option
void activateArcLookupTable(Graph * graph);
void deactivateArcLookupTable(Graph * graph);

////////////////////////////////////////////////////////////
// Short read marker functions 
////////////////////////////////////////////////////////////

ShortReadMarker *getShortReadMarkerAtIndex(ShortReadMarker * array,
					   IDnum index);

IDnum getShortReadMarkerID(ShortReadMarker * marker);

extern inline Coordinate getShortReadMarkerPosition(ShortReadMarker * marker);
extern inline void setShortReadMarkerPosition(ShortReadMarker * marker,
				       Coordinate position);

extern inline ShortLength getShortReadMarkerOffset(ShortReadMarker * marker);
extern inline void setShortReadMarkerOffset(ShortReadMarker * marker,
				     ShortLength offset);

////////////////////////////////////////////////////////////
// Gap marker functions 
////////////////////////////////////////////////////////////

GapMarker *getGap(Node * node, Graph * graph);
GapMarker *getNextGap(GapMarker * marker);
Coordinate getGapStart(GapMarker * marker);
Coordinate getGapFinish(GapMarker * marker);

////////////////////////////////////////////////////////////
// Graph functions
////////////////////////////////////////////////////////////

// Memory allocation
Graph *emptyGraph(IDnum sequenceCount, int wordLength);
void allocateNodeSpace(Graph * graph, IDnum nodeCount);
void addNodeToGraph(Graph * graph, Node * node);
Node *addEmptyNodeToGraph(Graph * graph, IDnum nodeID);
void destroyGraph(Graph * graph);

// Dimensions
IDnum nodeCount(Graph * graph);
IDnum sequenceCount(Graph * graph);
void renumberNodes(Graph * graph);
int getWordLength(Graph * graph);

// Element status 
void resetNodeStatus(Graph * graph);
void resetPassageMarkersStatus(Graph * graph);

// Arc mults
void reassessArcMultiplicities(Graph * graph);

// File IO
void displayGraph(Graph * graph);
Graph *importGraph(char *filename);
Graph *importSimplifiedGraph(char *filename);
void exportGraph(char *filename, Graph * graph, TightString ** sequences);
void exportDOTGraph(char *filename, Graph * graph);
Graph *readPreGraphFile(char *preGraphFilename, boolean * double_strand);

// Read starts
void activateReadStarts(Graph * graph);
boolean readStartsAreActivated(Graph * graph);
void createNodeReadStartArrays(Graph * graph);
void orderNodeReadStartArrays(Graph * graph);

// Insert lengths
void setInsertLengths(Graph * graph, Category cat, Coordinate insertLength,
		      Coordinate insertLength_std_dev);
Coordinate getInsertLength(Graph * graph, Category cat);
double getInsertLength_var(Graph * graph, Category cat);

// Gaps markers
void activateGapMarkers(Graph * graph);
void deactivateGapMarkers(Graph * graph);
void sortGapMarkers(Graph * graph);

void displayArcMemory();
void displayNodeMemory();

void checkPassageMarkersStatus(Graph * graph);
#endif
