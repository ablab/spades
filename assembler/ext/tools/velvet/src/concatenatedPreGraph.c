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
#include "preGraph.h"
#include "utility.h"

// Replaces two consecutive preNodes into a single equivalent preNode
// The extra memory is freed
static void concatenatePreNodes(IDnum preNodeAID, PreArcI oldPreArc,
				PreGraph * preGraph)
{
	IDnum preNodeBID = preNodeAID;
	IDnum currentPreNodeID, nextPreNodeID;
	PreArcI preArc = oldPreArc;
	Coordinate totalLength = 0;
	Coordinate arrayLength;
	Descriptor * descr, * ptr;
	int writeOffset = 0;
	int wordLength = getWordLength_pg(preGraph);
	Coordinate totalOffset = 0;

	//velvetLog("Concatenating nodes %li and %li\n", preNodeAID, preNodeBID);

	while(hasSinglePreArc_pg(preNodeBID, preGraph)
		       &&
		       hasSinglePreArc_pg(getOtherEnd_pg
					  (preArc, preNodeBID),
					  preGraph)
		       && !isLoop_pg(preArc) 
		       && getDestination_pg(preArc, preNodeBID) != preNodeAID) {

		totalLength += getPreNodeLength_pg(preNodeBID, preGraph);
		preNodeBID = getDestination_pg(preArc, preNodeBID);
		preArc = getPreArc_pg(preNodeBID, preGraph);
	}
	totalLength += getPreNodeLength_pg(preNodeBID, preGraph);
	totalLength += wordLength - 1;

	// Reference marker management
	if (referenceMarkersAreActivated_pg(preGraph)) {
		currentPreNodeID = preNodeAID;
		preArc = getPreArc_pg(currentPreNodeID, preGraph);
		for (currentPreNodeID = getDestination_pg(preArc, currentPreNodeID); currentPreNodeID != preNodeBID; currentPreNodeID = getDestination_pg(preArc, currentPreNodeID)) {
			concatenateReferenceMarkers_pg(preNodeAID, currentPreNodeID, preGraph, totalOffset);
			preArc = getPreArc_pg(currentPreNodeID, preGraph);
			totalOffset += getPreNodeLength_pg(currentPreNodeID, preGraph);
		}
		concatenateReferenceMarkers_pg(preNodeAID, currentPreNodeID, preGraph, totalOffset);
	}

	// Descriptor management (preNode)
	arrayLength = totalLength / 4;
	if (totalLength % 4)
		arrayLength++;
	descr = callocOrExit(arrayLength, Descriptor);
	ptr = descr;
	if (preNodeAID > 0) {
		currentPreNodeID = preNodeAID;
		appendDescriptors_pg(&ptr, &writeOffset, currentPreNodeID, preGraph, true);
		preArc = getPreArc_pg(currentPreNodeID, preGraph);
		currentPreNodeID = getDestination_pg(preArc, currentPreNodeID);
		while (currentPreNodeID != preNodeBID) {
			appendDescriptors_pg(&ptr, &writeOffset, currentPreNodeID, preGraph, false);
			preArc = getPreArc_pg(currentPreNodeID, preGraph);
			currentPreNodeID = getDestination_pg(preArc, currentPreNodeID);
		}
		appendDescriptors_pg(&ptr, &writeOffset, currentPreNodeID, preGraph, false);
	} else {
		currentPreNodeID = -preNodeBID;
		appendDescriptors_pg(&ptr, &writeOffset ,currentPreNodeID, preGraph, true);
		preArc = getPreArc_pg(currentPreNodeID, preGraph);
		currentPreNodeID = getDestination_pg(preArc, currentPreNodeID);
		while (currentPreNodeID != -preNodeAID) {
			appendDescriptors_pg(&ptr, &writeOffset ,currentPreNodeID, preGraph, false);
			preArc = getPreArc_pg(currentPreNodeID, preGraph);
			currentPreNodeID = getDestination_pg(preArc, currentPreNodeID);
		}
		appendDescriptors_pg(&ptr, &writeOffset ,currentPreNodeID, preGraph, false);
	}

	if (writeOffset != 0) 
		while (writeOffset++ != 4)
			(*ptr) >>= 2;

	setPreNodeDescriptor_pg(descr, totalLength - wordLength + 1, preNodeAID, preGraph); 

	// Correct preArcs
	for (preArc = getPreArc_pg(preNodeBID, preGraph); preArc != NULL_IDX;
	     preArc = getNextPreArc_pg(preArc, preNodeBID)) {
		if (getDestination_pg(preArc, preNodeBID) != -preNodeBID)
			createAnalogousPreArc_pg(preNodeAID,
						 getDestination_pg(preArc,
								   preNodeBID),
						 preArc, preGraph);
		else
			createAnalogousPreArc_pg(preNodeAID, -preNodeAID,
						 preArc, preGraph);
	}

	// Freeing gobbled preNode
	currentPreNodeID = -preNodeBID;
	while (currentPreNodeID != -preNodeAID) {
		preArc = getPreArc_pg(currentPreNodeID, preGraph);
		nextPreNodeID = getDestination_pg(preArc, currentPreNodeID);
		destroyPreNode_pg(currentPreNodeID, preGraph);
		currentPreNodeID = nextPreNodeID;
	}
}

// Detects sequences that could be simplified through concatentation
// Iterates till preGraph cannot be more simplified
// Useless preNodes are freed from memory and remaining ones are renumbered
void concatenatePreGraph_pg(PreGraph * preGraph)
{
	IDnum preNodeIndex;
	PreArcI preArc;
	PreNode *preNode;

	velvetLog("Concatenation...\n");

	for (preNodeIndex = 1; preNodeIndex < preNodeCount_pg(preGraph);
	     preNodeIndex++) {
		preNode = getPreNodeInPreGraph_pg(preGraph, preNodeIndex);

		if (preNode == NULL)
			continue;

		preArc = getPreArc_pg(preNodeIndex, preGraph);

		while (hasSinglePreArc_pg(preNodeIndex, preGraph)
		       &&
		       hasSinglePreArc_pg(getOtherEnd_pg
					  (preArc, preNodeIndex),
					  preGraph)) {
			if (isLoop_pg(preArc))
				break;
			concatenatePreNodes(preNodeIndex, preArc,
					    preGraph);
			preArc = getPreArc_pg(preNodeIndex, preGraph);
		}

		preArc = getPreArc_pg(-preNodeIndex, preGraph);

		while (hasSinglePreArc_pg(-preNodeIndex, preGraph)
		       &&
		       hasSinglePreArc_pg(getOtherEnd_pg
					  (preArc, -preNodeIndex),
					  preGraph)) {
			if (isLoop_pg(preArc))
				break;
			concatenatePreNodes(-preNodeIndex, preArc,
					    preGraph);
			preArc = getPreArc_pg(-preNodeIndex, preGraph);
		}
	}

	renumberPreNodes_pg(preGraph);
	velvetLog("Concatenation over!\n");
}

static boolean isEligibleTip(IDnum index, PreGraph * preGraph, Coordinate
			     cutoffLength)
{
	IDnum currentIndex = -index;
	Coordinate totalLength = 0;
	PreArcI activeArc = NULL_IDX;
	PreArcI arc;
	IDnum mult = 0;

	if (getPreArc_pg(index, preGraph) != NULL_IDX)
		return false;

	// Finding first tangle
	while (currentIndex != 0
	       && simplePreArcCount_pg(-currentIndex, preGraph) < 2
	       && simplePreArcCount_pg(currentIndex, preGraph) < 2
	       && (!referenceMarkersAreActivated_pg(preGraph) || !hasPreMarkers(currentIndex, preGraph))) {
		totalLength += getPreNodeLength_pg(currentIndex, preGraph);
		activeArc = getPreArc_pg(currentIndex, preGraph);
		currentIndex = getDestination_pg(activeArc, currentIndex);
	}

	// If too long
	if (totalLength >= cutoffLength)
		return false;

	// If isolated snippet:
	if (currentIndex == 0)
		return true;

	// Joined tips  
	if (simplePreArcCount_pg(-currentIndex, preGraph) < 2)
		return false;

	// If unique event
	if (getMultiplicity_pg(activeArc) == 1)
		return true;

	// Computing max arc
	for (arc = getPreArc_pg(-currentIndex, preGraph); arc != NULL_IDX;
	     arc = getNextPreArc_pg(arc, -currentIndex))
		if (getMultiplicity_pg(arc) > mult)
			mult = getMultiplicity_pg(arc);

	// Testing for minority
	return mult > getMultiplicity_pg(activeArc);
}

void clipTips_pg(PreGraph * preGraph)
{
	IDnum index;
	PreNode *current;
	boolean modified = true;
	Coordinate cutoffLength = getWordLength_pg(preGraph) * 2;
	IDnum counter = 0;

	velvetLog("Clipping short tips off preGraph\n");

	while (modified) {
		modified = false;
		for (index = 1; index <= preNodeCount_pg(preGraph);
		     index++) {
			current = getPreNodeInPreGraph_pg(preGraph, index);

			if (current == NULL)
				continue;

			if (isEligibleTip(index, preGraph, cutoffLength)
			    || isEligibleTip(-index, preGraph,
					     cutoffLength)) {
				counter++;
				destroyPreNode_pg(index, preGraph);
				modified = true;
			}
		}
	}

	concatenatePreGraph_pg(preGraph);
	velvetLog("%li tips cut off\n", (long) counter);
	velvetLog("%li nodes left\n", (long) preNodeCount_pg(preGraph));
}
