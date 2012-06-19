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
#include "roadMap.h"
#include "readSet.h"
#include "concatenatedPreGraph.h"
#include "utility.h"
#include "kmer.h"

#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

// Internal structure used to mark the ends of an Annotation
struct insertionMarker_st {
	Annotation *annot;
	boolean isStart;
};

Coordinate getInsertionMarkerPosition(InsertionMarker * marker)
{
	if (marker->isStart)
		return getStart(marker->annot);
	else
		return getFinish(marker->annot);
}

int compareInsertionMarkers(const void *A, const void *B)
{
	Coordinate Apos =
	    getInsertionMarkerPosition((InsertionMarker *) A);
	Coordinate Bpos =
	    getInsertionMarkerPosition((InsertionMarker *) B);

	if (Apos < Bpos)
		return -1;
	else if (Apos == Bpos)
		return 0;
	else
		return 1;
}

// Applies mergeSort to each insertion marker list (in order of position)
static void
orderInsertionMarkers(InsertionMarker ** insMarkers,
		      IDnum * markerCounters, RoadMapArray * rdmaps)
{
	IDnum sequenceIndex;
	IDnum sequenceCounter = rdmaps->length;

	puts("Ordering insertion markers");
	for (sequenceIndex = 1; sequenceIndex <= sequenceCounter;
	     sequenceIndex++) {
		qsort(insMarkers[sequenceIndex],
		      markerCounters[sequenceIndex],
		      sizeof(InsertionMarker), compareInsertionMarkers);
	}
}

// Creates insertion marker lists 
static void
setInsertionMarkers(RoadMapArray * rdmaps,
		    IDnum * markerCounters,
		    InsertionMarker ** veryLastMarker,
		    InsertionMarker ** insertionMarkers)
{
	IDnum sequenceCounter = rdmaps->length;
	IDnum sequenceIndex, sequenceIndex2;
	IDnum totalCount = 0;
	RoadMap *rdmap;
	Annotation *annot = rdmaps->annotations;
	InsertionMarker *nextMarker, *newMarker;
	IDnum annotIndex, lastAnnotIndex;
	InsertionMarker **insMarkers =
	    callocOrExit(rdmaps->length + 1, InsertionMarker *);
	// Counting insertion markers
	for (sequenceIndex = 1; sequenceIndex < sequenceCounter + 1;
	     sequenceIndex++) {
		//printf("Going through sequence %d\n", sequenceIndex);
		rdmap = getRoadMapInArray(rdmaps, sequenceIndex - 1);
		lastAnnotIndex = getAnnotationCount(rdmap);

		// Set insertion markers in previous sequences :

		for (annotIndex = 0; annotIndex < lastAnnotIndex;
		     annotIndex++) {
			if (getAnnotSequenceID(annot) > 0) {
				markerCounters[getAnnotSequenceID(annot)]
				    += 2;
			} else {
				markerCounters[-getAnnotSequenceID(annot)]
				    += 2;
			}
			totalCount += 2;
			annot = getNextAnnotation(annot);
		}
	}

	// Allocating space
	*insertionMarkers = callocOrExit(totalCount, InsertionMarker);
	*veryLastMarker = *insertionMarkers + totalCount;

	// Pointing each node to its space      
	nextMarker = *insertionMarkers;
	for (sequenceIndex = 1; sequenceIndex < sequenceCounter + 1;
	     sequenceIndex++) {
		insMarkers[sequenceIndex] = nextMarker;
		nextMarker = nextMarker + markerCounters[sequenceIndex];
		markerCounters[sequenceIndex] = 0;
	}

	// Filling up space with data
	annot = rdmaps->annotations;
	for (sequenceIndex = 1; sequenceIndex < sequenceCounter + 1;
	     sequenceIndex++) {
		//printf("Going through sequence %d\n", sequenceIndex);
		rdmap = getRoadMapInArray(rdmaps, sequenceIndex - 1);
		lastAnnotIndex = getAnnotationCount(rdmap);

		// Set insertion markers in previous sequences :

		for (annotIndex = 0; annotIndex < lastAnnotIndex;
		     annotIndex++) {
			sequenceIndex2 = getAnnotSequenceID(annot);
			if (sequenceIndex2 > 0) {
				newMarker =
				    insMarkers[sequenceIndex2] +
				    (markerCounters[sequenceIndex2])++;
				newMarker->annot = annot;
				newMarker->isStart = true;

				newMarker =
				    insMarkers[sequenceIndex2] +
				    (markerCounters[sequenceIndex2])++;
				newMarker->annot = annot;
				newMarker->isStart = false;
			} else {
				incrementAnnotationCoordinates(annot);

				newMarker =
				    insMarkers[-sequenceIndex2] +
				    (markerCounters[-sequenceIndex2])++;
				newMarker->annot = annot;
				newMarker->isStart = true;

				newMarker =
				    insMarkers[-sequenceIndex2] +
				    (markerCounters[-sequenceIndex2])++;
				newMarker->annot = annot;
				newMarker->isStart = false;
			}
			annot = getNextAnnotation(annot);
		}
	}

	orderInsertionMarkers(insMarkers, markerCounters, rdmaps);
	free(insMarkers);
}

// Counts how many preNodes are to be created to allocate appropriate memory
static void
countPreNodes(RoadMapArray * rdmaps, PreGraph * preGraph,
	      IDnum * markerCounters, InsertionMarker * insertionMarkers,
	      InsertionMarker * veryLastMarker)
{
	Annotation *annot = rdmaps->annotations;
	InsertionMarker *currentMarker = insertionMarkers;
	IDnum markerIndex, lastMarkerIndex;
	IDnum sequenceIndex;
	Coordinate currentPosition, nextStop;
	IDnum preNodeCounter = 0;
	RoadMap *rdmap;
	IDnum annotIndex, lastAnnotIndex;

	// Now that we have read all of the annotations, we go on to create the preNodes and tie them up
	for (sequenceIndex = 1;
	     sequenceIndex <= sequenceCount_pg(preGraph);
	     sequenceIndex++) {
		rdmap = getRoadMapInArray(rdmaps, sequenceIndex - 1);
		annotIndex = 0;
		lastAnnotIndex = getAnnotationCount(rdmap);
		markerIndex = 0;
		lastMarkerIndex = markerCounters[sequenceIndex];
		currentPosition = 0;


		while (annotIndex < lastAnnotIndex) {
			if (markerIndex == lastMarkerIndex
			    || getPosition(annot) <=
			    getInsertionMarkerPosition(currentMarker))
				nextStop = getPosition(annot);
			else
				nextStop =
				    getInsertionMarkerPosition
				    (currentMarker);

			if (currentPosition != nextStop) {
				preNodeCounter++;
				currentPosition = nextStop;
			}

			while (markerIndex < lastMarkerIndex
			       && getInsertionMarkerPosition(currentMarker)
			       == currentPosition) {
				currentMarker++;
				markerIndex++;
			}

			while (annotIndex < lastAnnotIndex
			       && getPosition(annot) == currentPosition) {
				annot = getNextAnnotation(annot);
				annotIndex++;
			}

		}

		while (markerIndex < lastMarkerIndex) {
			if (currentPosition ==
			    getInsertionMarkerPosition(currentMarker)) {
				currentMarker++;
				markerIndex++;
			} else {
				preNodeCounter++;
				currentPosition =
				    getInsertionMarkerPosition
				    (currentMarker);
			}
		}
	}

	allocatePreNodeSpace_pg(preGraph, preNodeCounter);
}

static void convertInsertionMarkers(InsertionMarker * insertionMarkers,
				    InsertionMarker * veryLastMarker,
				    IDnum * chains)
{
	InsertionMarker *marker;
	Annotation *annot;

	for (marker = insertionMarkers; marker != veryLastMarker; marker++) {
		annot = marker->annot;

		if (getAnnotSequenceID(annot) > 0) {
			if (marker->isStart) {
				if (getStartID(annot) == 0)
					setStartID(annot,
						   chains
						   [getAnnotSequenceID
						    (annot)]);
				else
					setStartID(annot,
						   getStartID(annot) + 1);
			}
		} else {
			if (marker->isStart)
				setStartID(annot, -getStartID(annot));
			else {
				if (getFinishID(annot) == 0)
					setFinishID(annot,
						    -chains
						    [-getAnnotSequenceID
						     (annot)]);
				else
					setFinishID(annot,
						    -getFinishID(annot) -
						    1);
			}
		}
	}

	free(insertionMarkers);
}

static void convertMarker(InsertionMarker * marker, IDnum nodeID)
{
	if (marker->isStart)
		setStartID(marker->annot, nodeID);
	else
		setFinishID(marker->annot, nodeID);
}

// Creates the preNode using insertion marker and annotation lists for each sequence
static void
createPreNodes(RoadMapArray * rdmaps, PreGraph * preGraph,
	       IDnum * markerCounters, InsertionMarker * insertionMarkers,
	       InsertionMarker * veryLastMarker, IDnum * chains,
	       char *sequenceFilename, int WORDLENGTH)
{
	Annotation *annot = rdmaps->annotations;
	IDnum latestPreNodeID;
	InsertionMarker *currentMarker = insertionMarkers;
	IDnum sequenceIndex;
	Coordinate currentPosition, nextStop;
	IDnum preNodeCounter = 1;
	FILE *file = fopen(sequenceFilename, "r");
	char line[50000];
	int lineLength = 50000;
	Coordinate readIndex;
	boolean tooShort;
	Kmer initialKmer;
	char c;
	RoadMap *rdmap;
	IDnum annotIndex, lastAnnotIndex;
	IDnum markerIndex, lastMarkerIndex;

	if (file == NULL) 
		exitErrorf(EXIT_FAILURE, true, "Could not read %s", sequenceFilename);
	// Reading sequence descriptor in first line
	if (!fgets(line, lineLength, file))
		exitErrorf(EXIT_FAILURE, true, "%s incomplete.", sequenceFilename);

	// Now that we have read all of the annotations, we go on to create the preNodes and tie them up
	for (sequenceIndex = 1;
	     sequenceIndex <= sequenceCount_pg(preGraph);
	     sequenceIndex++) {
		if (sequenceIndex % 100000 == 0)
			printf("Sequence %d / %d\n", sequenceIndex,
			       sequenceCount_pg(preGraph));

		while (line[0] != '>')
			if (!fgets(line, lineLength, file))
				exitErrorf(EXIT_FAILURE, true, "%s incomplete.", sequenceFilename);

		rdmap = getRoadMapInArray(rdmaps, sequenceIndex - 1);
		annotIndex = 0;
		lastAnnotIndex = getAnnotationCount(rdmap);
		markerIndex = 0;
		lastMarkerIndex = markerCounters[sequenceIndex];
		currentPosition = 0;

		// Reading first (k-1) nucleotides
		tooShort = false;
		clearKmer(&initialKmer);
		//printf("Initial kmer: ");
		for (readIndex = 0; readIndex < WORDLENGTH - 1;
		     readIndex++) {
			if (!isalpha((c = getc(file)))) {
				if (c == '>') {
					ungetc(c, file);
					tooShort = true;
					break;
				} else {
					continue;
				}
			}
			//printf("%c", c);      
			switch (c) {
			case 'A':
				pushNucleotide(&initialKmer, ADENINE);
				break;
			case 'C':
				pushNucleotide(&initialKmer, CYTOSINE);
				break;
			case 'G':
				pushNucleotide(&initialKmer, GUANINE);
				break;
			case 'T':
				pushNucleotide(&initialKmer, THYMINE);
				break;
			default:
				printf
				    ("Irregular sequence file: are you sure your Sequence and Roadmap file come from the same source?\n");
				fflush(stdout);
				abort();
			}
		}
		//puts("");

		if (tooShort) {
			//printf("Skipping short read.. %d\n", sequenceIndex);
			chains[sequenceIndex] = preNodeCounter;
			if (!fgets(line, lineLength, file))
				exitErrorf(EXIT_FAILURE, true, "%s incomplete.", sequenceFilename);
			continue;
		}

		latestPreNodeID = 0;

		while (annotIndex < lastAnnotIndex) {
			if (markerIndex == lastMarkerIndex
			    || getPosition(annot) <=
			    getInsertionMarkerPosition(currentMarker))
				nextStop = getPosition(annot);
			else {
				nextStop =
				    getInsertionMarkerPosition
				    (currentMarker);
			}

			if (currentPosition != nextStop) {
				addPreNodeToPreGraph_pg(preGraph,
							currentPosition,
							nextStop,
							file,
							&initialKmer,
							preNodeCounter);
				if (latestPreNodeID == 0) {
					chains[sequenceIndex] =
					    preNodeCounter;
				}
				latestPreNodeID = preNodeCounter++;
				currentPosition = nextStop;
			}

			while (markerIndex < lastMarkerIndex
			       && getInsertionMarkerPosition(currentMarker)
			       == nextStop) {
				convertMarker(currentMarker,
					      latestPreNodeID);
				currentMarker++;
				markerIndex++;
			}

			while (annotIndex < lastAnnotIndex
			       && getPosition(annot) == nextStop) {
				for (readIndex = 0;
				     readIndex <
				     getAnnotationLength(annot);
				     readIndex++) {
					c = getc(file);
					while (!isalpha(c))
						c = getc(file);

					//printf("(%c)", c);
					switch (c) {
					case 'A':
						pushNucleotide(&initialKmer, ADENINE);
						break;
					case 'C':
						pushNucleotide(&initialKmer, CYTOSINE);
						break;
					case 'G':
						pushNucleotide(&initialKmer, GUANINE);
						break;
					case 'T':
						pushNucleotide(&initialKmer, THYMINE);
						break;
					default:
						printf
						    ("Irregular sequence file: are you sure your Sequence and Roadmap file come from the same source?\n");
						fflush(stdout);
						exit(1);
					}
				}

				annot = getNextAnnotation(annot);
				annotIndex++;
			}

		}

		while (markerIndex < lastMarkerIndex) {
			if (currentPosition ==
			    getInsertionMarkerPosition(currentMarker)) {
				convertMarker(currentMarker,
					      latestPreNodeID);
				currentMarker++;
				markerIndex++;
			} else {
				nextStop =
				    getInsertionMarkerPosition
				    (currentMarker);
				addPreNodeToPreGraph_pg(preGraph,
							currentPosition,
							nextStop, file,
							&initialKmer,
							preNodeCounter);
				if (latestPreNodeID == 0)
					chains[sequenceIndex] =
					    preNodeCounter;
				latestPreNodeID = preNodeCounter++;
				currentPosition =
				    getInsertionMarkerPosition
				    (currentMarker);
			}
		}

		// End of sequence
		if (!fgets(line, lineLength, file) && sequenceIndex < sequenceCount_pg(preGraph))
			exitErrorf(EXIT_FAILURE, true, "%s incomplete.", sequenceFilename);
		//puts(" ");

		if (latestPreNodeID == 0)
			chains[sequenceIndex] = preNodeCounter;
	}

	free(markerCounters);
	fclose(file);

}

static void connectPreNodeToTheNext(IDnum * currentPreNodeID,
				    IDnum nextPreNodeID,
				    Coordinate * currentPosition,
				    PassageMarker ** latestPassageMarker,
				    IDnum sequenceIndex,
				    PreGraph * preGraph)
{
	if (nextPreNodeID == 0)
		return;

	if (*currentPreNodeID != 0)
		createPreArc_pg(*currentPreNodeID, nextPreNodeID,
				preGraph);

	*currentPreNodeID = nextPreNodeID;

	*currentPosition +=
	    getPreNodeLength_pg(*currentPreNodeID, preGraph);

}

static IDnum chooseNextInternalPreNode(IDnum currentPreNodeID,
				       IDnum sequenceIndex,
				       PreGraph * preGraph, IDnum * chains)
{
	if (currentPreNodeID >= preNodeCount_pg(preGraph))
		return 0;
	if (sequenceIndex >= sequenceCount_pg(preGraph))
		return currentPreNodeID + 1;
	if (currentPreNodeID + 1 < chains[sequenceIndex + 1])
		return currentPreNodeID + 1;
	return 0;
}

static void connectAnnotation(IDnum * currentPreNodeID, Annotation * annot,
			      Coordinate * currentPosition,
			      PassageMarker ** latestPassageMarker,
			      IDnum sequenceIndex, PreGraph * preGraph)
{
	IDnum nextPreNodeID = getStartID(annot);

	connectPreNodeToTheNext(currentPreNodeID, nextPreNodeID,
				currentPosition, latestPassageMarker,
				sequenceIndex, preGraph);

	while (*currentPreNodeID != getFinishID(annot)) {
		nextPreNodeID = (*currentPreNodeID) + 1;

		connectPreNodeToTheNext(currentPreNodeID, nextPreNodeID,
					currentPosition,
					latestPassageMarker, sequenceIndex,
					preGraph);
	}
}

// Threads each sequences and creates preArcs according to road map indications
static void connectPreNodes(RoadMapArray * rdmaps, PreGraph * preGraph,
			    IDnum * chains)
{
	Coordinate currentPosition, currentInternalPosition;
	IDnum sequenceIndex;
	Annotation *annot = rdmaps->annotations;
	IDnum currentPreNodeID, nextInternalPreNodeID;
	PassageMarker *latestPassageMarker;
	RoadMap *rdmap;
	IDnum annotIndex, lastAnnotIndex;

	for (sequenceIndex = 1;
	     sequenceIndex <= sequenceCount_pg(preGraph);
	     sequenceIndex++) {

		if (sequenceIndex % 100000 == 0)
			printf("Connecting %d / %d\n", sequenceIndex,
			       sequenceCount_pg(preGraph));

		rdmap = getRoadMapInArray(rdmaps, sequenceIndex - 1);
		annotIndex = 0;
		lastAnnotIndex = getAnnotationCount(rdmap);
		nextInternalPreNodeID = chooseNextInternalPreNode
		    (chains[sequenceIndex] - 1, sequenceIndex,
		     preGraph, chains);

		currentPosition = 0;
		currentInternalPosition = 0;
		currentPreNodeID = 0;
		latestPassageMarker = NULL;
		// Recursion up to last annotation
		while (annotIndex < lastAnnotIndex
		       || nextInternalPreNodeID != 0) {
			if (annotIndex == lastAnnotIndex
			    || (nextInternalPreNodeID != 0
				&& currentInternalPosition <
				getPosition(annot))) {
				connectPreNodeToTheNext(&currentPreNodeID,
							nextInternalPreNodeID,
							&currentPosition,
							&latestPassageMarker,
							sequenceIndex,
							preGraph);
				nextInternalPreNodeID =
				    chooseNextInternalPreNode
				    (currentPreNodeID, sequenceIndex,
				     preGraph, chains);
				currentInternalPosition +=
				    getPreNodeLength_pg(currentPreNodeID,
							preGraph);

			} else {
				connectAnnotation(&currentPreNodeID, annot,
						  &currentPosition,
						  &latestPassageMarker,
						  sequenceIndex, preGraph);
				annot = getNextAnnotation(annot);
				annotIndex++;
			}
		}
	}
}

// Post construction memory deallocation routine (of sorts, could certainly be optimized)
static void
cleanUpMemory(PreGraph * preGraph, RoadMapArray * rdmaps, IDnum * chains)
{
	// Killing off roadmaps
	free(rdmaps->annotations);
	free(rdmaps->array);
	free(rdmaps);

	// Finishing off the chain markers
	free(chains);
}

// The full monty, wrapped up in one function
PreGraph *newPreGraph_pg(RoadMapArray * rdmapArray, char *sequenceFilename)
{
	int WORDLENGTH = rdmapArray->WORDLENGTH;
	IDnum sequenceCount = rdmapArray->length;
	IDnum *markerCounters = callocOrExit(sequenceCount + 1, IDnum);
	IDnum *chains = callocOrExit(sequenceCount + 1, IDnum);
	InsertionMarker *insertionMarkers;
	InsertionMarker *veryLastMarker;

	PreGraph *preGraph =
	    emptyPreGraph_pg(sequenceCount, rdmapArray->WORDLENGTH, rdmapArray->double_strand);

	puts("Creating insertion markers");
	setInsertionMarkers(rdmapArray, markerCounters, &veryLastMarker,
			    &insertionMarkers);

	puts("Counting preNodes");
	countPreNodes(rdmapArray, preGraph, markerCounters,
		      insertionMarkers, veryLastMarker);

	printf("%d preNodes counted, creating them now\n",
	       preNodeCount_pg(preGraph));
	createPreNodes(rdmapArray, preGraph, markerCounters,
		       insertionMarkers, veryLastMarker, chains,
		       sequenceFilename, WORDLENGTH);

	puts("Adjusting marker info...");
	convertInsertionMarkers(insertionMarkers, veryLastMarker, chains);

	puts("Connecting preNodes");
	connectPreNodes(rdmapArray, preGraph, chains);

	puts("Cleaning up memory");
	cleanUpMemory(preGraph, rdmapArray, chains);
	puts("Concatenating preGraph");
	concatenatePreGraph_pg(preGraph);
	puts("Done creating preGraph");

	return preGraph;
}
