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
#include <limits.h>

#include "globals.h"
#include "graph.h"
#include "passageMarker.h"
#include "readSet.h"
#include "tightString.h"
#include "recycleBin.h"
#include "utility.h"
#include "kmer.h"

#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

typedef struct kmerOccurence_st KmerOccurence;
typedef struct kmerOccurenceTable_st KmerOccurenceTable;
typedef struct smallNodeList_st SmallNodeList;

// Internal structure used to mark the ends of an Annotation
struct kmerOccurence_st {
	Kmer kmer;
	Coordinate position;
	IDnum nodeID;
};

struct kmerOccurenceTable_st {
	KmerOccurence *kmerTable;
	IDnum *accelerationTable;
	IDnum kmerTableSize;
	short int accelerationShift;
	short int accelerationBits;
};

struct smallNodeList_st {
	Node *node;
	SmallNodeList *next;
};

static RecycleBin *smallNodeListMemory = NULL;
static SmallNodeList *nodePile = NULL;

#define BLOCKSIZE 1000

static SmallNodeList *allocateSmallNodeList()
{
	if (smallNodeListMemory == NULL)
		smallNodeListMemory =
		    newRecycleBin(sizeof(SmallNodeList), BLOCKSIZE);

	return allocatePointer(smallNodeListMemory);
}

static void deallocateSmallNodeList(SmallNodeList * smallNodeList)
{
	deallocatePointer(smallNodeListMemory, smallNodeList);
}

static void memorizeNode(Node * node)
{
	SmallNodeList *list = allocateSmallNodeList();
	list->node = node;
	list->next = nodePile;
	nodePile = list;
}

static void unlockMemorizedNodes()
{
	SmallNodeList *list;

	while (nodePile) {
		list = nodePile;
		nodePile = list->next;
		setSingleNodeStatus(list->node, false);
		deallocateSmallNodeList(list);
	}
}

int compareKmerOccurences(void const *A, void const *B)
{
	KmerOccurence *a = (KmerOccurence *) A;
	KmerOccurence *b = (KmerOccurence *) B;

	if (compareKmers(&(a->kmer), &(b->kmer)) < 0)
		return -1;
	else if (compareKmers(&(a->kmer), &(b->kmer)) > 0)
		return 1;
	else
		return 0;
}

static inline KmerKey keyInAccelerationTable(Kmer * kmer,
					  KmerOccurenceTable * table)
{
	KmerKey key = 0;
	Kmer copy;
	int i;

	copyKmers(&copy, kmer);
	for (i = 0; i < table->accelerationShift; i+= 2)
		popNucleotide(&copy);

	for (i = 0; i < table->accelerationBits; i+= 2) {
		key += ((KmerKey) popNucleotide(&copy)) << table->accelerationBits;
		key >>= 2;
	}
	
	return key;
}

static KmerOccurenceTable *referenceGraphKmers(char *preGraphFilename,
					       short int accelerationBits, Graph * graph, boolean double_strand)
{
	FILE *file = fopen(preGraphFilename, "r");
	const int maxline = MAXLINE;
	char line[MAXLINE];
	char c;
	int wordLength;
	Coordinate lineLength, kmerCount;
	Kmer word;
	Kmer antiWord;
	KmerOccurenceTable *kmerTable = NULL;
	KmerOccurence *kmerOccurences, *kmerOccurencePtr;
	Coordinate kmerOccurenceIndex;
	IDnum index;
	IDnum nodeID = 0;
	IDnum *accelPtr = NULL;
	KmerKey lastHeader = 0;
	KmerKey header;
	Nucleotide nucleotide;

	if (file == NULL)
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", preGraphFilename);

	// Count kmers
	printf("Scanning pre-graph file %s for k-mers\n",
	       preGraphFilename);

	// First  line
	if (!fgets(line, maxline, file))
		exitErrorf(EXIT_FAILURE, true, "PreGraph file incomplete");
	sscanf(line, "%*i\t%*i\t%i\n", &wordLength);

	// Initialize kmer occurence table:
	kmerTable = mallocOrExit(1, KmerOccurenceTable);
	if (accelerationBits > 2 * wordLength)
		accelerationBits = 2 * wordLength;

	if (accelerationBits > 32)
		accelerationBits = 32;

	if (accelerationBits > 0) {
		kmerTable->accelerationBits = accelerationBits;
		kmerTable->accelerationTable =
		    callocOrExit((((size_t) 1) << accelerationBits) + 1,
			   IDnum);
		accelPtr = kmerTable->accelerationTable;
		kmerTable->accelerationShift =
		    (short int) 2 *wordLength - accelerationBits;
	} else {
		kmerTable->accelerationBits = 0;
		kmerTable->accelerationTable = NULL;
		kmerTable->accelerationShift = 0;
	}

	// Read nodes
	if (!fgets(line, maxline, file))
		exitErrorf(EXIT_FAILURE, true, "PreGraph file incomplete");
	kmerCount = 0;
	while (line[0] == 'N') {
		lineLength = 0;
		while ((c = getc(file)) != EOF && c != '\n')
			lineLength++;
		kmerCount += lineLength - wordLength + 1;
		if (fgets(line, maxline, file) == NULL)
			break;
	}
	fclose(file);

	// Create table
	printf("%li kmers found\n", (long) kmerCount);
	kmerOccurences = callocOrExit(kmerCount, KmerOccurence);
	kmerOccurencePtr = kmerOccurences;
	kmerOccurenceIndex = 0;
	kmerTable->kmerTable = kmerOccurences;
	kmerTable->kmerTableSize = kmerCount;

	// Fill table
	file = fopen(preGraphFilename, "r");
	if (file == NULL)
		exitErrorf(EXIT_FAILURE, true, "Could not open %s", preGraphFilename);

	if (!fgets(line, maxline, file))
		exitErrorf(EXIT_FAILURE, true, "PreGraph file incomplete");

	// Read nodes
	if (!fgets(line, maxline, file))
		exitErrorf(EXIT_FAILURE, true, "PreGraph file incomplete");
	while (line[0] == 'N') {
		nodeID++;

		// Fill in the initial word : 
		clearKmer(&word);
		clearKmer(&antiWord);

		for (index = 0; index < wordLength - 1; index++) {
			c = getc(file);
			if (c == 'A')
				nucleotide = ADENINE;
			else if (c == 'C')
				nucleotide = CYTOSINE;
			else if (c == 'G')
				nucleotide = GUANINE;
			else if (c == 'T')
				nucleotide = THYMINE;
			else if (c == '\n')
				exitErrorf(EXIT_FAILURE, true, "PreGraph file incomplete");
			else
				nucleotide = ADENINE;
				

			pushNucleotide(&word, nucleotide);
			if (double_strand) {
#ifdef COLOR
				reversePushNucleotide(&antiWord, nucleotide);
#else
				reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
			}
		}

		// Scan through node
		index = 0;
		while((c = getc(file)) != '\n' && c != EOF) {
			if (c == 'A')
				nucleotide = ADENINE;
			else if (c == 'C')
				nucleotide = CYTOSINE;
			else if (c == 'G')
				nucleotide = GUANINE;
			else if (c == 'T')
				nucleotide = THYMINE;
			else
				nucleotide = ADENINE;

			pushNucleotide(&word, nucleotide);
			if (double_strand) {
#ifdef COLOR
				reversePushNucleotide(&antiWord, nucleotide);
#else
				reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
			}

			if (!double_strand || compareKmers(&word, &antiWord) <= 0) {
				copyKmers(&kmerOccurencePtr->kmer, &word);
				kmerOccurencePtr->nodeID = nodeID;
				kmerOccurencePtr->position =
				    index;
			} else {
				copyKmers(&kmerOccurencePtr->kmer, &antiWord);
				kmerOccurencePtr->nodeID = -nodeID;
				kmerOccurencePtr->position =
				    getNodeLength(getNodeInGraph(graph, nodeID)) - 1 - index;
			}

			kmerOccurencePtr++;
			kmerOccurenceIndex++;
			index++;
		}

		if (fgets(line, maxline, file) == NULL)
			break;
	}

	fclose(file);

	// Sort table
	qsort(kmerOccurences, kmerCount, sizeof(KmerOccurence),
	      compareKmerOccurences);

	// Fill up acceleration table
	if (kmerTable->accelerationTable != NULL) {
		*accelPtr = (IDnum) 0;
		for (kmerOccurenceIndex = 0;
		     kmerOccurenceIndex < kmerCount;
		     kmerOccurenceIndex++) {
			header =
			    keyInAccelerationTable(&kmerOccurences
						   [kmerOccurenceIndex].
						   kmer, kmerTable);
			while (lastHeader < header) {
				lastHeader++;
				accelPtr++;
				*accelPtr = kmerOccurenceIndex;
			}
		}

		while (lastHeader < (KmerKey) 1 << accelerationBits) {
			lastHeader++;
			accelPtr++;
			*accelPtr = kmerCount;
		}
	}

	return kmerTable;
}

static KmerOccurence *findKmerOccurenceInSortedTable(Kmer * kmer,
						     KmerOccurenceTable *
						     table)
{
	KmerOccurence *array = table->kmerTable;
	KmerKey key = keyInAccelerationTable(kmer, table);
	Coordinate leftIndex, rightIndex, middleIndex;

	if (table->accelerationTable != NULL) {
		leftIndex = table->accelerationTable[key];
		rightIndex = table->accelerationTable[key + 1];
	} else {
		leftIndex = 0;
		rightIndex = table->kmerTableSize;
	}

	while (true) {
		middleIndex = (rightIndex + leftIndex) / 2;

		if (leftIndex >= rightIndex)
			return NULL;
		else if (compareKmers(&(array[middleIndex]).kmer, kmer) == 0) 
			return &(array[middleIndex]);
		else if (leftIndex == middleIndex)
			return NULL;
		else if (compareKmers(&(array[middleIndex]).kmer, kmer) > 0)
			rightIndex = middleIndex;
		else
			leftIndex = middleIndex;
	}
}

static void ghostThreadSequenceThroughGraph(TightString * tString,
					    KmerOccurenceTable *
					    kmerOccurences, Graph * graph,
					    IDnum seqID, Category category,
					    boolean readTracking,
					    boolean double_strand,
					    boolean second_in_pair)
{
	Kmer word;
	Kmer antiWord;
	Coordinate readNucleotideIndex;
	KmerOccurence *kmerOccurence;
	int wordLength = getWordLength(graph);
	Nucleotide nucleotide;
	boolean reversed;

	Node *node;
	Node *previousNode = NULL;

	clearKmer(&word);
	clearKmer(&antiWord);

	// Neglect any read which will not be short paired
	if ((!readTracking && category % 2 == 0)
	    || category / 2 >= CATEGORIES)
		return;

	// Neglect any string shorter than WORDLENGTH :
	if (getLength(tString) < wordLength)
		return;

	// Verify that all short reads are reasonnably short
	if (getLength(tString) > USHRT_MAX) {
		printf("Short read of length %lli, longer than limit %i\n",
		       (long long) getLength(tString), SHRT_MAX);
		puts("You should better declare this sequence as long, because it genuinely is!");
		exit(1);
	}
	// Allocate memory for the read pairs
	if (!readStartsAreActivated(graph))
		activateReadStarts(graph);

	// Fill in the initial word : 
	for (readNucleotideIndex = 0;
	     readNucleotideIndex < wordLength - 1; readNucleotideIndex++) {
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);
		if (double_strand) {
#ifdef COLOR
			reversePushNucleotide(&antiWord, nucleotide);
#else
			reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
		}
	}

	// Go through sequence
	while (readNucleotideIndex < getLength(tString)) {
		// Shift word:
		nucleotide = getNucleotide(readNucleotideIndex++, tString);
		pushNucleotide(&word, nucleotide);
		if (double_strand) {
#ifdef COLOR
			reversePushNucleotide(&antiWord, nucleotide);
#else
			reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
		}

		// Search in table
		reversed = false;
		if (double_strand) {
			if (compareKmers(&word, &antiWord) <= 0) {
			    	kmerOccurence =
				findKmerOccurenceInSortedTable(&word,
							       kmerOccurences);
			} else { 
				kmerOccurence =
				       findKmerOccurenceInSortedTable(&antiWord,
				        kmerOccurences);
				reversed = true;
			}
		} else {
			if (!second_in_pair) {
			    	kmerOccurence =
				findKmerOccurenceInSortedTable(&word,
							       kmerOccurences);
			} else { 
				kmerOccurence =
				       findKmerOccurenceInSortedTable(&antiWord,
				        kmerOccurences);
				reversed = true;
			}
		}
		
		if (kmerOccurence) {
			if (!reversed) 
				node = getNodeInGraph(graph, kmerOccurence->nodeID);
			else
				node = getNodeInGraph(graph, -kmerOccurence->nodeID);
		} else {
			node = NULL;
			if (previousNode) 
				break;
		}

		previousNode = node;

		// Fill in graph
		if (node && !getNodeStatus(node)) {
			incrementReadStartCount(node, graph);
			setSingleNodeStatus(node, true);
			memorizeNode(node);
		}
	}

	unlockMemorizedNodes();
}

static void threadSequenceThroughGraph(TightString * tString,
				       KmerOccurenceTable * kmerOccurences,
				       Graph * graph,
				       IDnum seqID, Category category,
				       boolean readTracking,
				       boolean double_strand, 
				       boolean second_in_pair)
{
	Kmer word;
	Kmer antiWord;
	Coordinate readNucleotideIndex;
	Coordinate kmerIndex;
	KmerOccurence *kmerOccurence;
	int wordLength = getWordLength(graph);

	PassageMarker *marker = NULL;
	PassageMarker *previousMarker = NULL;
	Node *node;
	Node *previousNode = NULL;
	Coordinate coord;
	Coordinate previousCoord = 0;
	Nucleotide nucleotide;
	boolean reversed;

	clearKmer(&word);
	clearKmer(&antiWord);

	// Neglect any string shorter than WORDLENGTH :
	if (getLength(tString) < wordLength)
		return;

	// Fill in the initial word : 
	for (readNucleotideIndex = 0;
	     readNucleotideIndex < wordLength - 1; readNucleotideIndex++) {
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);
		if (double_strand) {
#ifdef COLOR
			reversePushNucleotide(&antiWord, nucleotide);
#else
			reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
		}
	}

	// Go through sequence
	while (readNucleotideIndex < getLength(tString)) {
		nucleotide = getNucleotide(readNucleotideIndex++, tString);
		pushNucleotide(&word, nucleotide);
		if (double_strand) {
#ifdef COLOR
			reversePushNucleotide(&antiWord, nucleotide);
#else
			reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
		}

		// Search in table
		reversed = false;
		if (double_strand) {
			if (compareKmers(&word, &antiWord) <= 0) {
			    	kmerOccurence =
				findKmerOccurenceInSortedTable(&word,
							       kmerOccurences);
			} else { 
				kmerOccurence =
				       findKmerOccurenceInSortedTable(&antiWord,
				        kmerOccurences);
				reversed = true;
			}
		} else {
			if (!second_in_pair) {
			    	kmerOccurence =
				findKmerOccurenceInSortedTable(&word,
							       kmerOccurences);
			} else { 
				kmerOccurence =
				       findKmerOccurenceInSortedTable(&antiWord,
				        kmerOccurences);
				reversed = true;
			}
		}
		
		if (kmerOccurence) {
			if (!reversed) {
				node = getNodeInGraph(graph, kmerOccurence->nodeID);
				coord = kmerOccurence->position;
			} else {
				node = getNodeInGraph(graph, -kmerOccurence->nodeID);
				coord = getNodeLength(node) - kmerOccurence->position - 1;
			}
		} else {
			node = NULL;
			if (previousNode) 
				break;
		}

		// Fill in graph
		if (node) {
			kmerIndex = readNucleotideIndex - wordLength;

			if (previousNode == node
			    && previousCoord == coord - 1) {
				if (category / 2 >= CATEGORIES) {
					setPassageMarkerFinish(marker,
							       kmerIndex +
							       1);
					setFinishOffset(marker,
							getNodeLength(node)
							- coord - 1);
				} else {
					incrementVirtualCoverage(node,
								 category /
								 2, 1);
					incrementOriginalVirtualCoverage
					    (node, category / 2, 1);
				}

			} else {
				if (category / 2 >= CATEGORIES) {
					marker =
					    newPassageMarker(seqID,
							     kmerIndex,
							     kmerIndex + 1,
							     coord,
							     getNodeLength
							     (node) -
							     coord - 1);
					transposePassageMarker(marker,
							       node);
					connectPassageMarkers
					    (previousMarker, marker,
					     graph);
					previousMarker = marker;
				} else {
					if (readTracking) {
						if (!getNodeStatus(node)) {
							addReadStart(node,
								     seqID,
								     coord,
								     graph,
								     kmerIndex);
							setSingleNodeStatus
							    (node, true);
							memorizeNode(node);
						} else {
							blurLastShortReadMarker
							    (node, graph);
						}
					}

					incrementVirtualCoverage(node,
								 category /
								 2, 1);
					incrementOriginalVirtualCoverage
					    (node, category / 2, 1);
				}

				createArc(previousNode, node, graph);
			}

			previousNode = node;
			previousCoord = coord;
		}
	}

	unlockMemorizedNodes();
}

static void fillUpGraph(ReadSet * reads,
			KmerOccurenceTable * kmerOccurences, Graph * graph,
			boolean readTracking, boolean double_strand)
{
	IDnum readIndex;
	Category category;
	boolean second_in_pair = false;

	resetNodeStatus(graph);

	for (readIndex = 0; readIndex < reads->readCount; readIndex++) {
		category = reads->categories[readIndex];
		ghostThreadSequenceThroughGraph(reads->
						tSequences[readIndex],
						kmerOccurences,
						graph, readIndex + 1,
						category, 
						readTracking, double_strand, second_in_pair);

		if (category % 2) 
			second_in_pair = (second_in_pair? false : true);
		else 
			second_in_pair = false;
	}

	createNodeReadStartArrays(graph);

	second_in_pair = false;
	for (readIndex = 0; readIndex < reads->readCount; readIndex++) {
		category = reads->categories[readIndex];

		if (readIndex % 100000 == 0)
			printf("Threading through reads %d / %d\n",
			       readIndex, reads->readCount);

		threadSequenceThroughGraph(reads->tSequences[readIndex],
					   kmerOccurences,
					   graph, readIndex + 1, category,
					   readTracking, double_strand, second_in_pair);

		if (category % 2) 
			second_in_pair = (second_in_pair? false : true);
		else 
			second_in_pair = false;
	}

	orderNodeReadStartArrays(graph);

	if (smallNodeListMemory != NULL)
		destroyRecycleBin(smallNodeListMemory);

	free(kmerOccurences->kmerTable);
	free(kmerOccurences->accelerationTable);
	free(kmerOccurences);
}

Graph *importPreGraph(char *preGraphFilename, ReadSet * reads,
		      boolean readTracking, short int accelerationBits)
{
	boolean double_strand = false;
	Graph *graph = readPreGraphFile(preGraphFilename, &double_strand);
	KmerOccurenceTable *kmerOccurences =
	    referenceGraphKmers(preGraphFilename, accelerationBits, graph, double_strand);
	fillUpGraph(reads, kmerOccurences, graph, readTracking, double_strand);

	return graph;
}
