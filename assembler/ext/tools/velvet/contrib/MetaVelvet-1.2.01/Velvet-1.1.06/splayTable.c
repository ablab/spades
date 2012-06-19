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
#include <time.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "globals.h"
#include "readSet.h"
#include "splay.h"
#include "tightString.h"
#include "utility.h"
#include "kmer.h"
#include "kmerOccurenceTable.h"
#include "recycleBin.h"

typedef struct mask_st Mask;

struct mask_st {
	Coordinate start;
	Coordinate finish;
	Mask* next;
}; 

static RecycleBin * maskMemory = NULL;

static Mask *allocateMask()
{
	if (maskMemory == NULL)
		maskMemory = newRecycleBin(sizeof(Mask), 10000);

	return (Mask *) allocatePointer(maskMemory);
}

static Mask * newMask(Coordinate position)
{
	Mask * mask = allocateMask();
	mask->start = position;
	mask->finish = position;
	mask->next = NULL;
	return mask;
}

#define HASH_BUCKETS_NB 16777216

#ifdef _OPENMP

#define NB_PUSH 32
#define BUFFER_SIZE 4096

static StringBuffer **annotationBuffer = NULL;
static StringBuffer **annotationBufferW = NULL;
static int *nbPush = NULL;
static boolean producing = 1;

static void initAnnotationBuffers(void)
{
	int n;
	int i;

	n = omp_get_max_threads();
	annotationBuffer = callocOrExit(n, StringBuffer*);
	annotationBufferW = callocOrExit(n, StringBuffer*);
	nbPush = callocOrExit(n, int);

	for (i = 0; i < n; i++)
	{
		annotationBuffer[i] = newStringBuffer(BUFFER_SIZE);
		annotationBufferW[i] = newStringBuffer(BUFFER_SIZE);
	}
}

static void destroyAnnotationBuffers(void)
{
	int n;
	int i;

	n = omp_get_max_threads();

	for (i = 0; i < n; i++)
	{
		destroyStringBuffer(annotationBuffer[i], 1);
		destroyStringBuffer(annotationBufferW[i], 1);
	}

	free(annotationBuffer);
	free(annotationBufferW);
	free(nbPush);
	annotationBuffer = NULL;
	annotationBufferW = NULL;
	nbPush = NULL;
}

static void pushBufferCommit(int thread)
{
	StringBuffer *tmp;
	char *s;

	s = annotationBufferW[thread]->str;
	do
	{
		#pragma omp flush(s)
	}
	while (*s);
	tmp = annotationBufferW[thread];
	annotationBufferW[thread] = annotationBuffer[thread];
	annotationBuffer[thread] = tmp;
	tmp = annotationBufferW[thread];
	#pragma omp flush(tmp)
}

static void pushBuffer(int thread)
{
	if (++nbPush[thread] == NB_PUSH)
	{
		nbPush[thread] = 0;
		pushBufferCommit(thread);
	}
}

static void writeBuffers(FILE *outFile, int nbThreads)
{
	int i;

	for (i = 0; i < nbThreads; i++)
	{
		StringBuffer *b;
		char *s;

		b = annotationBufferW[i];
		#pragma omp flush(b)
		s = b->str;
		#pragma omp flush(s)
		if (*s)
		{
			velvetFprintf(outFile, "%s", annotationBufferW[i]->str);
			resetStringBuffer(annotationBufferW[i]);
		}
	}
}

static void bufferWritter(FILE *outFile)
{
	int n;

	n = omp_get_max_threads();
	#pragma omp flush(producing)
	while (producing)
	{
		writeBuffers(outFile, n);
		#pragma omp flush(producing)
	}
	writeBuffers(outFile, n);
}

static void appendLine(char *line, int thread)
{
	appendStringBuffer(annotationBuffer[thread], line);
}
#else

#define BUFFER_SIZE 1024

StringBuffer *annotationBuffer = NULL;

static void appendLine(char *line, int thread)
{
	appendStringBuffer(annotationBuffer, line);
}
#endif

struct splayTable_st {
	SplayTree **table;
#ifdef _OPENMP
	omp_lock_t *tableLocks;
#endif
	KmerOccurenceTable *kmerOccurenceTable;
	int WORDLENGTH;
	boolean double_strand;
};

SplayTable *newSplayTable(int WORDLENGTH, boolean double_strand)
{
	SplayTable *splayTable = mallocOrExit(1, SplayTable);
	splayTable->WORDLENGTH = WORDLENGTH;
	splayTable->table = callocOrExit(HASH_BUCKETS_NB, SplayTree *);
	splayTable->kmerOccurenceTable = NULL;
	splayTable->double_strand = double_strand;
#ifdef _OPENMP
	splayTable->tableLocks = mallocOrExit(HASH_BUCKETS_NB, omp_lock_t);
	int i;
	#pragma omp parallel for
	for (i = 0; i < HASH_BUCKETS_NB; i++)
		omp_init_lock(splayTable->tableLocks + i);
	initSplayTreeMemory();
#endif
	return splayTable;
}

void destroySplayTable(SplayTable * splayTable)
{
	velvetLog("Destroying splay table\n");

	destroyAllSplayTrees();
	free(splayTable->table);
	destroyKmerOccurenceTable(splayTable->kmerOccurenceTable);
	free(splayTable);

	velvetLog("Splay table destroyed\n");
}

static KmerKey hash_kmer(Kmer * kmer)
{
#if KMER_LONGLONGS
	KmerKey key = kmer->longlongs[0];

#if KMER_LONGLONGS > 1
	key ^= kmer->longlongs[1];
#endif
#if KMER_LONGLONGS > 2
	key ^= kmer->longlongs[2];
#endif

	key = (~key) + (key << 21);
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8);
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4);
	key = key ^ (key >> 28);
	key = key + (key << 31);

	return key % HASH_BUCKETS_NB;
#elif KMER_LONGS
	KmerKey key = kmer->longs;

	key += ~(key << 15);
	key ^= (key >> 10);
	key += (key << 3);
	key ^= (key >> 6);
	key += ~(key << 11);
	key ^= (key >> 16);

	return key % HASH_BUCKETS_NB;

#elif KMER_INTS
	return kmer->ints % HASH_BUCKETS_NB;
#elif KMER_CHARS
	return kmer->chars % HASH_BUCKETS_NB;
#endif
}

static Coordinate getNearestHSPIndex(Coordinate position, IDnum * sequenceIDs, Coordinate sequenceLength) {
	Coordinate back_offset = -1;
	Coordinate front_offset = -1;

	for (back_offset = 1; position - back_offset > 0; back_offset++) 
		if (sequenceIDs[position - back_offset])
			break;

	for (front_offset = 1; position + front_offset < sequenceLength; front_offset++) 
		if (sequenceIDs[position + front_offset])
			break;

	if (back_offset == position && position + front_offset == sequenceLength) 
		return -1;
	else if (back_offset == position)
		return position + front_offset;
	else if (front_offset + position == sequenceLength)
		return position - back_offset;
	else
		return back_offset < front_offset? position - back_offset : position + front_offset;
}

static KmerOccurence * getMostAppropriateHit(Coordinate readCoord, Coordinate readLength, boolean direct, KmerOccurence * kmerOccurence, IDnum mapCount, IDnum * mapSequenceID, Coordinate * mapCoord, int wordLength) {
	KmerOccurence * current;
	KmerOccurence * best = NULL;
	Coordinate expectedPosition;
	Coordinate positionError;
	IDnum mapIndex;

	// If only one hit
	if (!getNextKmerOccurence(kmerOccurence))
		return kmerOccurence;

	// If multiple hits by unmapped read
	if (mapCount == 0)
		return NULL;

	// Compare cases
	for (current = kmerOccurence; current; current = getNextKmerOccurence(current)) {
		for (mapIndex = 0; mapIndex < mapCount; mapIndex++) {

			// If wrong sequence or unconsistent orientation	
			if ((direct && getKmerOccurenceNodeID(current) != mapSequenceID[mapIndex]) 
			    || (!direct && getKmerOccurenceNodeID(current) != -mapSequenceID[mapIndex])) 
				continue;

			// Compute where it is supposed to land on reference
			if (mapSequenceID[mapIndex] < 0)
				expectedPosition = mapCoord[mapIndex] + readLength - readCoord - 1;
			else 
				expectedPosition = mapCoord[mapIndex] + readCoord - wordLength + 1;
		
			// Compute positional error
			positionError = getKmerOccurencePosition(current) - expectedPosition;

			// If potential hit record
			if (positionError < 1 && positionError > -1) {
				if (best)
					// If competing hit, give up
					return NULL;
				else
					// Record current hit
					best = current;
			}
		}
	}

	return best;
}

static inline boolean
doFindOrInsertOccurenceInSplayTree(Kmer * kmer, IDnum * seqID,
				   Coordinate * position, SplayTable *table)
{
#ifdef _OPENMP
	const KmerKey kmerHash = hash_kmer(kmer);
	boolean ret;

	omp_set_lock(table->tableLocks + kmerHash);
	ret =  findOrInsertOccurenceInSplayTree(kmer, seqID, position,
						table->table + kmerHash);
	omp_unset_lock(table->tableLocks + kmerHash);

	return ret;
#else
	return findOrInsertOccurenceInSplayTree(kmer, seqID, position,
						&table->table[hash_kmer(kmer)]);
#endif
}


static boolean findOrInsertOccurenceInSplayTable(Kmer * kmer, IDnum * seqID,
					      Coordinate * position,
					      SplayTable * table, IDnum * sequenceIDs,
						 Coordinate * coords, Coordinate readIndex, Coordinate readLength, boolean direct)
{
	KmerOccurence * hit;
	Coordinate HSPIndex;

	// Check if previous anchor
	if (sequenceIDs && sequenceIDs[readIndex]) {
		if (direct)
			*seqID = sequenceIDs[readIndex];
		else 
			*seqID = -sequenceIDs[readIndex];
		if (sequenceIDs[readIndex] > 0) 
			*position = coords[readIndex] + readIndex;
		else
			*position = coords[readIndex] - readIndex + readLength - 1;

		return true;
	}
	else if (coords && coords[readIndex]) 
		// If in buffer zone:
		return doFindOrInsertOccurenceInSplayTree(kmer, seqID, position, table);

	// Look up first in reference sequence k-mers
	if (table->kmerOccurenceTable 
	    && (hit = findKmerInKmerOccurenceTable(kmer, table->kmerOccurenceTable))) {
		if (!getNextKmerOccurence(hit)) {
			*seqID = getKmerOccurenceNodeID(hit);
			*position = getKmerOccurencePosition(hit);
			return true;
		} else if ((HSPIndex = getNearestHSPIndex(*position, sequenceIDs, readLength)) > 0) {
	    		hit = getMostAppropriateHit(readIndex, readLength, direct, hit, 1, &(sequenceIDs[HSPIndex]), &(coords[HSPIndex]), table->WORDLENGTH);
			if (hit) {
				*seqID = getKmerOccurenceNodeID(hit);
				*position = getKmerOccurencePosition(hit);
				return true;
			}

		}
	} 

	// If not, go through the novel k-mers
	return doFindOrInsertOccurenceInSplayTree(kmer, seqID, position, table);
}

static void printAnnotations(IDnum *sequenceIDs, Coordinate * coords,
			     TightString * array, SplayTable * table,
			     FILE * file, boolean second_in_pair, IDnum seqID) 
{
	Coordinate readNucleotideIndex = 0;
	Coordinate writeNucleotideIndex = 0;
	Kmer word;
	Kmer antiWord;
	boolean annotationClosed = true;
	IDnum sequenceID;
	Coordinate coord;
	boolean found;
	Coordinate position = 0;
	Coordinate start = 0;
	Coordinate finish = 0;
	IDnum referenceSequenceID = 0;
	Nucleotide nucleotide;
	char lineBuffer[MAXLINE];
	TightString * tString = getTightStringInArray(array, seqID - 1); 
	int thread = 0;

	clearKmer(&word);
	clearKmer(&antiWord);

#ifdef _OPENMP
	thread = omp_get_thread_num();
#endif

	sprintf(lineBuffer, "ROADMAP %li\n", (long)seqID);
	appendLine(lineBuffer, thread);

	// Neglect any string shorter than WORDLENGTH :
	if (getLength(tString) < table->WORDLENGTH) {
#ifdef _OPENMP
		pushBuffer(thread);
#else
		velvetFprintf(file, "%s", annotationBuffer->str);
		resetStringBuffer(annotationBuffer);
#endif
		return;
	}

	// Fill in the initial word : 
	for (readNucleotideIndex = 0;
	     readNucleotideIndex < table->WORDLENGTH - 1;
	     readNucleotideIndex++) { 
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);
#ifdef COLOR
		reversePushNucleotide(&antiWord, nucleotide);
#else
		reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
	}

	while (readNucleotideIndex < getLength(tString)) {
		// Shift word:
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);
#ifdef COLOR
		reversePushNucleotide(&antiWord, nucleotide);
#else
		reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif

		sequenceID = seqID;
		coord = writeNucleotideIndex;

		if (table->double_strand) {
			if (compareKmers(&word, &antiWord) <= 0) {
				found =
				    findOrInsertOccurenceInSplayTable(&word,
								      &sequenceID,
								      &coord,
								      table,
								      sequenceIDs, 
								      coords,
								      readNucleotideIndex,
								      getLength(tString),
								      true);
			} else {
				sequenceID = -sequenceID;
				found =
				    findOrInsertOccurenceInSplayTable(&antiWord,
								      &sequenceID,
								      &coord,
								      table, 
								      sequenceIDs, 
								      coords,
								      readNucleotideIndex,
								      getLength(tString),
								      false);
				sequenceID = -sequenceID;
				if( sequenceID!=seqID && sequenceID<0 ){
				  printf("*** sequenceID = %d\n", sequenceID);
				  exit(-1);
				}
			}
		} else {
			if (!second_in_pair) {
				found =
				    findOrInsertOccurenceInSplayTable(&word,
								      &sequenceID,
								      &coord,
								      table,
								      sequenceIDs, 
								      coords,
								      readNucleotideIndex,
								      getLength(tString),
								      true);
			} else {
				sequenceID = -sequenceID;
				found =
				    findOrInsertOccurenceInSplayTable(&antiWord,
								      &sequenceID,
								      &coord,
								      table, 
								      sequenceIDs, 
								      coords,
								      readNucleotideIndex,
								      getLength(tString),
								      false);
				sequenceID = -sequenceID;
			}
		}

		if (!found) {
			writeNucleotideIndex++;
			if (!annotationClosed) {
				sprintf(lineBuffer, "%ld\t%lld\t%lld\t%lld\n",
					(long) referenceSequenceID, (long long) position,
					(long long) start, (long long) finish);
				appendLine(lineBuffer, thread);
			}
			annotationClosed = true;
		}
		// Otherwise create/complete annotation:
		else {
			// Forbidden k-mer
			if (sequenceID == 0) {
				break;
			}
			// Closed/inexistant annotation
			else if (annotationClosed) {
				referenceSequenceID = sequenceID;
				position = writeNucleotideIndex;
				start = finish = coord;

				if (referenceSequenceID > 0)
					finish++;
				else
					finish--;

				annotationClosed = false;
			}
			// Open annotation
			else if (sequenceID == referenceSequenceID
				 && coord == finish) {
				if (referenceSequenceID > 0)
					finish++;
				else
					finish--;
			}
			// Previous non corresponding annotation
			else {
				sprintf(lineBuffer, "%ld\t%lld\t%lld\t%lld\n",
					(long) referenceSequenceID, (long long) position,
					(long long) start, (long long) finish);
				appendLine(lineBuffer, thread);

				referenceSequenceID = sequenceID;
				position = writeNucleotideIndex;
				start = finish = coord;

				if (referenceSequenceID > 0)
					finish++;
				else
					finish--;
			}
		}

		readNucleotideIndex++;
	}

	if (!annotationClosed) {
		sprintf(lineBuffer, "%ld\t%lld\t%lld\t%lld\n",
			(long) referenceSequenceID, (long long) position,
			(long long) start, (long long) finish);
		appendLine(lineBuffer, thread);
	}
#ifdef _OPENMP
	pushBuffer(thread);
#else
	velvetFprintf(file, "%s", annotationBuffer->str);
	resetStringBuffer(annotationBuffer);
#endif
	return;
}

static void computeClearHSPs(TightString * array, FILE * seqFile, boolean second_in_pair, SplayTable * table, IDnum ** sequenceIDs, Coordinate ** coords, IDnum seqID) {
	Coordinate readNucleotideIndex = 0;
	Kmer word;
	Kmer antiWord;
	Kmer polyA;
	Nucleotide nucleotide;
	KmerOccurence * hit;
	char line[MAXLINE];
        char* start;
	char c;
	
	Coordinate mapCount = 0;
	Coordinate maxCount = 10;	
	IDnum * mapReferenceIDs = callocOrExit(maxCount, IDnum);
	Coordinate * mapCoords = callocOrExit(maxCount, Coordinate);
	long long_var;
	long long longlong_var;
	int penalty;
	TightString * tString;
	Coordinate length;

	clearKmer(&polyA);

#ifdef _OPENMP
	#pragma omp critical
	{
#endif
		// Get read ID:
		if (!fgets(line, MAXLINE, seqFile)) {
			puts("Incomplete Sequences file (computeHSPScores)");
#ifdef DEBUG
			abort();
#endif 
			exit(1);
		}
                start = strchr(line, '\t');
		tString = getTightStringInArray(array, seqID - 1);
		length = getLength(tString);
		*sequenceIDs = callocOrExit(length, IDnum);
		*coords = callocOrExit(length, Coordinate);

		// Parse file for mapping info
		while (seqFile && (c = getc(seqFile)) != EOF) {
			if (c == '>')
				break;

			fgets(line, MAXLINE, seqFile);

			if (c == 'M') {
				sscanf(line,"\t%li\t%lli\n", &long_var, &longlong_var);
				mapReferenceIDs[mapCount] = (IDnum) long_var;
				mapCoords[mapCount] = (Coordinate) longlong_var;

				if (++mapCount == maxCount) {
					maxCount *= 2;
					mapReferenceIDs = reallocOrExit(mapReferenceIDs, maxCount, IDnum);
					mapCoords = reallocOrExit(mapCoords, maxCount, Coordinate);
				}
			}
		}
#ifdef _OPENMP
	}
#endif

	// First pass for unambiguous hits
	// Fill in the initial word : 
	clearKmer(&word);
	clearKmer(&antiWord);
	for (readNucleotideIndex = 0;
	     readNucleotideIndex < table->WORDLENGTH - 1;
	     readNucleotideIndex++) { 
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);
#ifdef COLOR
		reversePushNucleotide(&antiWord, nucleotide);
#else
		reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
	}

	// Kill silly poly-T beginnings
	while (readNucleotideIndex < getLength(tString) && (compareKmers(&antiWord, &polyA) == 0 || compareKmers(&word, &polyA) == 0)) {
		nucleotide = getNucleotide(readNucleotideIndex++, tString);
		pushNucleotide(&word, nucleotide);
#ifdef COLOR
		reversePushNucleotide(&antiWord, nucleotide);
#else
		reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
	}

	while (readNucleotideIndex < getLength(tString)) {
		// Shift word:
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);

#ifdef COLOR
		reversePushNucleotide(&antiWord, nucleotide);
#else
		reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif

		if (table->double_strand) {
			if (compareKmers(&word, &antiWord) <= 0) {
				hit = findKmerInKmerOccurenceTable(&word, table->kmerOccurenceTable); 

				if (hit && (hit = getMostAppropriateHit(readNucleotideIndex, getLength(tString), true, hit, mapCount, mapReferenceIDs, mapCoords, table->WORDLENGTH))) 
					(*sequenceIDs)[readNucleotideIndex] = getKmerOccurenceNodeID(hit);
			} else {
				hit = findKmerInKmerOccurenceTable(&antiWord, table->kmerOccurenceTable); 

				if (hit && (hit = getMostAppropriateHit(readNucleotideIndex, getLength(tString), false, hit, mapCount, mapReferenceIDs, mapCoords, table->WORDLENGTH))) 
					(*sequenceIDs)[readNucleotideIndex] = -getKmerOccurenceNodeID(hit);
			}
		} else {
			if (!second_in_pair) {
				hit = findKmerInKmerOccurenceTable(&word, table->kmerOccurenceTable); 

				if (hit && (hit = getMostAppropriateHit(readNucleotideIndex, getLength(tString), true, hit, mapCount, mapReferenceIDs, mapCoords, table->WORDLENGTH))) 
					(*sequenceIDs)[readNucleotideIndex] = getKmerOccurenceNodeID(hit);
			} else {
				hit = findKmerInKmerOccurenceTable(&antiWord, table->kmerOccurenceTable); 

				if (hit && (hit = getMostAppropriateHit(readNucleotideIndex, getLength(tString), false, hit, mapCount, mapReferenceIDs, mapCoords, table->WORDLENGTH))) 
					(*sequenceIDs)[readNucleotideIndex] = -getKmerOccurenceNodeID(hit);
			}
		}

		if ((*sequenceIDs)[readNucleotideIndex]) {
			if ((*sequenceIDs)[readNucleotideIndex] > 0) 
				(*coords)[readNucleotideIndex] = getKmerOccurencePosition(hit) - readNucleotideIndex;
			else
				(*coords)[readNucleotideIndex] = getKmerOccurencePosition(hit) + readNucleotideIndex - getLength(tString) + 1;
		}
	
		// Barrier to flip-flopping
		if ((*sequenceIDs)[readNucleotideIndex - 1] != 0
		    && ((*sequenceIDs)[readNucleotideIndex] != (*sequenceIDs)[readNucleotideIndex - 1]
			|| (*coords)[readNucleotideIndex] != (*coords)[readNucleotideIndex - 1])) {
			// Break in continuity... skip k positions 
			(*sequenceIDs)[readNucleotideIndex] = 0;
			(*coords)[readNucleotideIndex] = -1;
			readNucleotideIndex++;

			for (penalty = 0; penalty < table->WORDLENGTH  - 1 && readNucleotideIndex < getLength(tString); penalty++) {
				nucleotide = getNucleotide(readNucleotideIndex, tString);
				pushNucleotide(&word, nucleotide);

#ifdef COLOR
				reversePushNucleotide(&antiWord, nucleotide);
#else
				reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
				(*sequenceIDs)[readNucleotideIndex] = 0;
				(*coords)[readNucleotideIndex] = -1;
				readNucleotideIndex++;
			}
		} else
			readNucleotideIndex++;
		
	}

	free(mapReferenceIDs);
	free(mapCoords);
}

void inputSequenceIntoSplayTable(TightString * array,
				 SplayTable * table,
				 FILE * file, FILE * seqFile,
				 boolean second_in_pair,
				 IDnum seqID)
{
	IDnum * sequenceIDs = NULL;
	Coordinate * coords = NULL;

	// If appropriate, get the HSPs on reference sequences
	if (table->kmerOccurenceTable) 
		computeClearHSPs(array, seqFile, second_in_pair, table, &sequenceIDs, &coords, seqID);

	// Go through read, eventually with annotations
	printAnnotations(sequenceIDs, coords, array, table, file, second_in_pair, seqID);

	// Clean up
	if (sequenceIDs) {
		free(sequenceIDs);
		free(coords);
	}
}

void inputReferenceIntoSplayTable(TightString * tString,
				 SplayTable * table, FILE * file, IDnum seqID, Mask * mask)
{
	IDnum currentIndex;
	Coordinate readNucleotideIndex = 0;
	Coordinate kmerIndex = 0;
	Kmer word;
	Kmer antiWord;
	Nucleotide nucleotide;
	Mask * currentMask = mask;
#ifdef _OPENMP
	char lineBuffer[MAXLINE];
#endif

	clearKmer(&word);
	clearKmer(&antiWord);

	currentIndex = seqID;
#ifdef _OPENMP
	sprintf(lineBuffer, "ROADMAP %li\n", (long)currentIndex);
	appendLine(lineBuffer, omp_get_thread_num());
#else
	velvetFprintf(file, "ROADMAP %li\n", (long)currentIndex);
#endif

	// Neglect any string shorter than WORDLENGTH :
	if (getLength(tString) < table->WORDLENGTH) {
		return;
	}

	// Fill in the initial word : 
	for (readNucleotideIndex = 0;
	     readNucleotideIndex < table->WORDLENGTH - 1;
	     readNucleotideIndex++) { 
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);
		if (table->double_strand) {
#ifdef COLOR
			reversePushNucleotide(&antiWord, nucleotide);
#else
			reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
		}
	}

	while (readNucleotideIndex < getLength(tString)) {
		// Shift word:
		nucleotide = getNucleotide(readNucleotideIndex, tString);
		pushNucleotide(&word, nucleotide);

		if (table->double_strand) {
#ifdef COLOR
			reversePushNucleotide(&antiWord, nucleotide);
#else
			reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif
		}

		// Check for gap masks:
		if (currentMask && currentMask->start - table->WORDLENGTH + 1 <= readNucleotideIndex) {
			while(currentMask && currentMask->finish + table->WORDLENGTH - 1 < readNucleotideIndex)
				currentMask = currentMask->next;

			if (currentMask && currentMask->finish + table->WORDLENGTH - 1 >= readNucleotideIndex) {
				readNucleotideIndex++;
				kmerIndex++;
				continue;
			}
		}

		// Record k-mer
		if (table->double_strand) {
			if (compareKmers(&word, &antiWord) <= 0)
				recordKmerOccurence(&word, currentIndex, 
						    kmerIndex,
						    table->kmerOccurenceTable);
			else
				recordKmerOccurence(&antiWord, -currentIndex, 
							kmerIndex,
							table->kmerOccurenceTable);
		} else {
			recordKmerOccurence(&word, currentIndex, 
					    kmerIndex,
					    table->kmerOccurenceTable);
		}
		readNucleotideIndex++;
		kmerIndex++;
	}

	return;
}

static Coordinate countReferenceKmers(ReadSet * reads, int wordLength) {
	IDnum readIndex;
	Coordinate length = 0;
	

	for (readIndex = 0; readIndex < reads->readCount && reads->categories[readIndex] == REFERENCE; readIndex++)	
	{
		Coordinate tmpLength = getLength(getTightStringInArray(reads->tSequences, readIndex));
		if (tmpLength >= wordLength)
			length += tmpLength - wordLength + 1;
	}

	return length;
}

Mask ** scanReferenceSequences(FILE * file, IDnum referenceSequenceCount) {
	Mask ** referenceMasks = callocOrExit(referenceSequenceCount, Mask*);
	IDnum index;
	char line[MAXLINE];
	char c;

	// Search sequences for masks
	for (index = 0; index < referenceSequenceCount; index++) {
		Mask * current = NULL;
		Coordinate position = 0;
		boolean openMask = false;

		// Read through header
		fgets(line, MAXLINE, file);

		// Read through sequence
		while ((c = getc(file))) {
			if (c == EOF || c == '>')
				break;
			else if (c == '\r' || c == '\n')
				continue;
			else if (c == 'n' || c == 'N') {
				if (openMask)
					current->finish++;
				else if (referenceMasks[index] == NULL) {
					referenceMasks[index] = newMask(position);
					current = referenceMasks[index];
				} else {
					current->next = newMask(position);
					current = current->next;		
				}
				openMask = true;
				position++;
			} else {
				openMask = false;
				position++;
			}
		}
	}

	return referenceMasks;
}

void inputSequenceArrayIntoSplayTableAndArchive(ReadSet * reads,
						SplayTable * table,
						char *filename, char* seqFilename)
{
	IDnum index;
	IDnum sequenceCount = reads->readCount;
	TightString *array;
	FILE *outfile = fopen(filename, "w");
	FILE *seqFile = NULL;
	IDnum kmerCount;
	IDnum referenceSequenceCount = 0;
	struct timeval start, end, diff;

	// DEBUG 
	Mask ** referenceMasks;

	if (outfile == NULL)
		exitErrorf(EXIT_FAILURE, true, "Couldn't write to file %s", filename);
	else
		velvetLog("Writing into roadmap file %s...\n", filename);

	// Count reference sequences
	for (index = 0; index < reads->readCount && reads->categories[index] == REFERENCE; index++)
		referenceSequenceCount++;

	velvetFprintf(outfile, "%ld\t%ld\t%i\t%hi\n", (long) sequenceCount, (long) referenceSequenceCount, table->WORDLENGTH, (short) table->double_strand);

	if (reads->tSequences == NULL)
		convertSequences(reads);

	gettimeofday(&start, NULL);
	array = reads->tSequences;

#ifdef _OPENMP
	if (omp_get_max_threads() == 1)
	{
		omp_set_num_threads(2);
		omp_set_nested(0);
	}
	else
		omp_set_nested(1);
	initAnnotationBuffers();
#else
	annotationBuffer = newStringBuffer(BUFFER_SIZE);
#endif

	if (referenceSequenceCount && (kmerCount = countReferenceKmers(reads, table->WORDLENGTH)) > 0) {
		table->kmerOccurenceTable = newKmerOccurenceTable(24 , table->WORDLENGTH);
		allocateKmerOccurences(kmerCount, table->kmerOccurenceTable);
		seqFile = fopen(seqFilename, "r");
		
		if (seqFile == NULL)
			exitErrorf(EXIT_FAILURE, true, "Couldn't write to file %s", seqFilename);
		else
			velvetLog("Reading mapping info from file %s\n", seqFilename);

		// Skip through reference headers quickly
		referenceMasks = scanReferenceSequences(seqFile, referenceSequenceCount);

#ifdef _OPENMP
		producing = 1;
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				bufferWritter(outfile);
			}
			#pragma omp section
			{
				#pragma omp parallel for
#endif
				for (index = 0; index < referenceSequenceCount; index++)
					inputReferenceIntoSplayTable(getTightStringInArray(array, index),
								     table, outfile, index + 1, referenceMasks[index]);

#ifdef _OPENMP
				for (index = omp_get_max_threads() - 1; index >= 0; index--)
					pushBufferCommit(index);
				producing = 0;
				#pragma omp flush(producing)
			}
		}
#endif

		if (maskMemory)
			destroyRecycleBin(maskMemory);
		sortKmerOccurenceTable(table->kmerOccurenceTable);
	}

	velvetLog("Inputting sequences...\n");

#ifdef _OPENMP
	producing = 1;
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			bufferWritter(outfile);
		}
		#pragma omp section
		{
			#pragma omp parallel for
#endif
			for (index = referenceSequenceCount; index < sequenceCount; index++)
			{
				boolean second_in_pair;

				// Progress report on screen
				if (index % 1000000 == 0) {
					velvetLog("Inputting sequence %li / %li\n",
						  (long)index, (long)sequenceCount);
					fflush(stdout);
				}

				// Test to make sure that all the reference reads are before all the other reads
				if (reads->categories[index] == REFERENCE) {
					velvetLog("Reference sequence placed after a non-reference read!\n");
					velvetLog(">> Please re-order the filenames in your command line so as "
						  "to have the reference sequence files before all the others\n");
#ifdef DEBUG 
					abort();
#endif 
					exit(0);
				}
				second_in_pair = reads->categories[index] % 2 && isSecondInPair(reads, index);

				// Hashing the reads
				inputSequenceIntoSplayTable(array, table,
							    outfile, seqFile,
							    second_in_pair, index + 1);
			}
#ifdef _OPENMP
			for (index = omp_get_max_threads() - 1; index >= 0; index--)
				pushBufferCommit(index);
			producing = 0;
			#pragma omp flush(producing)
		}
	}
	destroyAnnotationBuffers();
#else
	destroyStringBuffer(annotationBuffer, 1);
#endif

	gettimeofday(&end, NULL);
	timersub(&end, &start, &diff);
	velvetLog(" === Sequences loaded in %ld.%06ld s\n", diff.tv_sec, diff.tv_usec);

	fclose(outfile);
	if (seqFile)
		fclose(seqFile);

	//free(reads->tSequences);
	//reads->tSequences = NULL;
	//destroyReadSet(reads);
	velvetLog("Done inputting sequences\n");
}
