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
#include <time.h>

#include "globals.h"
#include "readSet.h"
#include "splay.h"
#include "tightString.h"
#include "crc.h"
#include "utility.h"
#include "kmer.h"

struct splayTable_st {
	SplayTree **table;
	IDnum lastIndex;
	int WORDLENGTH;
};

SplayTable *newSplayTable(int WORDLENGTH)
{
	SplayTable *splayTable = mallocOrExit(1, SplayTable);
	splayTable->WORDLENGTH = WORDLENGTH;
	splayTable->table = callocOrExit(CRC_HASH_BUCKETS, SplayTree *);
	splayTable->lastIndex = 0;
	return splayTable;
}

void destroySplayTable(SplayTable * splayTable)
{
	puts("Destroying splay table");

	destroyAllSplayTrees();
	free(splayTable->table);
	free(splayTable);

	puts("Splay table destroyed");
}

static int hash_kmer(Kmer * kmer)
{
	return crc32_v((char *) kmer, KMER_BYTE_SIZE);
}

static boolean findOrInsertOccurenceInSplayTable(Kmer * kmer, IDnum * seqID,
						 Coordinate * position,
						 SplayTable * table)
{
	if (table == NULL) {
		puts("NULL table!");
		exit(1);
	}

	return findOrInsertOccurenceInSplayTree(kmer, seqID, position,
						&table->
						table[hash_kmer(kmer)]);
}

void inputSequenceIntoSplayTable(TightString * tString,
				 SplayTable * table, FILE * file, boolean double_strand, boolean second_in_pair)
{
	IDnum currentIndex;
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

	clearKmer(&word);
	clearKmer(&antiWord);

	table->lastIndex++;

	currentIndex = table->lastIndex;
	fprintf(file, "ROADMAP %d\n", currentIndex);

	// Neglect any string shorter than WORDLENGTH :
	if (getLength(tString) < table->WORDLENGTH) {
		destroyTightString(tString);
		return;
	}
	// Fill in the initial word : 
	for (readNucleotideIndex = 0;
	     readNucleotideIndex < table->WORDLENGTH - 1;
	     readNucleotideIndex++) { 
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

		sequenceID = currentIndex;
		coord = writeNucleotideIndex;

		if (double_strand) {
			if (compareKmers(&word, &antiWord) <= 0) {
				found =
				    findOrInsertOccurenceInSplayTable(&word,
								      &sequenceID,
								      &coord,
								      table);
			} else {
				sequenceID = -sequenceID;
				found =
				    findOrInsertOccurenceInSplayTable(&antiWord,
								      &sequenceID,
								      &coord,
								      table);
				sequenceID = -sequenceID;
			}
		} else {
			if (!second_in_pair) {
				found =
				    findOrInsertOccurenceInSplayTable(&word,
								      &sequenceID,
								      &coord,
								      table);
			} else {
				sequenceID = -sequenceID;
				found =
				    findOrInsertOccurenceInSplayTable(&antiWord,
								      &sequenceID,
								      &coord,
								      table);
				sequenceID = -sequenceID;
			}
		}

		if (!found) {
			writeNucleotideIndex++;
			if (!annotationClosed)
				fprintf(file, "%ld\t%lld\t%lld\t%lld\n",
					(long) referenceSequenceID, (long long) position,
					(long long) start, (long long) finish);
			annotationClosed = true;
		}
		// Other wise create/complete annotation:
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
				fprintf(file, "%ld\t%lld\t%lld\t%lld\n",
					(long) referenceSequenceID, (long long) position,
					(long long) start, (long long) finish);

				referenceSequenceID = sequenceID;
				position = writeNucleotideIndex;
				start = finish = coord;

				if (referenceSequenceID > 0)
					finish++;
				else
					finish--;
			}
		}
	}

	if (!annotationClosed)
		fprintf(file, "%ld\t%lld\t%lld\t%lld\n",
			(long) referenceSequenceID, (long long) position,
			(long long) start, (long long) finish);

	destroyTightString(tString);
	return;
}

void inputSequenceArrayIntoSplayTableAndArchive(ReadSet * reads,
						SplayTable * table,
						char *filename, 
						boolean double_strand)
{
	IDnum index;
	IDnum sequenceCount = reads->readCount;
	TightString **array;
	FILE *outfile = fopen(filename, "w");
	boolean second_in_pair = false;

	if (outfile == NULL)
		exitErrorf(EXIT_FAILURE, true, "Couldn't write to file %s", filename);
	else
		printf("Writing into roadmap file %s...\n", filename);

	fprintf(outfile, "%ld\t%i\t%hi\n", (long) sequenceCount, table->WORDLENGTH, (short) double_strand);

	if (reads->tSequences == NULL)
		convertSequences(reads);

	array = reads->tSequences;

	puts("Inputting sequences...");
	for (index = 0; index < sequenceCount; index++) {
		if (index % 100000 == 0) {
			printf("Inputting sequence %d / %d\n", index,
			       sequenceCount);
			fflush(stdout);
		}
		inputSequenceIntoSplayTable(array[index], table, outfile, double_strand, second_in_pair);

		if (reads->categories[index] % 2) 
			second_in_pair = (second_in_pair? false : true);
		else 
			second_in_pair = false;
	}

	fclose(outfile);

	free(reads->tSequences);
	reads->tSequences = NULL;
	destroyReadSet(reads);
	puts("Done inputting sequences");
}

void inputMaskIntoSplayTable(TightString * tString, SplayTable * table)
{
	Coordinate readNucleotideIndex = 0;
	Kmer word;
	Kmer antiWord;
	IDnum sequenceID = 0;
	Coordinate coord = 0;
	Nucleotide nucleotide;

	clearKmer(&word);
	clearKmer(&antiWord);

	// Neglect any string shorter than WORDLENGTH :
	if (getLength(tString) < table->WORDLENGTH) {
		destroyTightString(tString);
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
		nucleotide = getNucleotide(readNucleotideIndex++, tString);
		pushNucleotide(&word, nucleotide);
#ifdef COLOR
		reversePushNucleotide(&antiWord, nucleotide);
#else
		reversePushNucleotide(&antiWord, 3 - nucleotide);
#endif

		if (compareKmers(&word, &antiWord) <= 0)
			findOrInsertOccurenceInSplayTable(&word,
							  &sequenceID,
							  &coord, table);
		else
			findOrInsertOccurenceInSplayTable(&antiWord,
							  &sequenceID,
							  &coord, table);
	}

	destroyTightString(tString);
	return;
}

void inputMaskArrayIntoSplayTable(ReadSet * reads, SplayTable * table)
{
	IDnum index;
	IDnum sequenceCount = reads->readCount;
	TightString **array;

	if (reads->tSequences == NULL)
		convertSequences(reads);

	array = reads->tSequences;

	puts("Loading masks...");
	for (index = 0; index < sequenceCount; index++)
		inputMaskIntoSplayTable(array[index], table);

	free(reads->tSequences);
	reads->tSequences = NULL;
	destroyReadSet(reads);
	puts("Done loading masks");
}
