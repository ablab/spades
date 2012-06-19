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
#include <string.h>
#include <stdio.h>

#include "globals.h"
#include "utility.h"

typedef unsigned char Codon;

struct tString_st {
	Codon *sequence;
	Coordinate length;
	Coordinate arrayLength;
};

static const Nucleotide Adenine = 0;
static const Nucleotide Cytosine = 1;
static const Nucleotide Guanine = 2;
static const Nucleotide Thymine = 3;

// Binary 11111100
static const Codon FILTER0 = ~((Codon) 3);
// Binary 11110011
static const Codon FILTER1 = ~(((Codon) 3) << 2);
// Binary 11001111
static const Codon FILTER2 = ~(((Codon) 3) << 4);
// Binary 00111111
static const Codon FILTER3 = (Codon) ~ ((Codon) 3 << 6);

//
// Adds a number into the Codon pointed to by codonPtr at the desired
// position (0, 1, 2, or 3);
//
void
writeNucleotideNumber(Nucleotide nucleotide, Codon * codonPtr,
		      Coordinate position)
{
	if (position == 3) {
		*codonPtr &= FILTER3;
		*codonPtr += nucleotide << 6;
	} else if (position == 2) {
		*codonPtr &= FILTER2;
		*codonPtr += nucleotide << 4;
	} else if (position == 1) {
		*codonPtr &= FILTER1;
		*codonPtr += nucleotide << 2;
	} else if (position == 0) {
		*codonPtr &= FILTER0;
		*codonPtr += nucleotide;
	}
}

//
// Adds a nucleotide into the Codon pointed to by codonPtr at the desired
// position (0, 1, 2, or 3);
//
Nucleotide charToNucleotide(char c)
{
	switch (c) {
	case 'A':
		return Adenine;
	case 'C':
		return Cytosine;
	case 'G':
		return Guanine;
	case 'T':
		return Thymine;
	case 'a':
		return Adenine;
	case 'c':
		return Cytosine;
	case 'g':
		return Guanine;
	case 't':
		return Thymine;
	case '\n':
		return '\n';
	default:
		return Adenine;
	}

}

//
// Adds a nucleotide into the Codon pointed to by codonPtr at the desired
// position (0, 1, 2, or 3);
//
void writeNucleotide(Nucleotide nucleotide, Codon * codonPtr, int position)
{
	int nucleotideNum;

	switch (nucleotide) {
	case 'A':
		nucleotideNum = Adenine;
		break;
	case 'C':
		nucleotideNum = Cytosine;
		break;
	case 'G':
		nucleotideNum = Guanine;
		break;
	case 'T':
		nucleotideNum = Thymine;
		break;
	case 'a':
		nucleotideNum = Adenine;
		break;
	case 'c':
		nucleotideNum = Cytosine;
		break;
	case 'g':
		nucleotideNum = Guanine;
		break;
	case 't':
		nucleotideNum = Thymine;
		break;
	default:
		nucleotideNum = Adenine;
	}

	writeNucleotideNumber(nucleotideNum, codonPtr, position);
}

//
// Creates a tightString from a tradionnal string of A,T,G, and C of length size
//
TightString *newTightStringFromString(char *sequence)
{
	TightString *newTString = mallocOrExit(1, TightString);

	int size = (int) strlen(sequence);
	int arrayLength = size / 4;
	int index;

	if (size % 4 > 0)
		arrayLength++;

	newTString->length = size;
	newTString->arrayLength = arrayLength;
	newTString->sequence = callocOrExit(arrayLength, Codon);

	for (index = 0; index < arrayLength; index++)
		newTString->sequence[index] = 0;

	for (index = 0; index < size; index++)
		writeNucleotide(sequence[index],
				&(newTString->sequence[index / 4]),
				index % 4);

	free(sequence);
	return newTString;
}

//
// Creates a tightString from an array of normal strings
//
TightString **newTightStringArrayFromStringArray(char **sequences,
						 IDnum sequenceCount)
{
	IDnum sequenceIndex;
	TightString **tStringArray =
	    mallocOrExit(sequenceCount, TightString *);

	for (sequenceIndex = 0; sequenceIndex < sequenceCount;
	     sequenceIndex++)
		tStringArray[sequenceIndex] =
		    newTightStringFromString(sequences[sequenceIndex]);

	free(sequences);
	return tStringArray;
}

char readNucleotide(Nucleotide nucleotide)
{
	switch (nucleotide) {
	case 0:
		return 'A';
	case 1:
		return 'C';
	case 2:
		return 'G';
	case 3:
		return 'T';
	}

	return '?';
}

char *readTightString(TightString * tString)
{
	Coordinate index, index4;
	char *string;
	Codon codon;

	if (tString == NULL || tString->length == 0) {
		string = callocOrExit(5, char);
		strcpy(string, "VOID");
		return string;
	}

	string = callocOrExit(tString->length + 1, char);

	for (index = 0; index < tString->length / 4; index++) {
		index4 = index << 2;
		codon = tString->sequence[index];
		string[index4] = readNucleotide(codon & 3);
		string[index4 + 1] = readNucleotide((codon & 12) >> 2);
		string[index4 + 2] = readNucleotide((codon & 48) >> 4);
		string[index4 + 3] = readNucleotide((codon & 192) >> 6);
	}

	index4 = index << 2;
	codon = tString->sequence[index];

	switch (tString->length % 4) {
	case 3:
		string[index4 + 3] = readNucleotide((codon & 192) >> 6);
	case 2:
		string[index4 + 2] = readNucleotide((codon & 48) >> 4);
	case 1:
		string[index4 + 1] = readNucleotide((codon & 12) >> 2);
	case 0:
		string[index4] = readNucleotide(codon & 3);
	}

	string[tString->length] = '\0';

	return string;
}

Nucleotide getNucleotide(Coordinate nucleotideIndex, TightString * tString)
{
	Codon codon = tString->sequence[nucleotideIndex / 4];

	switch (nucleotideIndex % 4) {
	case 3:
		return (codon & 192) >> 6;
	case 2:
		return (codon & 48) >> 4;
	case 1:
		return (codon & 12) >> 2;
	case 0:
		return (codon & 3);
	}

	return '?';
}

void readTightStringFragment(TightString * tString, Coordinate start,
			     Coordinate finish, char *string)
{
	Coordinate index;
	Coordinate inFinish = finish;

	if (start >= tString->length) {
		string[0] = '\0';
		return;
	}

	if (inFinish > tString->length)
		inFinish = tString->length;

	for (index = start; index < inFinish; index++) {
		string[index - start] =
		    readNucleotide(getNucleotide(index, tString));
	}

	string[inFinish - start] = '\0';
}

char getNucleotideChar(Coordinate nucleotideIndex, TightString * tString)
{
	Codon codon;

	codon = tString->sequence[nucleotideIndex / 4];

	switch (nucleotideIndex % 4) {
	case 3:
		return readNucleotide((codon & 192) >> 6);
	case 2:
		return readNucleotide((codon & 48) >> 4);
	case 1:
		return readNucleotide((codon & 12) >> 2);
	case 0:
		return readNucleotide((codon & 3));
	}

	return '?';
}

char getInverseNucleotideChar(Coordinate nucleotideIndex,
			      TightString * tString)
{
	Codon codon = tString->sequence[nucleotideIndex / 4];

	switch (nucleotideIndex % 4) {
#ifndef COLOR
	case 3:
		return readNucleotide(3 - ((codon & 192) >> 6));
	case 2:
		return readNucleotide(3 - ((codon & 48) >> 4));
	case 1:
		return readNucleotide(3 - ((codon & 12) >> 2));
	case 0:
		return readNucleotide(3 - ((codon & 3)));
#else
	case 3:
		return readNucleotide(((codon & 192) >> 6));
	case 2:
		return readNucleotide(((codon & 48) >> 4));
	case 1:
		return readNucleotide(((codon & 12) >> 2));
	case 0:
		return readNucleotide(((codon & 3)));
#endif
	}

	return '?';
}

TightString *newTightString(Coordinate length)
{
	Coordinate arrayLength = length / 4;
	Coordinate index;
	TightString *newTString = mallocOrExit(1, TightString);
	if (length % 4 > 0)
		arrayLength++;

	newTString->length = length;
	newTString->arrayLength = arrayLength;
	newTString->sequence = callocOrExit(arrayLength, Codon);

	for (index = 0; index < arrayLength; index++)
		newTString->sequence[index] = 0;

	return newTString;
}

void
writeNucleotideAtPosition(Nucleotide nucleotide, Coordinate position,
			  TightString * tString)
{
	if (position >= tString->length)
		return;

	writeNucleotideNumber(nucleotide,
			      &tString->sequence[position / 4],
			      position % 4);
}

void trimTightString(TightString * tString, Coordinate length)
{
	Coordinate newArrayLength = length / 4;
	if (length % 4 == 0)
		newArrayLength++;

	tString->length = length;
	tString->arrayLength = newArrayLength;
	tString->sequence = 
	    reallocOrExit(tString->sequence, newArrayLength, Codon);
}

Coordinate getLength(TightString * tString)
{
	return tString->length;
}

TightString **concatenateTightStringArrays(TightString ** array1,
					   TightString ** array2,
					   IDnum size1, IDnum size2)
{
	TightString **unionArray;
	IDnum index;

	if (array1 == NULL)
		return array2;

	if (array2 == NULL)
		return array1;

	unionArray =
	    reallocOrExit(array1, size1 + size2, TightString *);

	for (index = 0; index < size2; index++)
		unionArray[size1 + index] = array2[index];

	free(array2);

	return unionArray;
}

void destroyTightString(TightString * tString)
{
	free(tString->sequence);
	free(tString);
}

void destroyTightStringArray(TightString ** array, IDnum sequenceCount)
{
	IDnum index;
	for (index = 0; index < sequenceCount; index++)
		destroyTightString(array[index]);
	free(array);
}

void setTightStringLength(TightString * tString, Coordinate length)
{
	Coordinate newArrayLength = length / 4;
	if (length % 4 > 0)
		newArrayLength++;

	if (newArrayLength > tString->arrayLength) {
		tString->sequence =
		    reallocOrExit(tString->sequence,
			    newArrayLength, Codon);
		tString->arrayLength = newArrayLength;
	}

	tString->length = length;
}

// Shortens reads to a fixed size (good for Solexa where errors are markedly towards the end)
void trimTightStringArray(TightString ** array, IDnum sequenceCount,
			  Coordinate length)
{
	IDnum index;

	for (index = 0; index < sequenceCount; index++)
		trimTightString(array[index], length);
}

void trimTightStringArraySanger(TightString ** array, IDnum sequenceCount,
				Coordinate min, Coordinate max)
{
	IDnum index;

	for (index = 0; index < sequenceCount; index++) {
		if (getLength(array[index]) > max)
			trimTightString(array[index], max);
		else if (getLength(array[index]) < min)
			trimTightString(array[index], 0);
	}
}

void clipTightString(TightString * tString, Coordinate start,
		     Coordinate finish)
{
	Coordinate position;
	Coordinate newLength = finish - start;

	for (position = 0; position < newLength; position++)
		writeNucleotideAtPosition(getNucleotide
					  (position + start, tString),
					  position, tString);

	trimTightString(tString, newLength);
}

Nucleotide getNucleotideFromString(Coordinate nucleotideIndex,
				   char *string)
{
	char letter = string[nucleotideIndex];

	switch (letter) {
	case 'A':
		return Adenine;
	case 'C':
		return Cytosine;
	case 'G':
		return Guanine;
	case 'T':
		return Thymine;
	default:
		return Adenine;
	}
}

void exportTightString(FILE * outfile, TightString * sequence, IDnum index)
{
	Coordinate start, finish;
	char str[100];

	if (sequence == NULL)
		return;

	fprintf(outfile, ">SEQUENCE_%ld_length_%lld\n", (long) index,
		(long long) getLength(sequence));

	start = 0;
	while (start <= getLength(sequence)) {
		finish = start + 60;
		readTightStringFragment(sequence, start, finish, str);
		fprintf(outfile, "%s\n", str);
		start = finish;
	}

	fflush(outfile);
}

void exportSequenceArray(char *filename, TightString ** array,
			 IDnum sequenceCount)
{
	IDnum index;
	FILE *outfile = fopen(filename, "w+");

	if (outfile == NULL) {
		puts("Couldn't open file, sorry");
		return;
	} else
		printf("Writing into file: %s\n", filename);

	for (index = 0; index < sequenceCount; index++) {
		exportTightString(outfile, array[index], index);
	}

	fclose(outfile);

	puts("Done");
}
