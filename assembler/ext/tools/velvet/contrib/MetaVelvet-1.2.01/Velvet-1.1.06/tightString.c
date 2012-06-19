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
#include "tightString.h"
#include "utility.h"

typedef unsigned char Codon;

struct tString_st {
	Codon *sequence;
	IDnum length;
}  ATTRIBUTE_PACKED;

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

TightString *getTightStringInArray(TightString * tString,
				   IDnum	 position)
{
	return tString + position;
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
	case 'N':
		nucleotideNum = Adenine;
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
	case 'n':
		nucleotideNum = Adenine;
	default:
		nucleotideNum = Adenine;
	}

	writeNucleotideNumber(nucleotideNum, codonPtr, position);
}

static void fillTightStringWithString(TightString * tString,
				      char *sequence,
				      Codon *newSequence)
{
	int index;

	tString->sequence = newSequence;
	for (index = 0; index < tString->length; index++)
		writeNucleotide(sequence[index],
				&(newSequence[index / 4]),
				index % 4);
	free(sequence);
}

//
// Creates a tightString from an array of normal strings
//
TightString *newTightStringArrayFromStringArray(char **sequences,
						IDnum sequenceCount,
						char **tSeqMem)
{
	IDnum sequenceIndex;
	Codon *tmp;
	TightString *tStringArray = mallocOrExit(sequenceCount, TightString);
	Coordinate totalLength = 0;
	int arrayLength;

	for (sequenceIndex = 0; sequenceIndex < sequenceCount; sequenceIndex++)
	{
		tStringArray[sequenceIndex].length = strlen (sequences[sequenceIndex]);
		arrayLength = (tStringArray[sequenceIndex].length + 3) / 4;
		totalLength += arrayLength;
	}
	*tSeqMem = callocOrExit (totalLength, char);
	tmp = (Codon*)*tSeqMem;
	for (sequenceIndex = 0; sequenceIndex < sequenceCount; sequenceIndex++)
	{
		fillTightStringWithString (&tStringArray[sequenceIndex],
					   sequences[sequenceIndex],
					   tmp);
		arrayLength = (tStringArray[sequenceIndex].length + 3) / 4;
		tmp += arrayLength;
	}

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

TightString *newTightString(Coordinate length)
{
	Coordinate arrayLength = length / 4;
	Coordinate index;
	TightString *newTString = mallocOrExit(1, TightString);
	if (length % 4 > 0)
		arrayLength++;

	newTString->length = length;
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

IDnum getLength(TightString * tString)
{
	return tString->length;
}

void destroyTightString(TightString * tString)
{
	free(tString->sequence);
	free(tString);
}

void setTightStringLength(TightString * tString, Coordinate length)
{
	Coordinate arrayLength = tString->length / 4;
	Coordinate newArrayLength = length / 4;
	if (tString->length % 4 > 0)
		arrayLength++;
	if (length % 4 > 0)
		newArrayLength++;

	if (newArrayLength > arrayLength) {
		tString->sequence =
		    reallocOrExit(tString->sequence,
			    newArrayLength, Codon);
		arrayLength = newArrayLength;
	}

	tString->length = length;
}

void exportTightString(FILE * outfile, TightString * sequence, IDnum index)
{
	Coordinate start, finish;
	char str[100];

	if (sequence == NULL)
		return;

	velvetFprintf(outfile, ">SEQUENCE_%ld_length_%lld\n", (long) index,
		(long long) getLength(sequence));

	start = 0;
	while (start <= getLength(sequence)) {
		finish = start + 60;
		readTightStringFragment(sequence, start, finish, str);
		velvetFprintf(outfile, "%s\n", str);
		start = finish;
	}

	fflush(outfile);
}
