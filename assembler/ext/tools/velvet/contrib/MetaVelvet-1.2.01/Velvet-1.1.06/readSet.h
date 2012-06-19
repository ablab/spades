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
#ifndef _READSET_H_
#define _READSET_H_

struct readSet_st {
	char **sequences;
	TightString *tSequences;
	char **labels;
	char *tSeqMem;
	Quality **confidenceScores;
	Probability **kmerProbabilities;
	IDnum *mateReads;
	Category *categories;
	unsigned char *secondInPair;
	IDnum readCount;
};

ReadSet *newReadSet();

ShortLength *getSequenceLengths(ReadSet * reads, int wordLength);

void convertSequences(ReadSet * rs);

ReadSet *importReadSet(char *filename);

// The overall argument parser and file reader for the hash function
void parseDataAndReadFiles(char * filename, int argc, char **argv, boolean * double_strand, boolean * noHash);

void logInstructions(int argc, char **argv, char *directory);

// Read pairing info
void createReadPairingArray(ReadSet* reads);
int pairedCategories(ReadSet * reads);
boolean isSecondInPair(ReadSet * reads, IDnum index);
void detachDubiousReads(ReadSet * reads, boolean * dubiousReads);

void destroyReadSet(ReadSet * reads);
#endif
