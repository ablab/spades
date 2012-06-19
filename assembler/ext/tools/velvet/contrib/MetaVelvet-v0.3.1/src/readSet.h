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
	TightString **tSequences;
	char **labels;
	Quality **confidenceScores;
	Probability **kmerProbabilities;
	IDnum *mateReads;
	Category *categories;
	IDnum readCount;
};

ReadSet *newReadSet();
ReadSet *newReadSetAroundTightStringArray(TightString ** array,
					  IDnum length);

Coordinate *getSequenceLengths(ReadSet * reads, int wordLength);
Coordinate *getSequenceLengthsFromFile(char *filename, int wordLength);

void concatenateReadSets(ReadSet * A, ReadSet * B);

void convertSequences(ReadSet * rs);

void convertConfidenceScores(ReadSet * rs, int WORDLENGTH);

void categorizeReads(ReadSet * reads, Category category);
void simplifyReads(ReadSet * reads);

// Exports a .sed script allowing to transform internal IDs to the original ones
void exportIDMapping(char *filename, ReadSet * reads);

ReadSet *importReadSet(char *filename);
void exportReadSet(char *filename, ReadSet * reads);

// The overall argument parser and file reader for the hash function
void parseDataAndReadFiles(char * filename, int argc, char **argv, boolean * double_strand);

void logInstructions(int argc, char **argv, char *directory);

// Read pairing info
void createReadPairingArray(ReadSet* reads);
boolean pairUpReads(ReadSet * reads, Category cat);
void detachDubiousReads(ReadSet * reads, boolean * dubiousReads);

void destroyReadSet(ReadSet * reads);
#endif
