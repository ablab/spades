/***************************************************************************
 * Title:          ElandToFastq.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

using namespace std;

ssize_t AdvancePastTab(char *&pos, char*end) {
	char *start = pos;
	// pass the current tab.
	while (pos < end and
				 *pos != '\t')
		++pos;
	++pos;
	return start - pos;
}

ssize_t AdvanceWhitespace(char *&pos, char*end) {
	char *start = pos;
	while (pos < end and 
				 (*pos == ' ' or *pos == '\t'))
		++pos;
	return pos - start;
}
ssize_t AdvanceNonWhitespace(char *&pos, char*end) {
	char *start = pos;
	while (pos < end and
				 *pos != ' ' and 
				 *pos != '\t')
		++pos;
	return pos - start;
}

void PrintUsage() {
		cout << "usage: elandToFastq eland.txt reads.fastq  [-printNotPure]" << endl;
		cout << "   Translate an eland file to fastq.  " << endl;
		cout << "   -printNotPure  Print all reads, even those that do not" <<endl
				 << "                  pass the purity check." << endl;
}

int main(int argc, char* argv[]) {

	string elandInName, fastqOutName;

	ssize_t onlyPure = 1;
	if (argc < 3) {
		PrintUsage();
		exit(0);
	}

	elandInName  = argv[1];
	fastqOutName = argv[2];
	int argi = 3;
	while (argi < argc) {
		if (strcmp(argv[argi], "-printNotPure") == 0) {
			onlyPure = 0;
		}
		else {
			PrintUsage();
			cout << "bad option: " << argv[argi]<<endl;
		}
		++argi;
	}
	//UNUSED// ssize_t seqIndex = 0;
	
	ifstream elandIn;
	ofstream fastqOut;

	openck(elandInName, elandIn, std::ios::in);
	openck(fastqOutName, fastqOut, std::ios::out);

	char *title = (char *) NULL;
	ssize_t titleLength = 0;
	string line;
	char *seqPtr = (char *) NULL, *qualPtr = (char *) NULL;
	ssize_t seqLength = 0;
	string elandLine;

	/*

		A description of the eland output.
   1.  Lane
   2. Tile
   3. X Coordinate of cluster
   4. Y Coordinate of cluster
   5. Index string (Blank for a non-indexed run)
   6. Read number (1 or 2 for paired-read analysis, blank for a single-read analysis)
   7. Read
   8. Quality string--In symbolic ASCII format (ASCII character code = quality value + 64)
   9. Match chromosome--Name of chromosome match OR code indicating why no match resulted
  10. Match Contig--Gives the contig name if there is a match and the match chromosome is split into contigs (Blank if no match found)
  11. Match Position--Always with respect to forward strand, numbering starts at 1 (Blank if no match found)
  12. Match Strand--"F" for forward, "R" for reverse (Blank if no match found)
  13. Match Descriptor--Concise description of alignment (Blank if no match found)
          * A numeral denotes a run of matching bases
          * A letter denotes substitution of a nucleotide: For a 35 base read, "35" denotes an exact match and "32C2" denotes substitution of a "C" at the 33rd position 
  14. Single-Read Alignment Score--Alignment score of a single-read match, or for a paired read, alignment score of a read if it were treated as a single read. Blank if no match found; any scores less than 4 should be considered as aligned to a repeat
  15. Paired-Read Alignment Score--Alignment score of a paired read and its partner, taken as a pair. Blank if no match found; any scores less than 4 should be considered as aligned to a repeat
  16. Partner Chromosome--Name of the chromosome if the read is paired and its partner aligns to another chromosome (Blank for single-read analysis)
  17. Partner Contig--Not blank if read is paired and its partner aligns to another chromosome and that partner is split into contigs (Blank for single-read analysis)
  18. Partner Offset--If a partner of a paired read aligns to the same chromosome and contig, this number, added to the Match Position, gives the alignment position of the partner (Blank for single-read analysis)
  19. Partner Strand--To which strand did the partner of the paired read align? "F" for forward, "R" for reverse (Blank if no match found, blank for single-read analysis)
  20. Filtering--Did the read pass quality filtering? "Y" for yes, "N" for no 
	*/
	ssize_t numNotFiltered = 0;
	ssize_t totalEntries = 0;
	while(elandIn) {
		getline(elandIn, elandLine);
		if (elandLine.size() == 0)
			break;

		++totalEntries;

		// compute the title length;
		//UNUSED// ssize_t p;
		ssize_t elandLineLength = elandLine.size();
		char *elandBegin = (char*) elandLine.c_str();
		char *elandEnd = elandBegin + elandLineLength;
		char *elandPtr = elandBegin;
		//UNUSED// char *titlePtr = elandBegin;

		ssize_t i;
		for (i = 0; i < 7; i++) {
			AdvancePastTab(elandPtr, elandEnd); 
		}
		AdvanceNonWhitespace(elandPtr, elandEnd);
		ssize_t elandTitleLength = elandPtr - elandBegin;
		if (elandTitleLength + 2 > titleLength) {
			if (titleLength != 0)
				delete[] title;
			titleLength = elandTitleLength + 2;
			title = new char[titleLength];
		}
		memcpy(&title[1], elandBegin, elandTitleLength);
		i = 0;
		ssize_t tabNumber = 0;
		while (i < titleLength) {
			if (title[i] == '\t') {
				// handle read pairing printing.
				if (tabNumber == 6) {
					if (i < titleLength-1 and 
							(title[i+1] == '1' or 
							 title[i+1] == '2')) {
						title[i] = '/';
					}
				}
				else {
					title[i] = ':';
				}
				++tabNumber;
			}
			i++;
		}
		

		title[elandTitleLength+1] = 0;
		// Move forward to the beginning of the sequence
		AdvanceWhitespace(elandPtr, elandEnd);
		char *elandSeqBegin = elandPtr;
		AdvanceNonWhitespace(elandPtr, elandEnd);
		ssize_t elandSeqLength = elandPtr - elandSeqBegin;
		if (elandSeqLength + 1 > seqLength) {
			if (seqLength != 0) {
				delete[] seqPtr;
				delete[] qualPtr;
			}
			seqLength = elandSeqLength + 1;
			seqPtr  = new char[elandSeqLength + 1];
			qualPtr = new char[elandSeqLength + 1];
		}
		memcpy(seqPtr, elandSeqBegin, elandSeqLength);
		seqPtr[seqLength] = 0;

		AdvanceWhitespace(elandPtr, elandEnd);
		elandSeqBegin = elandPtr;
		AdvanceNonWhitespace(elandPtr, elandEnd);
		assert(elandPtr - elandSeqBegin == elandSeqLength);
		memcpy(qualPtr, elandSeqBegin, elandSeqLength);

		// Skip all the alignment columns, and get to the
		// purity filter.
		for (i = 0; i < 12; i++) 
			AdvancePastTab(elandPtr, elandEnd);
		assert(elandPtr < elandEnd);
		
		if (!onlyPure or *elandPtr == 'Y') {
			title[0] = '@';
			fastqOut << title << endl;
			fastqOut << seqPtr << endl;
			title[0] = '+';
			fastqOut << title << endl;
			fastqOut << qualPtr << endl;
			numNotFiltered++;
		}

	}
	if (titleLength > 0)
		delete[] title;
	if (seqLength > 0)
		delete[] seqPtr;

	cout << "Of " << totalEntries << " reads, " << numNotFiltered 
			 << " passed the purity filter. " << endl;
	return 0;
}
