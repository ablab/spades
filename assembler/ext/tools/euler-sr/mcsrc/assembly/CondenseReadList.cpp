/***************************************************************************
 * Title:          CondenseReadList.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "SimpleSequence.h"
#include "SeqReader.h"
#include "DNASequence.h"
#include "NumericHashedSpectrum.h"
#include "utils.h"
#include "IntegralTupleStatic.h"

using namespace std;

class IndexedRead : public IntegralTuple {
public:
	vector<string> readList;
	vector<string> titleList;
	vector<ssize_t> readListCount;
	static string read;
	static string title;
	ssize_t IncrementMult() {
		ssize_t i;
		for (i = 0; i < readList.size(); i++) {
			if (readList[i] == read) {
				//				cout << read << " " << readListCount[i] << endl;
				return ++readListCount[i];
			}			
		}
		readList.push_back(read);
		titleList.push_back(title);
		readListCount.push_back(1);
		return 1;
	}		
};

string IndexedRead::read;
string IndexedRead::title;

typedef	TupleHash<IndexedRead,11,15> IndexedReadHash;

void PrintUsage() {
	std::cout << "usage: integralCountSpectrum reads tupleSize spectrum [-printCount] [-binary]" << std::endl;
	std::cout << "     -skipGapped   Do not count sequences with 'N' or '.' in them." << std::endl;
}

int main(int argc, char* argv[]) {
	
	std::string readsFileName, readsOutFileName;
	int tupleSize;
	if (argc < 3) {
		PrintUsage();
		exit(1);
	}
	int argi = 1;
	readsFileName = argv[argi++];
	tupleSize     = atoi(argv[argi++]);
	readsOutFileName = argv[argi++];

	DNASequence read;
	DNASequence simpleRead;
	DNASequence readRC;
	ssize_t r = 0; 

	IndexedReadHash spectrum;

	std::ifstream readsIn;
	std::ofstream readsOut;

	openck(readsFileName, readsIn, std::ios::in);
	openck(readsOutFileName, readsOut, std::ios::out);
	//	IndexedRead::tupleSize = tupleSize;
	IndexedRead::SetTupleSize(tupleSize);
	while(SeqReader::GetSeq(readsIn, read, SeqReader::noConvert)) {
		simpleRead.seq = read.seq;
		simpleRead.length = read.length;

		ssize_t seqIsGapped = 0;
		ssize_t curPos;
		for (curPos = 0; curPos < read.length; curPos++ ){ 
			if (unmasked_nuc_index[read.seq[curPos]] >= 4) {
				seqIsGapped = 1;
				break;
			}
		}
		if (seqIsGapped) {
			continue;
		}

		
		ssize_t p;
		ssize_t newTupleStored = 0;
		ssize_t nextInvalid = -1;

		IndexedRead indexedRead;
		indexedRead.read.assign((const char*) simpleRead.seq, simpleRead.length);
		indexedRead.title = read.namestr;
		
		indexedRead.StringToTuple(&simpleRead.seq[0]);
		if (spectrum.hashTable.Store(indexedRead, spectrum.hashFunction, spectrum.updateFunction)) {
			// The first created the entry, this actually stores the sequence.
			spectrum.hashTable.Store(indexedRead, spectrum.hashFunction, spectrum.updateFunction);
		}
	}
	IndexedReadHash::HashTable::iterator end, it;

	ssize_t i;
	ssize_t chainPos = 0;
	SpectrumHash::HashTable::Page *page;
	SpectrumHash::HashTable::Data *data, *pageEnd;
	for (i = 0; i < spectrum.hashTable.size; i++ ) {
		end = spectrum.hashTable.table[i].End();
		it  = spectrum.hashTable.table[i].Begin();
		
		page = spectrum->hashTable.table[i].head;
	
		while (page != NULL) {
			pageEnd = &(*page).page[page->size];
			ssize_t r;
			for (data = &(*page).page[0]; data != pageEnd; ++data) {
				titlestrm << readIt.titleList[r] << " mult=" << readIt.readListCount[r];
				readsOut << ">" << titlestrm.str() << endl;
				readsOut << readIt.readList[r] << endl;

			}
			page = page->next;
		}
		while (it != end) {
			IndexedRead readIt = *it;
			ssize_t r;
			cout  << i << " " << chainPos << " " << readIt.readList.size() << endl;
			++chainPos;
			for (r = 0; r < readIt.readList.size(); r++ ) {
				std::stringstream titlestrm;
			}
			++it;
		}
	}	
}
		
