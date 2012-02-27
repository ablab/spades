/***************************************************************************
 * Title:          CountIntegralTuples.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "SeqReader.h"
#include "DNASequence.h"
#include "SeqUtils.h"
#include "utils.h"
#include "IntegralTupleStatic.h"


ssize_t  ValidTuple(unsigned char* seq, ssize_t pos, ssize_t seqLength, ssize_t tupleLength);
void AdvanceToValid(unsigned char *seq, ssize_t curPos, ssize_t seqLength, ssize_t tupleLength,
										CountedIntegralTuple &tuple, ssize_t &nextPos);

void PrintUsage() {
		std::cout << "usage: countIntegralTuples tupleList tupleSize seqFile [-skipGapped]" << std::endl;
}
int main(int argc, char* argv[]) {

	std::string tupleListName, outTupleListName;
	std::string seqFileName;
	int tupleSize;
	if (argc < 4) {
		PrintUsage();
		exit(1);
	}
	tupleListName = argv[1];
	tupleSize     = atoi(argv[2]);
	seqFileName   = argv[3];
	int argi = 4;
	ssize_t skipGapped = 0;

	while (argi < argc) {
		if (strcmp(argv[argi], "-skipGapped") == 0) {
			skipGapped = 1;
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
		}
		argi++;
	}

	std::ifstream listIn, seqIn;

	openck(tupleListName, listIn, std::ios::in | std::ios::binary);

	// TODO: tupleSize was added to the file format late; can remove it from command-line
	ssize_t tupleSize_SSZT;
	listIn.read((char*) &tupleSize_SSZT, sizeof(tupleSize_SSZT));
	if (tupleSize != tupleSize_SSZT) {
		std::cerr << "Warning: command-line parameter tupleSize=" << tupleSize
							<< " disagrees with binary file tupleSize=" << tupleSize_SSZT << std::endl;
	}
	CountedIntegralTuple::SetTupleSize(tupleSize);

	ssize_t nTuples;
	listIn.read((char*) &nTuples, sizeof(ssize_t));
	std::cout << "reading: " << nTuples << " tuples." << std::endl;
	std::cout << "they will require: " << sizeof(CountedIntegralTuple) * nTuples << " bytes." << std::endl;
	_SZT_ i;
	CountedIntegralTupleList tupleList;
	tupleList.resize(nTuples);

	for (i = 0; i < nTuples; i++) {
		if (i % 100000 == 0) 
			std::cout << i << std::endl;
		listIn.read((char *) &(tupleList[i].tuple), sizeof(LongTuple));
		//		LongTuple tuple;
		//		listIn.read((char*) &tuple, sizeof(LongTuple));
		//		tupleList[i].tuple = tuple;
	}
	
	openck(seqFileName, seqIn, std::ios::in);

	DNASequence seq, seqRC;
 	CountedIntegralTuple tuple;
	//	tuple.tupleSize = tupleSize;
	//UNUSED// ssize_t nextPos;
	CountedIntegralTupleList::iterator listIt;
	CountedIntegralTuple nextTuple;
	ssize_t seqIndex = 0;
	DNASequence *seqs[2];
	seqs[0] = &seq;
	seqs[1] = &seqRC;
	while(SeqReader::GetSeq(seqIn, seq, SeqReader::noConvert)) {
		if (seqIndex % 100000 == 0 && seqIndex) {
			std::cout << seqIndex << std::endl;
		}
		++seqIndex;
		ssize_t curPos;

		if (seq.length < tupleSize)
			continue;
		
		
		if (skipGapped) {
			ssize_t seqIsGapped = 0;
			for (curPos = 0; curPos < seq.length; curPos++ ){ 
				if (unmasked_nuc_index[seq.seq[curPos]] >= 4) {
					seqIsGapped = 1;
					break;
				}
			}
			if (seqIsGapped) {
				continue;
			}
		}
				
		MakeRC(seq, seqRC);
		ssize_t s;
		DNASequence *seqPtr;
		for (s = 0; s < 2; s++ ) {
			seqPtr = seqs[s];
			if (ValidTuple(seqPtr->seq, 0, seqPtr->length, IntegralTuple::tupleSize)) {
				tuple.StringToTuple(&seqPtr->seq[0]);
				curPos = 0;
			}
			else {
				// find the last invalid nuc in the first tupleSize positions.
				ssize_t lastInvalidPos = 0;
				ssize_t p;
				for (p = 0; p < IntegralTuple::tupleSize; p++) {
					if (unmasked_nuc_index[seqPtr->seq[p]] >= 4)
						lastInvalidPos = p;
				}
				AdvanceToValid(seqPtr->seq, lastInvalidPos, 
											 seqPtr->length, IntegralTuple::tupleSize, tuple, curPos);
			}

			while (curPos + IntegralTuple::tupleSize <= seqPtr->length) {
				// assert that the next position
				// Better find this tuple			
				// 
				listIt = std::lower_bound(tupleList.begin(), tupleList.end(), tuple);
				assert(listIt != tupleList.end());
				(*listIt).IncrementMult();
			
				// the next one to store is from cur pos + 1 to curPos + tupleSize
				// make sure this position is ok to store
				if (curPos + IntegralTuple::tupleSize == seqPtr->length)
					break;

				if (unmasked_nuc_index[seqPtr->seq[curPos + IntegralTuple::tupleSize]] >= 4) {
					AdvanceToValid(seqPtr->seq, curPos + IntegralTuple::tupleSize,
												 seqPtr->length, IntegralTuple::tupleSize, tuple, curPos);
				}
				else {
					
					ForwardNuc(tuple,  unmasked_nuc_index[seqPtr->seq[curPos + IntegralTuple::tupleSize]], nextTuple);
				}
				++curPos;
				tuple.CopyTuple(nextTuple);
				//				tuple.tuple = nextTuple.tuple;
			}
		}
	}
	// There was an extra count for each tuple (due to the constructor), get rid of that.
	for (i = 0; i < nTuples; i++ ) 
		tupleList[i].count--;

	listIn.close();
	std::ofstream listOut;
	openck(tupleListName, listOut,std::ios::out | std::ios::binary);
	tupleSize_SSZT = tupleSize;
	listOut.write((char*) &tupleSize_SSZT, sizeof(tupleSize_SSZT));
	listOut.write((char*) &nTuples, sizeof(ssize_t));
	listOut.write((char*) &tupleList[0], sizeof(CountedIntegralTuple)*nTuples);

	return 0;
}

ssize_t ValidTuple(unsigned char* seq, ssize_t pos, ssize_t seqLength, ssize_t tupleLength) {
	ssize_t p;
	if (pos + tupleLength > seqLength)
		return 0;

	for (p = pos; p < pos + tupleLength ; p++) {
		if (numeric_nuc_index[seq[p]] >= 4)
			return 0;
	}

	return 1;
}

void AdvanceToValid(unsigned char *seq, ssize_t curPos, ssize_t seqLength, ssize_t tupleLength,
										CountedIntegralTuple &tuple, ssize_t &nextPos) {
	// assert that curPos is invalid.
	ssize_t pi;
	ssize_t nextInvalidPos = curPos;
	ssize_t p = curPos + 1;
	//UNUSED// ssize_t validFound = 0;
	do {
		p = nextInvalidPos + 1;
		// If there is no room for the next valid pos, bail out.
		if (p + tupleLength > seqLength) {
			nextPos = seqLength + 1;
			return;
		}
		
		for( pi = p; pi < p + tupleLength; pi++ ) {
			if (numeric_nuc_index[seq[pi]] >= 4) {
				nextInvalidPos = pi;
			}
		}
	} while (p < nextInvalidPos and p + tupleLength <= seqLength);
	if (p + tupleLength >= seqLength) {
		nextPos = seqLength + 1;
	}
	else {
		nextPos = p;
	}
}

