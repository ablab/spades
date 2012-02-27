/***************************************************************************
 * Title:          SimulateMatePairsFromReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqReader.h"
#include "utils.h"
#include "ParseTitle.h"
#include "SeqUtils.h"
#include "SimpleStats.h"
#include <map>
#include <sstream>

class ReadMap {
public:
	ssize_t start, end;
	ssize_t strand;
	std::string title;
};


class SortReadMap {
public:
	ssize_t operator()(const ReadMap &a, const ReadMap &b) {
		return a.start < b.start;
	}
};
ssize_t PickRead(std::vector<ReadMap> &readMapList, ssize_t startPos, ssize_t length);
void PrintUsage() {
	std::cout << "usage: simMPR refSeq readsFile mapFile outputFile " << std::endl;
	std::cout << "       [readSize, cloneSize, cloneStdev, count]" << std::endl;
}

int main(int argc, char *argv[]) {
	std::string readsFile, posFile, refSeqFile, mapFile, outputFile;
	if (argc < 5) {
		PrintUsage();
		exit(1);
	}
	refSeqFile= argv[1];
	readsFile = argv[2];
	mapFile   = argv[3];
	outputFile= argv[4];

	DNASequence refSeq;
	SeqReader::GetSeq(refSeqFile, refSeq, SeqReader::noConvert);
	std::cout << "ref seq length: " << refSeq.length << std::endl;

	std::vector<ssize_t> readLengths, cloneLengths, cloneStddevs, counts;
	int argi = 5;
	while (argi < argc) {
		if (argi + 4 > argc) {
			PrintUsage();
			exit(1);
		}
		readLengths.push_back(atoi(argv[argi++]));
		cloneLengths.push_back(atoi(argv[argi++]));
		cloneStddevs.push_back(atoi(argv[argi++]));
		counts.push_back(atoi(argv[argi++]));
	}
	
	if (readLengths.size() < 1) {
		std::cout << "Usage: simMatePairsFromReads readsFile posFile outputFile "
							<< "[readsize clonesize clonestddev len]"
							<< std::endl;
		exit(1);
	}
	std::ofstream matesOut;
	openck(outputFile, matesOut, std::ios::out);
	
	std::ifstream readsIn, mapIn;
	openck(readsFile, readsIn, std::ios::in);
	openck(mapFile, mapIn, std::ios::in);

	DNASequenceList reads;
	ReadDNASequences(readsFile, reads);
	//	std::cout << "got reference of " << reads.size()  << std::endl;
	ReadMap readMap;
	std::vector<ReadMap> readMapList;
	while(mapIn) {
		mapIn >> readMap.start >> readMap.end >> readMap.title >> readMap.strand;
		readMapList.push_back(readMap);
	}

	SortReadMap sortReadMap;
	std::sort(readMapList.begin(), readMapList.end(), sortReadMap);

	std::map<std::string, ssize_t> readNameToIndex;
	ssize_t r;
	std::string name;

	// Create the map from fasta names to index
	for (r = 0; r < reads.size(); r++ ) {
		ParseTitle(reads[r].namestr, name);
		readNameToIndex[name] = r;
	}

	std::ofstream ruleOut, pairFile;
	
	std::string nameRuleName("name.rul");
	std::string pairFileName = outputFile + ".pair";

	openck(nameRuleName, ruleOut, std::ios::out);
	openck(pairFileName, pairFile, std::ios::out);

  char cloneKeys[12][2] = {{'a', 'b'}, {'c', 'd'}, {'e', 'f'}, {'g', 'h'}, 
													 {'i', 'j'}, {'k', 'l'}, {'m', 'n'}, {'o', 'p'},
													 {'q', 'r'}, {'u', 'v'}, {'w', 'x'}, {'y', 'z'}};
	

	ssize_t i;
	ruleOut<<"		/*			Names are: Library_Plate_well.primer	" << std::endl;
	ruleOut<<"				primer  library plate   length range	*/" << std::endl;
	ruleOut<<""<< std::endl;
	ruleOut<<"Single reads:                   s       ALL     ALL" << std::endl;
 
	for (i = 0; i < cloneLengths.size(); i++) {
		ruleOut<<"Double-barreled reads:          " 
					 << cloneKeys[i][0] <<" " << cloneKeys[i][1]<<"     ALL     ALL     "
					 << cloneLengths[i] + readLengths[i]*2 - cloneStddevs[i]<< " " 
					 << cloneLengths[i] + readLengths[i]*2 + cloneStddevs[i]<< std::endl;
	}
	ruleOut.close();

	SeedRandom();

	
	ssize_t mateLib;
	ssize_t read;
	//UNUSED// ssize_t cloneStart;
	ssize_t mateLength, readLength, mateDev;
	ssize_t leftStart, rightStart;
	ssize_t leftReadIndex = -1, rightReadIndex = -1;
	//UNUSED// ssize_t leftOffset, rightOffset;
	std::stringstream namestrm;
	ssize_t readIndex = 0;
	//UNUSED+// ssize_t rightDNAIndex;
	ssize_t leftDNAIndex ;
	ssize_t leftMapIndex, rightMapIndex;
	DNASequence readRC;
	DNASequence leftRead, rightRead;
	leftRead._ascii = 1;
	rightRead._ascii = 1;
	readRC._ascii = 1;
	ssize_t randMateLength;
	ssize_t numTries= 0;
	ssize_t cloneIndex = 0;
	for (mateLib = 0; mateLib < readLengths.size(); mateLib++) {
		mateLength = cloneLengths[mateLib];
		readLength = readLengths[mateLib];
		mateDev    = cloneStddevs[mateLib];
		for (read = 0; read < counts[mateLib] ; read++) {
			randMateLength = readLength + (ssize_t) NormalRandom(double(mateLength), double(mateDev));

			leftMapIndex = -1; rightMapIndex = -1; numTries = 0;
			while(leftMapIndex == -1 or rightMapIndex == -1) {
				if (numTries > 100) {
					std::cout << "WARNING, tried to find a valid read start 100 times but failed."<<std::endl;
					std::cout << "something is probably wrong with the input." << std::endl;
					exit(1);
				}
				leftStart    = Random(refSeq.length - randMateLength - 2 * readLength);
				leftMapIndex = PickRead(readMapList, leftStart, readLength);
				/*
					std::cout << "picked: " << leftMapIndex << " " << readMapList[leftMapIndex].start << " " 
					<< reads[leftMapIndex].length  << " " << leftStart << std::endl;
				*/
				rightStart = leftStart + randMateLength + readLength;
				rightMapIndex = PickRead(readMapList, rightStart, readLength);
				++numTries;
			}

			leftDNAIndex = -1;
			assert (leftMapIndex != -1 and rightMapIndex != -1);
			ssize_t readPos;
			if (readNameToIndex.find(readMapList[leftMapIndex].title) != readNameToIndex.end()) {

				// maybe the map was off a little bit with the read ends (indels and such)

				leftReadIndex = readNameToIndex[readMapList[leftMapIndex].title];
				readPos = leftStart - readMapList[leftMapIndex].start;				
				
				if (readPos < 0) 
					readPos = 0;

				else if ( readPos + readLength >= reads[leftReadIndex].length) 
					readPos = reads[leftReadIndex].length - readLength;


				leftRead.seq = &reads[leftReadIndex].seq[readPos];
				leftRead.length = readLength;
				namestrm.str("");
				namestrm << "sim_1_" << cloneIndex << "_" << cloneKeys[mateLib][0];
				leftRead.namestr = namestrm.str();
				if (readMapList[leftMapIndex].strand == 0) {

					leftRead.PrintSeq(matesOut);
				}
				else {
					MakeRC(leftRead, readRC);
					readRC.namestr= leftRead.namestr;
					readRC.PrintSeq(matesOut);
					readRC.Reset();
				}
				matesOut << std::endl;
			}
			readIndex++;
			rightReadIndex = -1;
			if (readNameToIndex.find(readMapList[rightMapIndex].title) != readNameToIndex.end()) {
				rightReadIndex = readNameToIndex[readMapList[rightMapIndex].title];

				// maybe the map was off a little bit with the read ends (indels and such)
				readPos = rightStart - readMapList[rightMapIndex].start;				
				
				if (readPos < 0) 
					readPos = 0;
				else if ( readPos + readLength >= reads[rightReadIndex].length) 
					readPos = reads[rightReadIndex].length - readLength;

				rightRead.seq = &reads[rightReadIndex].seq[readPos];
				rightRead.length = readLength;
				namestrm.str("");
				namestrm << "sim_1_" << cloneIndex << "_" << cloneKeys[mateLib][1];
				rightRead.namestr = namestrm.str();
				if (readMapList[rightMapIndex].strand == 1) {
					rightRead.PrintSeq(matesOut);
				}
				else {
					MakeRC(rightRead, readRC);
					readRC.namestr= rightRead.namestr;
					readRC.PrintSeq(matesOut);
					readRC.Reset();
				}
				matesOut << std::endl;
			}
			readIndex++;
			if (leftReadIndex != -1 and rightReadIndex != -1) {
					pairFile << "sim_1_" << cloneIndex << " " << readIndex-2
									 << " " << readIndex-1 << " " << mateLib << std::endl;
      }
			cloneIndex++;
		}
	}
}

ssize_t PickRead(std::vector<ReadMap> &readMapList, ssize_t startPos, ssize_t readLength) {

	ssize_t begin, end, cur;
	begin = 0; end = readMapList.size();
	cur = (end + begin) / 2;
	// Find a read that covers 'startPos'
	while (cur < end and 
				 readMapList[cur].start !=  startPos) {
		if (cur < end and 
				readMapList[cur].start < startPos and 
				readMapList[cur].end > startPos+readLength)
			break;

		if (startPos < readMapList[cur].start) {
			end = cur;
		}
		else if (startPos > readMapList[cur].start) {
			if (cur == begin)
				begin = cur+ 1;
			else
				begin = cur;
		}			
		cur = (end + begin)/2;
	}
	if (cur >= end or 
			readMapList[cur].start > startPos  or
			readMapList[cur].end <= startPos + readLength) {
		// no valid start found.
		return -1;
	}

	while (cur > 0 and readMapList[cur-1].start < startPos and 
				 readMapList[cur-1].end > startPos + readLength)
		cur--;

	ssize_t last = cur;
	while (last < readMapList.size()-1 and 
				 startPos > readMapList[last].start and
				 startPos + readLength < readMapList[last].end) {
		last++;
	}
	//	std::cout << "got range: " << cur << " " << last << std::endl;
	ssize_t length = last - cur;
	ssize_t readIndex = Random(length) + cur;
	assert(readMapList[cur].start <= startPos and 
				 readMapList[cur].end   >= startPos + readLength);
	assert(readMapList[readIndex].start <= startPos and 
				 readMapList[readIndex].end > startPos + readLength);
	//	std::cout << "got read index: " << readIndex <<  " (" << length << " " << readIndex - cur << ")" << std::endl;
	return readIndex;
}

