/***************************************************************************
 * Title:          AlignedReadExtender.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqUtils.h"
#include "SeqReader.h"
#include "DNASequence.h"
#include "utils.h"
#include "SimpleStats.h"
#include "IntegralTupleStatic.h"

#include <string>
#include <vector>

using namespace std;

class Map {
public:
	ssize_t pos;
	ssize_t dist;
	ssize_t strand;
	ssize_t readIndex;
	ssize_t index;
	ssize_t mateIndex;
};

class SortByReadIndex {
public:
	ssize_t operator()(const Map&lhs, const Map &rhs) {
		if (lhs.readIndex < rhs.readIndex) return 1;
		else return 0;
	}
};
class SortByPos {
public:
	ssize_t operator()(const Map&lhs, const Map &rhs) {
		if (lhs.pos < rhs.pos) return 1;
		else return 0;
	}
};

class SortByMapIndex {
public:
	ssize_t operator()(const Map &lhs, const Map &rhs) {
		if (lhs.index < rhs.index) return 1;
		else return 0;
	}
};

void PrintUsage() {
	cout << "usage: alignedReadExtender genome reads map errorProfile readLength coverage mateSpan readsOut" << endl;
	cout << "   -maxDist d  Do not extend reads if they are more than 'd' away from the reference seq." <<endl;

}

int main(int argc, char* argv[]) {

	string genomeFileName, readsFileName, mapFileName, errorProfileName, readsOutName;
	if (argc < 9) {
		PrintUsage();
		exit(1);
	}
	int argi = 0;
	genomeFileName = argv[++argi];
	readsFileName  = argv[++argi];
	mapFileName    = argv[++argi];
	errorProfileName=argv[++argi];
	ssize_t readLength = atoi(argv[++argi]);
	ssize_t coverage = atoi(argv[++argi]);
	ssize_t mateSpan = atoi(argv[++argi]);
	readsOutName = argv[++argi];
	ssize_t maxDist = -1;
	++argi;
	while (argi < argc) {
		if (strcmp(argv[argi], "-maxDist") == 0) {
			maxDist = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			cout << "Bad option: " << argv[argi] << endl;
			exit(1);
		}
		++argi;
	}

	//
	// Open all necessary files so that we can bail early
	// if one doesn't exist.
	
	ifstream mapIn, profileIn, readsIn;
	openck(mapFileName, mapIn, ios::in);
	openck(errorProfileName, profileIn, std::ios::in);
	openck(readsFileName, readsIn, std::ios::in);
	ofstream readsOut;
	openck(readsOutName, readsOut, std::ios::out);

	DNASequence genome;
	SeqReader::GetSeq(genomeFileName, genome, SeqReader::noConvert);

	//
	// Count the lines in the map file
	//

	string line;
	ssize_t numLines = 0;
	while (getline(mapIn, line)) {
		++numLines;
	}
	// Reset the map file.
	mapIn.close(); mapIn.clear();

	//
	// Get the error profile (the number of bases to
	// extend each read).
	//
	vector<double> errorProfile;
	while(profileIn){ 
		double errProf;
		if (!(profileIn >> errProf))
			break;
		errorProfile.push_back(errProf);
	}
	profileIn.close();

	//
	// Read the map file.
	//
	openck(mapFileName, mapIn, ios::in);
	Map map;
	vector<Map> mappedReads;
	ssize_t error = 0;
	ssize_t mapIndex = 0;
	while(mapIn) {
		if (!(mapIn >> map.pos)) { error = 1; break; }
		if (!(mapIn >> map.dist)) { error = 2; break; }
		if (!(mapIn >> map.strand)) { error = 3; break; }
		if (!(mapIn >> map.readIndex)) { error = 4; break; }
		map.index = mapIndex;
		map.mateIndex = -1;
		++mapIndex;
		mappedReads.push_back(map);
	}

	// 
	// Assign mate pairs.
	// 
	cout << "Read: " << mappedReads.size() << " reads." << endl;
	sort(mappedReads.begin(), mappedReads.end(), SortByPos());

	double pSample;
	if ((double(genome.length) * double(coverage) / double(readLength)) > mappedReads.size()) {
		cout << "WARNING: There are only enough mapped reads for a coverage of " << 
			(double(mappedReads.size()) * readLength) / genome.length  << ". This will "
				 << "be less than the specified coverage of " << coverage << endl;
		pSample = 1.1;
	}
	else {
		pSample = (((double(genome.length) * coverage) / readLength ) /  mappedReads.size());
		cout << "Sampling at a rate of " << pSample << endl;
	}

	// Sample reads to reach the approximate coverage.

	ssize_t i;
	ssize_t numMapped = mappedReads.size();
	ssize_t numNoPairFound = 0;
	ssize_t numPairFound = 0;
	mateSpan += readLength + 2*errorProfile.size();
	for (i = 0; i < numMapped; i++ ){ 
		// Skip reads if the alignment distance is too large.
		if (maxDist >= 0 and mappedReads[i].dist > maxDist)
			continue;

		// Fix accidental map problems
		if (mappedReads[i].pos < 0)
			continue;
		if (Uniform() > pSample) {
			continue;
		}
		// Don't bother with this read if it is already included as 
		// a mate pair.
		if (mappedReads[i].mateIndex > 0) {
			continue;
		}

		// Search for the first not-taken mate for this read.
		ssize_t m = i + 1;

		// Advance until the appropriate distance is found.
		while (m < numMapped and 
					 mappedReads[m].pos - mappedReads[i].pos < mateSpan*0.90) 
			m++;

		// Advance while mate-reads are already found, or they are too divergent.
		while (m < numMapped and
					 (mappedReads[m].mateIndex != -1 or
						(maxDist >= 0 and mappedReads[m].dist > maxDist)))
			m++;
		
		if (m >= numMapped or 
				mappedReads[m].pos - mappedReads[i].pos > mateSpan*2) {
			numNoPairFound++;
			continue;
		}

		// Assign mate-pairs;
		mappedReads[m].mateIndex = mappedReads[i].index;
		mappedReads[i].mateIndex = mappedReads[m].index;
		numPairFound++;
	}
	
	cout << "found " << numPairFound << " pairs and " << numNoPairFound << " unpaired." << endl;

	// Sort the readsy by their position in the map file.  This way 
	// they are also in order of reads.  The mapped reads are 
	// a strict subset of the reads.

	sort(mappedReads.begin(), mappedReads.end(), SortByMapIndex());
	cout << "done re-ordering." << endl;
	//
	// Now print the reads that are marked by a mate.
	//
	ssize_t readIndex = 0;

	DNASequence read, mateRC;

	string seq;
	string readStr;
	readIndex = -1;
	mapIndex = -1;
	DNASequence beginMate, endMate;
	DNASequence readRC;

	string defaultTitle = "SLXA_SIM_0:0:0:0:";

	while (mapIndex < numMapped) {
		++mapIndex;
		if (mapIndex % 10000 == 0) {
			cout << mapIndex << endl;
		}
		while (mappedReads[mapIndex].mateIndex == -1) mapIndex++;
		
		if (mapIndex >= numMapped)
			break;

		while (readIndex != mappedReads[mapIndex].readIndex) {
			readIndex++;
			SeqReader::GetSeq(readsIn, read, SeqReader::noConvert);
		}
		
		if (mappedReads[mapIndex].mateIndex > mapIndex) {
			// Have the correct read and sequence here.
			// append to the read the genome sequence.
			
			if (mappedReads[mapIndex].strand == 0) {
				readStr.assign((char*) read.seq, read.length);
			}
			else {
				readRC.Reset(0);
				MakeRC(read, readRC);
				readStr.assign((char*) readRC.seq, readRC.length);
			}
			readStr.append((const char*) &genome.seq[mappedReads[mapIndex].pos + read.length], 
										 errorProfile.size());
				
			// Now corrupt the read.
			for (i = 0; i < errorProfile.size(); i++ ){ 
				if (Uniform() < errorProfile[i]) {
					readStr[read.length + i] = tolower(MutateNuc(readStr[read.length + i]));
				}
			}
			
			// Output the read.
			stringstream titleStrm;
			titleStrm << defaultTitle << mappedReads[mapIndex].index << "/1";
			
			readsOut << ">" << titleStrm.str() << endl;
			readsOut << readStr << endl;
		}
		else {
			// 
			// The mate index is before the mapIndex (current mapped read).
			// That means this read should be output in the reverse
			// orientation, and should use the index of the mate in the title.
			//
			//UNUSED// ssize_t mateIndex = mappedReads[mapIndex].mateIndex;
			if (mappedReads[mapIndex].strand == 0) 
				readStr.assign((char*) read.seq, read.length);
			else {
				// Make sure the read is being extended in the forward direction.
				readRC.Reset(0);
				MakeRC(read, readRC);
				readStr.assign((char*) readRC.seq, readRC.length);
			}
			readStr.append((const char*) &genome.seq[mappedReads[mapIndex].pos + read.length], 
										 errorProfile.size());

			DNASequence mateRead;
			mateRead.seq = (unsigned char*) readStr.c_str();
			mateRead.length = read.length + errorProfile.size();
			MakeRC(mateRead, readRC);
			// Now corrupt the read.
			for (i = 0; i < errorProfile.size(); i++ ){ 
				if (Uniform() < errorProfile[i]) {
					readRC.seq[read.length + i] = tolower(MutateNuc(readRC.seq[read.length + i]));
				}
			}
			stringstream titlestrm;
			titlestrm << defaultTitle << mappedReads[mapIndex].mateIndex << "/2";
			readRC.namestr = titlestrm.str();
			readRC.PrintlnSeq(readsOut);
		}
	}
	return 0;
}
