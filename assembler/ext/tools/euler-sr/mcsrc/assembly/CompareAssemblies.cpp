/***************************************************************************
 * Title:          CompareAssemblies.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "utils.h"
#include "Spectrum.h"
#include "DeBruijnGraph.h"
#include "HashedSpectrum.h"

#include "IntegralTupleStatic.h"

#include <vector>
#include <iostream>

class Enum {
public:
	ssize_t pos;
	ssize_t strand;
	ssize_t index;
	ssize_t queryPos;
};
bool operator<(const Enum &a, const Enum &b) {
	if (a.strand != b.strand) return a.strand < b.strand;
	else return a.pos < b.pos;
}

class OrderByQuery {
	public:
	ssize_t operator()(Enum a, Enum b) {
		return (a.queryPos < b.queryPos);
	}
};

int main(int argc, char* argv[]) {
	ssize_t minStripLength = 30;
	std::string refName, contigsName;
	if (argc < 3) {
		std::cout << "usage: cmpAssemblies refFile contigsFile " << std::endl;
		std::cout << "    -minLength minLength" << std::endl;
		std::cout << "    -onlyBad " << std::endl;
		std::cout << "    -maskLower - Mask off lower case letters." << std::endl;
		std::cout << "    -tupleSize t Use 't' as a tuple size." << std::endl;
		std::cout << "    -printPos  - Print the position of a map. " << std::endl;
		exit(1);
	}
	refName = argv[1];
	contigsName = argv[2];
	int tupleSize = 20;
	ssize_t minLength = 0;
	int argi = 3;
	ssize_t printOnlyBad = 0;
	ssize_t maskLower = 0;
	ssize_t printPos = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-minLength") == 0) {
			minLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-printPos") == 0){ 
			printPos = 1;
		}
		else if (strcmp(argv[argi], "-onlyBad") == 0){ 
			printOnlyBad = 1;
		}
		else if (strcmp(argv[argi], "-maskLower") == 0) {
			maskLower = 1;
		}
		else if (strcmp(argv[argi], "-tupleSize") == 0) {
			tupleSize = atoi(argv[++argi]);
		}
		else {
			std::cout << "bad option: " << argv[argi]<< std::endl;
			exit(1);
		}
		++argi;
	}

	HashValueFunctor calcHashValue;
	calcHashValue.hashLength = 10;
	std::string contigPosName = contigsName + ".pos";

	//UNUSED// ssize_t s;
	//UNUSED// ssize_t storeSeq;
	//UNUSED// ssize_t seqNumber = 0;

	DNASequence ref, refRC;
	std::ifstream refIn;
	openck(refName, refIn, std::ios::in);
	if (maskLower) 
		SeqReader::MaskRepeats();
	SeqReader::GetSeq(refIn, ref, SeqReader::noConvert);
	
	MakeRC(ref, refRC);
	SimpleSequenceList refStrands;
	refStrands.resize(2);
	refStrands[0].seq = ref.seq;
	refStrands[1].seq = refRC.seq;
	refStrands[0].length = refStrands[1].length = ref.length;
	std::vector<CountedReadPos> refPositions;
	ssize_t numUnmasked = 0;
	ssize_t i, m;
	m = -1;
	for (i = 0; i < tupleSize; i++) {
		if (!IsUnmasked(ref.seq[i])) {
			m = i;
		}
	}
	for (i = 0; i < ref.length - tupleSize; i++) {
		// if there is a masked nucleotide that covers
		// part of the current tuple, record that.
		if (!IsUnmasked(ref.seq[i+tupleSize-1])) {
			m = i + tupleSize - 1;
		}
		// Record an unmasked tuple if none of it is masked.
		if (i > m) {
			numUnmasked++;
		}
	}
	refPositions.resize(numUnmasked*2);
	m = -1;
	for (i = 0; i < tupleSize; i++) {
		if (!IsUnmasked(ref.seq[i])) {
			m = i;
		}
	}
	ssize_t curTuple = 0;
	for (i = 0; i < ref.length - tupleSize + 1; i++) {
		if (!IsUnmasked(ref.seq[i+tupleSize-1])) {
			m = i+tupleSize-1;
		}
		if (i > m) {
			refPositions[curTuple].read = 0;
			refPositions[curTuple].pos = i;
			refPositions[curTuple + numUnmasked + 1].read = 1;
			refPositions[curTuple + numUnmasked + 1].pos = i;
			++curTuple;
		}
	}
	
  CompareTuples<SimpleSequenceList> comp;
  comp.sequencesPtr = &refStrands;
  comp.length = tupleSize;

  std::sort(refPositions.begin(), refPositions.end(), comp);

	ssize_t prev = 0;
	ssize_t cur = 1;
	ssize_t count = 1;
	while (cur < refPositions.size()) {
		while (cur < refPositions.size() and
					 comp(refPositions[prev],refPositions[cur]) == 0) cur++;
		count = cur - prev;
		for (i = prev; i < cur; i++) refPositions[i].count = count;
		prev = cur;
	}

	SimpleSequenceList contigs;
	ReadSimpleSequences(contigsName, contigs);

	CountedReadPos readPos;
	calcHashValue.sequences = &contigs;
	ssize_t c, p;
	unsigned char *tuple;
	ssize_t index;
	ssize_t nMapped;
	std::vector<Enum> enumeration;
	OrderByQuery orderByQuery;
	for (c = 0; c < contigs.size(); c++ ){
		if (contigs[c].length < minLength)
			continue;

		/*
		enumeration.resize(contigs[c].length-tupleSize+1);
		for (i = 0; i < contigs[c].length -tupleSize+1; i++ ) {
			enumeration[i].strand = -1;
			enumeration[i].pos    = -1;
			enumeration[i].index  = i;
		}
		*/
		nMapped = 0;
		Enum  en;
		enumeration.clear();
		for (p = 0; p < contigs[c].length - tupleSize + 1; p++) { 
			tuple = &contigs[c].seq[p];
			if ((index = LocateTuple(refStrands, refPositions, tupleSize, (char*) tuple)) >= 0) {
				if (refPositions[index].count == 1) {
					/*		enumeration[p].strand = refPositions[index].read;
								enumeration[p].pos = refPositions[index].pos;
					*/
					en.strand = refPositions[index].read;
					en.pos = refPositions[index].pos;
					en.index = nMapped;
					en.queryPos = p;
					enumeration.push_back(en);
					nMapped++;
				}
			}
		}
		std::sort(enumeration.begin(), enumeration.end());
		ssize_t nStrips = 1;
		ssize_t stripStart = 0;
		std::vector<ssize_t> stripLengths, stripStarts, stripEnds, stripStrands;
		std::vector<ssize_t> qryStart, qryEnd;
		p = 0;

		// Compress small strips
		ssize_t prevNStrips;
		ssize_t curNStrips;

		curNStrips = enumeration.size();
		prevNStrips = 0;
		ssize_t iter = 0;

		while (curNStrips != prevNStrips) {
			p = 0;
			ssize_t numSmallRemoved = 0;
			while (p < enumeration.size() and enumeration[p].pos == -1) p++;
			stripStart = p;
			ssize_t stripLength;
			ssize_t r;
			for (p = p + 1; p < enumeration.size(); p++) {
				if (enumeration[p].index != enumeration[p-1].index + 1 or 
						enumeration[p].strand != enumeration[p-1].strand) {
					stripLength = p - stripStart;
					if (stripLength < minStripLength) {
						for (r = stripStart; r < p; r++) {
							enumeration[r].index = -1;
							enumeration[r].pos   = -1;
							enumeration[r].strand = -1;
							enumeration[r].queryPos = -1;
							numSmallRemoved++;
						}
					}
					stripStart = p;
				}
			}
			// Mask the last strip if it is small.
			stripLength = p - stripStart;
			if (stripLength < minStripLength) {
				for (r = stripStart; r < p; r++) {
					enumeration[r].index = -1;
					enumeration[r].pos   = -1;
					enumeration[r].strand = -1;
					enumeration[r].queryPos = -1;
					numSmallRemoved++;
				}
			}

			std::sort(enumeration.begin(), enumeration.end(), orderByQuery);
			p = 0;
			ssize_t cur = 0;
			
			while (p < enumeration.size() and enumeration[p].index == -1) p++;
			while (p < enumeration.size()) {
				enumeration[cur] = enumeration[p];
				enumeration[cur].index = cur;
				cur++;
				p++;
			}
			enumeration.resize(cur);
			std::sort(enumeration.begin(),enumeration.end());
			prevNStrips = curNStrips;
			curNStrips = enumeration.size();
			
			//			std::cout << "iter: " << iter << " prev: " << prevNStrips << " cur: " << curNStrips << std::endl;
			iter++;
		}
			
		nStrips = 1;
		p = 0;
		while (p < enumeration.size() and enumeration[p].pos == -1) p++;
		stripStart = p;
		for (p = p + 1; p < enumeration.size(); p++) {
			if (enumeration[p].index != enumeration[p-1].index + 1 or 
					enumeration[p].strand != enumeration[p-1].strand) {
				stripLengths.push_back(p - stripStart);
				stripStarts.push_back(enumeration[stripStart].pos);
				stripEnds.push_back(enumeration[p-1].pos);
				stripStrands.push_back(enumeration[stripStart].strand);
				qryStart.push_back(enumeration[stripStart].queryPos);
				qryEnd.push_back(enumeration[p-1].queryPos);
				nStrips++;
				stripStart = p;
			}
		}
		if (enumeration.size() > 0) {
			stripLengths.push_back(p - stripStart);
			stripStarts.push_back(enumeration[stripStart].pos);
			stripEnds.push_back(enumeration[p-1].pos);
			stripStrands.push_back(enumeration[stripStart].strand);
			qryStart.push_back(enumeration[stripStart].queryPos);
			qryEnd.push_back(enumeration[p-1].queryPos);
		}
		if (printPos) {
			if (nStrips > 1) {
				if (stripStrands[0] == 0) {
					std::cout << stripStarts[0] << " " << stripEnds[0] << " " 
										<< stripEnds[0] - stripStarts[0] + 1 << " " << stripStrands[0] << std::endl;
				}
				else {
					std::cout << ref.length - stripEnds[0] << " " << ref.length - stripStarts[0]
										<< " " << stripEnds[0] - stripStarts[0] + 1 << " " << stripStrands[0] << std::endl;
				}
			}
		}
		else if ((printOnlyBad and nStrips > 1) or !printOnlyBad) {
			std::cout << "match len refStart refEnd refStrand qryStart qryEnd" << std::endl; 
			std::cout << "Mapped " << c << " " << contigs[c].length << " "
								<< nMapped  << " positions with " 
								<< nStrips << " strips" << std::endl;
			ssize_t s;
			for (s = 0; s < stripLengths.size(); s++ ) {
				std::cout << stripLengths[s] << "\t" << stripStarts[s] << "\t" << stripEnds[s] << "\t" 
									<< stripStrands[s] << "\t" << qryStart[s] << "\t" << qryEnd[s] << std::endl;
			}
			for (p = 1; p < enumeration.size(); p++) {
				ssize_t qryLength, refLength;
				if (enumeration[p].index == enumeration[p-1].index + 1 and
						enumeration[p].strand == enumeration[p-1].strand) {
					// In the same strip, compute distance between anchors
					refLength = enumeration[p].pos - enumeration[p-1].pos;
					qryLength = enumeration[p].queryPos - enumeration[p-1].queryPos;
					if (abs(refLength - qryLength) > 100) {
						std::cout << "gap of length: " << refLength - qryLength << std::endl;
					}
				}
			}
		}
	}
}

