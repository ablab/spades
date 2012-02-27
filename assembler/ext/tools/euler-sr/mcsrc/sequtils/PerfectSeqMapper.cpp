/***************************************************************************
 * Title:          PerfectSeqMapper.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/19/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "compatibility.h"

class MatchTuple {
public:
	char *tuple;

	ssize_t index;

	// HAS BITFIELD
	//	unsigned int pos:31;
	//	unsigned int strand:1;
	_UINT_ pos:(_INT_BITS_-1);
	bool strand:1;

	static int tupleSize;
	MatchTuple() {
		//		if (tupleSize != 0) {
		//			tuple = new char[tupleSize];
			//			std::cout << "allocaing tuple " << (unsigned long) tuple << std::endl;
			//		}
			//		else {
			//			tuple = NULL;
			//		}
	}
	void Store(char *seq, ssize_t seqPos, ssize_t seqStrand, ssize_t seqIndex) {
		//UNUSED// ssize_t i;
		//		std::cout << "storing " << seqPos << " " << seqStrand << " " << seqIndex << std::endl;
		/*
		for (i = 0; i < tupleSize; i++) {
			
			tuple[i] = tolower(seq[seqPos + i]);
			}*/
		tuple = &seq[seqPos];
		pos = seqPos;
		strand = seqStrand;
		index = seqIndex;
	}
	bool operator<(const MatchTuple &cmp) const {
		return strncmp(tuple, cmp.tuple, tupleSize);
	}
	bool LessThan(const char *tuplePtr) const {
		return  strncmp(tuple, tuplePtr, tupleSize) < 0;
	}
	bool Equals(const char *tuplePtr) const {
		return strncmp(tuple, tuplePtr, tupleSize) == 0;
	}
};

class CompareMatchTuples {
public:
	ssize_t operator()(const MatchTuple &a, const MatchTuple & b) {
		return strncmp(a.tuple, b.tuple, a.tupleSize) < 0;
	}
};

void First(std::vector<MatchTuple> &tupleList, char *tuple, ssize_t &index) {
	while (index > 0 and tupleList[index-1].Equals(tuple))
		index--;
}

void Last(std::vector<MatchTuple> &tupleList, char*tuple, ssize_t &index) {
	while(index < tupleList.size() and tupleList[index].Equals(tuple))
		index++;
}

ssize_t FindTuple(std::vector<MatchTuple> &tupleList, char* tuple, ssize_t &index) {
	ssize_t begin, end, cur;
	begin = 0; end = tupleList.size();
	cur = (begin+ end) / 2;
	while (cur < end) {
		if (tupleList[cur].Equals(tuple)) {
			index = cur;
			return 1;
		}
		if (tupleList[cur].LessThan(tuple)) {
			begin = cur + 1;
			cur = (begin + end) / 2;
		}
		else {
			end = cur;
			cur = (begin + end) / 2;
		}
	}
	if (cur < tupleList.size() and
			tupleList[cur].Equals(tuple)) {
		index = cur;
		return 1;
	}
	return 0;
}

int MatchTuple::tupleSize = 0;
int main(int argc, char *argv[]) {

	std::string refSeqFile, seqListFile;
	int tupleSize;
	if (argc < 4) {
		std::cout <<" usage: mapper refSeq seqList tupleSize" << std::endl;
		exit(1);
	}

	refSeqFile = argv[1];
	seqListFile = argv[2];
	tupleSize = MatchTuple::tupleSize = atoi(argv[3]);

	DNASequence refSeq;
	//	std::cout << "opening " << refSeqFile << std::endl;
	SeqReader::GetSeq(refSeqFile, refSeq);
	
	//	std::cout << "storing pos list" << std::endl;
	std::vector<MatchTuple> refTupleList;
	//	refTupleList.resize((refSeq.length - tupleSize + 1)* 2 );
	ssize_t i;
	for (i = 0; i < refSeq.length - tupleSize + 1; i++) {
		//		refTupleList[i] = new MatchTuple;
		refTupleList.push_back(MatchTuple());
		refTupleList[i].Store((char*)refSeq.seq, i, 0, i);
	}
	DNASequence refSeqRC;
	MakeRC(refSeq, refSeqRC);
	ssize_t cur = i;
	for (i = 0; i < refSeq.length - tupleSize + 1; i++) {
		//		refTupleList[cur] = new MatchTuple;
		refTupleList.push_back(MatchTuple());
		refTupleList[cur].Store((char*) refSeqRC.seq, i, 1, cur);
		cur++;
	}
	//	std::cout << cur << " " << refTupleList.size() << " sorting " << std::endl;
	std::sort(refTupleList.begin(), refTupleList.end(), CompareMatchTuples());
		
	//	std::cout << "done." << std::endl;

	// Compress the tuple list;
	/*
	cur = 0;
	for (i = 1; i < refTupleList.size(); i++) {
		if (!refTupleList[i].Equals(refTupleList[cur].tuple)) {
			cur++;
			refTupleList[cur] = refTupleList[i];
		}
	}
	std::cout << "removed " << refTupleList.size() - cur << " non-unique tuples " << std::endl;
	refTupleList.resize(cur);
	*/

	DNASequenceList sequences;
	ReadDNASequences(seqListFile, sequences);
	ssize_t s, p;
	// Now map the sequences.
	ssize_t index;
	ssize_t first, last;
	ssize_t tupIndex;
	ssize_t isMatch;
	char *sbjctSeq;
	for (s = 0; s < sequences.size(); s++) { 
		if (sequences[s].length  >= tupleSize) {
			if (FindTuple(refTupleList, (char*) sequences[s].seq, index)) {
				first = index;
				last = index+1;
				assert(last <= refTupleList.size());
				First(refTupleList, (char*) sequences[s].seq, first);
				Last(refTupleList, (char*) sequences[s].seq, last);
				for (tupIndex= first; tupIndex < last; tupIndex++) {
					isMatch = 1;
					if (sequences[s].length + refTupleList[tupIndex].pos < refSeq.length) {
						if (refTupleList[tupIndex].strand == 0)
							sbjctSeq = (char*) refSeq.seq;
						else
							sbjctSeq = (char*) refSeqRC.seq;

						for (p = tupleSize; p < sequences[s].length and isMatch; p++ ){

							if (sequences[s].seq[p] != sbjctSeq[refTupleList[tupIndex].pos+ p])
								isMatch = 0;
						}
						
						if (isMatch) {
							if (refTupleList[tupIndex].strand == 0) {
								std::cout << refTupleList[tupIndex].pos << " " << sequences[s].length << " " 
													<< refTupleList[tupIndex].strand << std::endl;
							}
							else {
								std::cout << refSeqRC.length - refTupleList[tupIndex].pos - tupleSize + 1 
													<< " " << sequences[s].length << " " 
													<< refTupleList[tupIndex].strand << std::endl;
							}
						}
					}
				} // End ref seq list
			}			
		}
	}
}


