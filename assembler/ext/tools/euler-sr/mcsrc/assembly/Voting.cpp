/***************************************************************************
 * Title:          Voting.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifdef _OPENMP
#include <omp.h>
#endif


#include "compatibility.h"
#include "Voting.h"

void InitVotingMatrix(DNASequence &read, IntMatrix &votes) {
	if (votes.size() < read.length) {
		CreateMatrix(votes, read.length, 9);
	}
	else {
		//UNUSED+// ssize_t j;
		ssize_t i ;
		for (i = 0; i < votes.size(); i++ ){
			std::fill(votes[i].begin(), votes[i].end(), 0);
		}
	}
}

void InitSolidVector(DNASequence &read, IntVector &solid) {
	if (read.length > solid.size()) {
		solid.resize(read.length);
	}
	std::fill(solid.begin(), solid.end(), 0);
}

ssize_t PrepareSequence(DNASequence &read) {
	ssize_t p;
	for (p = 0; p < read.length; p++ ){ 
		read.seq[p] = toupper(read.seq[p]);
		if (!(read.seq[p] == 'A' ||
					read.seq[p] == 'C' ||
					read.seq[p] == 'T' || 
					read.seq[p] == 'G'))
			return 0;
	}
	return 1;
}


ssize_t SolidifySequence(DNASequence &read, 
										 CountedIntegralTupleDict &spectrum, int tupleSize, ssize_t minMult,
										 IntMatrix &votes, IntVector &solid,
										 ssize_t minVotes, Stats &stats, ssize_t numSearch, 
										 ssize_t DoDeletion, ssize_t DoInsertion, ssize_t earlyEnd) {
	
	ssize_t s;
	std::deque<ssize_t> history;
	//UNUSED// ssize_t changeMade;
	ssize_t startPos, fixPos;
	fixPos = -1;
	ssize_t iter = 0;
	ssize_t numFixed = 0;
	do {
		//		std::cout << endl << endl << read.namestr << " iter: " << iter << std::endl;
		if (fixPos > 0)
			startPos = fixPos;
		else 
			startPos = 0;

		for (s = 1; s <= numSearch; s++) {
			InitVotingMatrix(read, votes);
			InitSolidVector(read, solid);
			VoteSequence(read, spectrum, tupleSize, minMult, startPos, 
									 votes, solid, s, 0, DoDeletion, DoInsertion, earlyEnd, history);
			//			PrintMatrix(votes, std::cout, 2);
			if (FixSequence(read, spectrum, votes, solid, tupleSize, 
											minMult, minVotes, stats, s, fixPos)) {
				++numFixed;
				return numFixed;
			}
		}
		//		std::cout << "fp: " << fixPos << std::endl;
		++iter;
	} while (fixPos > 0);
	return 0;
}

ssize_t CheckSolid(DNASequence &seq, CountedIntegralTupleDict &spectrum, int tupleSize, ssize_t minMult) {
	ssize_t p;
	CountedIntegralTuple tuple;
	ssize_t index;
	for (p = 0; p < seq.length - tupleSize +1; p++ ) {
		if (tuple.StringToTuple(&(seq.seq[p]))) {
			if (((index = spectrum.DictLookupBinaryTuple(tuple)) == -1) or
					minMult > spectrum.list[index].count ) {
				return 0;
			}
		}
		else {
			return 0;
		}
	}
	return 1;
}


void VoteHistory(IntMatrix &votes, std::deque<ssize_t> &history) {
	ssize_t histPos, histMut;
	std::deque<ssize_t>::iterator histIt;
	// 
	for (histIt = history.begin(); histIt != history.end(); histIt++) {
		histPos = *histIt;
		++histIt;
		histMut = *histIt;
		votes[histPos][unmasked_nuc_index[histMut]]++;
	}
}

ssize_t VoteSequence(DNASequence &seq, CountedIntegralTupleDict &spectrum, int tupleSize,
								 ssize_t minMult, ssize_t startPos,
								 IntMatrix &votes, IntVector &solid, 
								 ssize_t numSearch, ssize_t runFast,
								 ssize_t checkInsertions, ssize_t checkDeletions, ssize_t earlyEnd,
								 std::deque<ssize_t> &history) {
	DNASequence dnaseq;
	dnaseq.seq    = seq.seq;
	dnaseq.length = seq.length;
	dnaseq._ascii = 1;
	ssize_t p;
	CountedIntegralTuple tempTuple, mutatedTuple;
	ssize_t tupleIndex;
	ssize_t end;
	if (earlyEnd) {
		end=  min(seq.length, (_SSZT_) earlyEnd) - tupleSize + 1; // TODO: remove cast
	}
	else {
		end = seq.length - tupleSize + 1;
	}
	static char nextNuc[256];

 	nextNuc['G'] = 'A';
	nextNuc['A'] = 'C';
	nextNuc['C'] = 'T';
	nextNuc['T'] = 'G';													
	nextNuc['g'] = 'A';
	nextNuc['a'] = 'C';
	nextNuc['c'] = 'T';
	nextNuc['t'] = 'G';
	ssize_t numCast = 0;
	for (p = 0; p < seq.length; p++) {
		assert(seq.seq[p] >= (ssize_t) 'A' and
					 seq.seq[p] <= (ssize_t) 'Z');
	}
	
	for (p = startPos; p < end; p++ ) {
		
		if (tempTuple.StringToTuple(&seq.seq[p])) {
			std::string tupleString, mutatedString;
			//			tempTuple.ToString(tupleString);
			
			if ((tupleIndex = spectrum.DictLookupBinaryTuple(tempTuple)) != -1 and
					spectrum.list[tupleIndex].count >= minMult ) {
			solid[p] = 1;
			}
			else {
				// Cast votes for mutations
				ssize_t vp;
				//UNUSED// ssize_t nucIndex;
				unsigned char mutNuc;
				//UNUSED// unsigned char un;
				ssize_t mut;
				//UNUSED// ssize_t histPos, histMut;
				for (vp = 0; vp < tupleSize - 1; vp++) {
					mutNuc = nextNuc[seq.seq[p + vp]];
					//UNUSED// char start = seq.seq[p + vp];
					MutateTuple(tempTuple, (unsigned char) numeric_nuc_index[mutNuc], vp, mutatedTuple);

					for (mut = 0; mut < 3; mut++ ) {
						if ((tupleIndex = spectrum.DictLookupBinaryTuple(mutatedTuple) != -1)) {
							if (spectrum.list[tupleIndex].count >= minMult) {
								votes[vp + p][unmasked_nuc_index[mutNuc]]++;
								++numCast;
							}
						}
						mutNuc = nextNuc[mutNuc];
						MutateTuple(tempTuple, numeric_nuc_index[mutNuc], vp, mutatedTuple);
					}
				}
			}
		}
	}
	return numCast;
}

ssize_t FixSequence(DNASequence &seq, 
								//								T_Spectrum &spectrum,
								CountedIntegralTupleDict &spectrum,
								IntMatrix &votes, IntVector &alreadySolid,
								int tupleSize, ssize_t minMult, ssize_t voteThreshold, Stats &stats, ssize_t numSearch,
								ssize_t &fixPos) {
	// numSearch is the number of mutations to search for to fix a read.
	
	// At first, no changes are made to the sequence
  fixPos = 0;
	ssize_t p, m;
	ssize_t numAboveThreshold = 0;
	ssize_t maxVotes = 0;
	ssize_t allGood  = 1;
	for (p = 0; p < seq.length - tupleSize + 1; p++ ) {
		if (alreadySolid[p] == 0) {
			allGood = 0;
			break;
		}
	}
	if (allGood) {
		// no need to fix this sequence
		/*		std::cout << "seq: " << seq.namestr << " is all good" << std::endl;*/
		return 1;
	}

	ssize_t s;
	//	PrintMatrix(votes, std::cout, 3);
	for (p = 0; p < votes.size(); p++) { 
		for (m = 0; m < votes[p].size(); m++) {
			if (votes[p][m] >= voteThreshold) 
				numAboveThreshold++;
			if (votes[p][m] >= maxVotes) {
				// Make room for the next vote
				maxVotes = votes[p][m];
			}
		}
	}
	
	//	std::cout << "max votes: " << maxVotes << std::endl;
	// Make sure there aren't multiple possible fixes
	std::vector<ssize_t> maxPos, maxMod;
	std::vector<ssize_t> tiePos;
	ssize_t numTies = -1;
	for (p = 0; p < votes.size(); p++) { 
		for (m = 0; m < votes[p].size(); m++) {
			if (votes[p][m] == maxVotes) {
				numTies++;
				maxPos.push_back(p);
				maxMod.push_back(m);
			}
		}
	}
	ssize_t mod, pos;
	
	if (numAboveThreshold > 0 ) {
		if (numTies < numSearch or 
				(maxPos.size() > 1 and maxPos[0] != maxPos[1])) {
			// Found at least one change to the sequence
			unsigned char prev, cur;
			for (s = 0; s < numSearch and s < maxPos.size(); s++) {
				mod = maxMod[s];
				pos = maxPos[s];
				fixPos = pos;
				if (mod < 4) {
					prev = seq.seq[pos];
					cur = nuc_char[mod];
					//					cout << "making modification " << seq.seq[pos] << " " << pos << " " << nuc_char[mod] << endl;
					seq.seq[pos] = nuc_char[mod];
					stats.numMut++;
				}
				else if (mod == 4) {
					ssize_t i;
					for (i = pos; i < seq.length; i++ ) {
						seq.seq[i] = seq.seq[i+1];
					}
					seq.length--;
					stats.numDel++;
				}
				else if (mod > 5) {
					seq.length++;
					unsigned char* newSeq = new unsigned char[seq.length];
					assert(pos > 0);
					memcpy(newSeq, seq.seq, pos);
					ssize_t i;
					for (i = pos + 1; i < seq.length; i++) {
						newSeq[i] = seq.seq[i-1];
					}
					newSeq[pos] = nuc_char[mod-5];
					delete[] seq.seq;
					seq.seq = newSeq;
					stats.numIns++;
				}
			}
			ssize_t solidRes = CheckSolid(seq, spectrum, tupleSize, minMult);

			return solidRes;
		} 
		else {
			// Bookkeeping for sequences that have ties.
			stats.numTies++;
			return 0;
		}
	}
	else {

		stats.numNotSolid++;
		return 0;
	}
}

ssize_t TrimSequence(DNASequence &seq, CountedIntegralTupleDict &spectrum,	int tupleSize, 
									ssize_t minMult, ssize_t &seqStart, ssize_t &seqEnd, ssize_t maxTrim, ssize_t printMap) {
	ssize_t i;
	seqStart = 0;
	//	typename T_Spectrum::TupleType tempTuple;
	CountedIntegralTuple tuple;
	ssize_t lastInvalid = -1;
	for (i = 0; (i < seq.length - tupleSize + 1) and (i < tupleSize) ; i++ ){ 
		if (nucToIndex[seq.seq[i]] >= 4)
			lastInvalid = i;
	}

	for (i = 0; i < seq.length - tupleSize + 1; i++ ){ 
		
		if ((i > lastInvalid) and
				tuple.StringToTuple(&seq.seq[i]) and
				spectrum.DictLookupBinaryTuple(tuple) != -1) {
			break;
		}
		if (nucToIndex[seq.seq[i+tupleSize-1]] >= 4)
			lastInvalid = i + tupleSize - 1;
		// Not solid yet, advance
		seqStart++;
	}
	seqEnd = seq.length;
	
	for (i = seqStart + 1; i < seq.length - tupleSize; i++ ) {
		if (nucToIndex[seq.seq[i+tupleSize-1]] >= 4) {
			// don't include the last tuple.
			i--;
			break;
		}
		if (tuple.StringToTuple(&seq.seq[i])) {
			if (spectrum.DictLookupBinaryTuple(tuple) == -1) {
				break;
			}
		}
	}

	// 
	// i is set to the position after the last solid tuple.  Since seqEnd - seqStart 
	// should be the length of the solid sequence this includes the length
	// of the tuple.
	//
	seqEnd = i + tupleSize;

	if (seqStart > maxTrim)
		return 0;
	
	if (seq.length - seqEnd > maxTrim)
		return 0;

	if (printMap) {
		std::stringstream mapstrm;
		mapstrm << " mapstart= "<< seqStart << " mapend=" << seqEnd;
		seq.namestr += mapstrm.str();
	}
	
	ssize_t s;
	ssize_t newLength = seqEnd - seqStart;
	for (s = 0; s < newLength; s++ ) {
		seq.seq[s] = seq.seq[s + seqStart];
	}
	seq.length = newLength;
	return 1;
}
	
ssize_t FindSolidSubsequence(DNASequence &seq, CountedIntegralTupleDict &spectrum, 
												 int tupleSize, ssize_t minMult, ssize_t &seqStart, ssize_t &seqEnd) {
	ssize_t i;
	ssize_t solidSubsequence = 1;
	CountedIntegralTuple tuple;
	for (i = seqStart; i < seqEnd - tupleSize + 1; i++) {
		if (!tuple.StringToTuple(&seq.seq[i]) or
				spectrum.DictLookupBinaryTuple(tuple) == -1) {
			solidSubsequence = 0;
			break;
		}
	}
	return solidSubsequence;
}
