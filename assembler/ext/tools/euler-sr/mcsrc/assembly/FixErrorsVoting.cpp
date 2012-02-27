/***************************************************************************
 * Title:          FixErrorsVoting.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "compatibility.h"
#include "SimpleSequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "utils.h"
#include "hash/HashUtils.h"
#include "Tuple.h"
#include "ListSpectrum.h"
#include "StringMultTuple.h"
#include "StringTuple.h"
#include "BitSpectrum.h"
#include "BufferedSeqReader.h"
#include "IntegralTupleStatic.h"



#include <vector>
#include <iostream>
#include <ext/hash_map>
#include <map>
#include <deque>

char nextNuc[256];

void PrintUsage() {
	std::cout << "fixErrorsVoting   Fix errors in reads using spectral alignment with " << std::endl
						<< "                  a voting approach instead of a dynamic programming. " << std::endl
						<< "                  This has the benefit of being able to fix point " << std::endl
						<< "                  mutations or single indels in very small reads (<25 nt),"<< std::endl
						<< "                  assuming there is no more than one or two errors" << std::endl
						<< "                  per read." << std::endl;
	std::cout << "   Usage: fixErrorsVoting seqFile spectrumFile tupleSize outputFile [options] " << std::endl
						<< "     -minMult  m  Only consider tuples with multiplicity above m"<<std::endl
						<< "                    to be solid." << std::endl
						<< "     -minVotes v  Require at least 'v' separate votes to fix any position."<<std::endl
						<< "                  A vote is cast for a position p, nucleotide n, if a"<<std::endl
						<< "                  change at (p,n) makes a tuple t change from below m" <<std::endl
						<< "                  to above." << std::endl;
	std::cout << "     -maxTrim  x  Trim at most x nucleotides off the ends of a read. If "<<std::endl
						<< "                  more than x nucleotides need to be trimmed, the read is unfixable."
						<< std::endl;
	std::cout << "     -deletions   Search for single deletions. " << std::endl;
	std::cout << "     -insertions  Search for single insertions. " << std::endl;
	std::cout << "     -search s    Try up to 's' changes to make a tuple solid. " << std::endl;
	std::cout << "     -compare file For benchmarking purposes, the correct reads are given "<< std::endl
						<< "                  in file 'file'.  Read those and compare the results."<<std::endl;
	std::cout << "     -discardFile file Print all reads that do not pass the threshold to 'file'" 
						<< std::endl;
	std::cout << "     -map mapFile Print portions of retained reads to 'mapFile'."<< std::endl;
	std::cout << "     -spectrum [concise|full].  Use either a concise or full spectrum." <<std::endl
						<< "                    The concise spectrum must be on words of size less than 16" << std::endl
						<< "                    and resets all multiplicities greater than 3 to 3."<<std::endl
						<< "                    The full spectrum may be of any length, and stores " << std::endl
						<< "                    the exact multiplicity of all k-mers." << std::endl;
}

class SolidStats {
public:
	ssize_t nIns;
	ssize_t nDel;
	ssize_t nMut;
	ssize_t nTies;
	ssize_t nNotSolid;
	SolidStats() {
		nIns = nDel = nMut = 0;
		nTies = 0;
		nNotSolid = 0;
	}
};

void PrintMap(DNASequence &seq, ssize_t start, ssize_t end, std::ofstream &out);

template<typename T_Spectrum>
void TrimSequence(DNASequence &seq, T_Spectrum &spectrum,
									_INT_ maxTrim, ssize_t &seqStar, ssize_t &seqEnd);

template<typename T_Spectrum> 
ssize_t SolidSubsequence(DNASequence &seq, T_Spectrum &spectrum,
										 int tupleSize, ssize_t &seqStart, ssize_t &seqEnd);

template <typename T_Spectrum>
ssize_t FixSequence(DNASequence &seq, T_Spectrum &spectrum,
								IntMatrix &votes, IntVector &alreadySolid,
								int tupleSize, ssize_t voteThreshold, SolidStats &stats, ssize_t searchSize, ssize_t &changeMade);

ssize_t PrepareSequence(DNASequence &read);
void InitVotingMatrix(DNASequence &read, IntMatrix &votes);
void InitSolidVector(DNASequence &read, IntVector &solid);

template <typename T_Spectrum>
void VoteSequence(DNASequence &seq, T_Spectrum &spectrum, int tupleSize, ssize_t startPos,
									IntMatrix &votes, IntVector &solid, 
									ssize_t numSearch,
									ssize_t checkInsertions, ssize_t checkDeletions,
									std::deque<ssize_t> &history);

template <typename T_Spectrum>
ssize_t SolidifySequence(DNASequence &read, T_Spectrum &spectrum, int tupleSize,
										 IntMatrix &votes, IntVector &solid,
										 ssize_t minVotes, SolidStats &stats, ssize_t numSearch, ssize_t DoDeletion, ssize_t DoInsertion);

template <typename T_Spectrum>
ssize_t CheckSolid(DNASequence &seq, T_Spectrum &spectrum, int tupleSize);

int main(int argc, char* argv[]) {
 	nextNuc['G'] = 'A';
	nextNuc['A'] = 'C';
	nextNuc['C'] = 'T';
	nextNuc['T'] = 'G';

	int argi;
	if (argc < 4) {
		PrintUsage();
		exit(1);
	}
	argi = 1;
	std::string readsFile = argv[argi++];
	std::string spectrumFileName = argv[argi++];
	int tupleSize = atoi(argv[argi++]);
	std::string outputFileName  = argv[argi++];
	std::string discardFileName = "";
	std::string spectrumType = "full";
	std::string compareFile;
	std::string mapFileName;
	ssize_t printMap;
	ssize_t minMult, minVotes;
	ssize_t doInsertion, doDeletion;
	ssize_t doTrim;
	ssize_t maxTrim;
	ssize_t printAll;
	ssize_t maxMods;
	printMap    = 0;
	doTrim      = 0;
	doInsertion = 0;
	doDeletion  = 0;
	minMult     = 3;
	minVotes    = 4;
	maxTrim     = -1;
	printAll    = 0;
	maxMods     = 99999;
	ssize_t numSearch = 1;
	std::string reportFileName = FormReportName(readsFile);
	std::ofstream reportOut;
	openck(reportFileName, reportOut, std::ios::app, reportOut);
	BeginReport(argc, argv, reportOut);
	while (argi < argc) {
    if (strcmp(argv[argi], "-minMult") == 0 ) {
      ++argi;
      minMult = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-minVotes") == 0 ) {
      ++argi;
      minVotes = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-discardFile") == 0 ) {
      ++argi;
      discardFileName = argv[argi];
    }
		else if (strcmp(argv[argi], "-compare") == 0 ) {
			compareFile = argv[++argi];
		}
		else if (strcmp(argv[argi], "-search") == 0) {
			numSearch = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-deletions") == 0) {
			doDeletion = 1;
		}
		else if (strcmp(argv[argi], "-insertions") == 0) {
			doInsertion = 1;
		}
		else if (strcmp(argv[argi], "-trim") == 0) {
			doTrim = 1;
		}
		else if (strcmp(argv[argi], "-map") == 0) {
			printMap = 1;
			//			mapFileName = argv[++argi];
		}
		else if (strcmp(argv[argi], "-maxTrim") == 0) {
			doTrim = 1;
			maxTrim = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-spectrumType") == 0){ 
			spectrumType = argv[++argi];
		}
		else if (strcmp(argv[argi], "-printAll") == 0){ 
			printAll = 1;
		}
		else if (strcmp(argv[argi], "-maxMods") == 0) {
			maxMods = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
			exit(0);
		}
		argi++;
	}


	doDeletion = 0;
	doInsertion = 0;
  std::ifstream seqIn;
	//	BufferedSeqReader<1000> seqReader; 
	//	seqReader.Init(readsFile);
	//  openck(readsFile, seqIn, std::ios::in, reportOut);

  std::ofstream seqOut, discardOut, mapOut;
  openck(outputFileName, seqOut, std::ios::out, reportOut);
	openck(readsFile, seqIn, std::ios::in, reportOut);
	
  if (discardFileName != "") 
    openck(discardFileName, discardOut, std::ios::out, reportOut);

	/*	if (printMap) {
		openck(mapFileName, mapOut, std::ios::out, reportOut);
	}
	*/
	//  std::cout << "getting reads "; std::cout.flush();

	//	DNASequenceList reads, compare;
	//	ReadDNASequences(readsFile, reads, reportOut);
	
	//	if (compareFile != "") {
	//		ReadDNASequences(compareFile, compare, reportOut);
	//	}
	Spectrum<StringTuple>* spectrumPtr = (ListSpectrum<StringTuple>*) NULL;
	ListSpectrum<StringMultTuple> listSpectrum;

	// Configure the bit spectrum
	BitSpectrum<StringTuple> bitSpectrum(tupleSize);
	bitSpectrum.FindOnlySolid();

	// Configure the type of spectrum used.
	if (spectrumType == "full")
		spectrumPtr = (ListSpectrum<StringTuple>*) &listSpectrum; // TODO: Fix compiler "warning: dereferencing type-punned pointer will break strict-aliasing rules"
	else if (spectrumType == "concise") 
		spectrumPtr = &bitSpectrum;
	else {
		std::cout << "ERROR, the spectrum type must either be full or concise." << std::endl;
		std::cout << "you specified " << spectrumType << std::endl;
	}

	//  spectrumPtr->tupleSize = tupleSize;
	std::cout << "reading spectrum." << std::endl;
	spectrumPtr->Read(spectrumFileName, minMult);
	std::cout << "done." << std::endl;

	ssize_t r;
	IntMatrix votes;
	IntVector solid;
	SolidStats stats;
	ssize_t numDiscarded = 0;
	std::stringstream mapstrm;
	DNASequence read;
	//	for (r = 0; r < reads.size(); r++ ){
	ssize_t numFixed = 0;
	r = -1;
	while (SeqReader::GetSeq(seqIn, read, SeqReader::noConvert)) {
		if (r % 1000 == 0 and r > 0) {
			std::cout << ".";
			std::cout.flush();
		}
		if (r % 50000 == 0 and r > 0)
			std::cout << " " << r << " " << numFixed << " " << numDiscarded << std::endl;
		r++;
		ssize_t discardSeq = 0;
		if (!PrepareSequence(read)) {
			//			std::cout << "sequence : " << read.namestr << " is bad!" << std::endl;
			discardSeq = 1;
		}
		else {
			DNASequence original = read;
			ssize_t numChanges = 0;
			if ((numChanges = 
					 SolidifySequence(read, *spectrumPtr, tupleSize,
														votes, solid, minVotes, stats, numSearch, doDeletion, doInsertion)) != 0) {
				if (numChanges > maxMods) {
					std::cout << read.namestr << " " << numChanges << std::endl;
					discardSeq = 1;
					read = original;
				}
				else {
					if (printMap) {
						mapstrm.str("");
						mapstrm << " mapstart=0 mapend="<< read.length;
						read.namestr += mapstrm.str();
						//						PrintMap(read, 0, read.length, mapOut);
					}
					mapstrm.str("");
					mapstrm << " index=" << r;
					read.namestr += mapstrm.str();
					read.PrintlnSeq(seqOut);
				}
			}
			else {
				// Try trimming the sequence to get something that works.
				ssize_t trimStart, trimEnd;
				// Find the locations of the first solid positions.
				if (doTrim) {
					TrimSequence(read, *spectrumPtr, tupleSize,
											 trimStart, trimEnd);
					// If there is space for one solid tuple (trimStart < trimEnd - ts+1)
					// and the subsequence between the trimmed ends is ok, print the
					// trimmed coordinates.
					DNASequence fixedSeq;
					//UNUSED// int ss = SolidSubsequence(read, *spectrumPtr, tupleSize, trimStart, trimEnd);
					if ((maxTrim == -1 or 
							 (trimStart < maxTrim  and
								(read.length - trimEnd) < maxTrim)) and
							(trimStart < trimEnd - tupleSize + 1 and
							 SolidSubsequence(read, *spectrumPtr, tupleSize,
																trimStart, trimEnd))) {
						// Part of the sequence works, print the trimmed region
						if (printMap) {
							mapstrm.str("");
							mapstrm << " mapstart=" << trimStart << " mapend=" << trimEnd;
							read.namestr += mapstrm.str();
							//							PrintMap(read, trimStart, trimEnd, mapOut);
						}
						fixedSeq.CopyDetails(read);
						fixedSeq.Copy(read, trimStart, trimEnd);
						mapstrm.str("");
						mapstrm << " index=" << r;
						fixedSeq.namestr += mapstrm.str();
						fixedSeq.PrintlnSeq(seqOut);
						discardSeq = 0;
					}
					else {
						// Either too much was trimmed,
						// or the sequence was not solid between the trimming.
						// Consider this sequence unfixable.
						discardSeq = 1;  
					}
				}
				else {
					if (printAll) {
						read.PrintSeq(seqOut);
						seqOut << std::endl;
					}
					else {
						discardSeq = 1;
					}
				}
			}
		}
		
		if (discardSeq) {
			++numDiscarded;
			if (discardFileName !="") {
				read.PrintSeq(discardOut);
				discardOut << std::endl;
			}
		}
		else {
			++numFixed;
		}
	}

	std::cout << std::endl << "stats: " << stats.nMut << " " 
						<< stats.nDel << " " <<stats.nIns << std::endl;
	std::cout << stats.nTies << " " << stats.nNotSolid << std::endl;
	std::cout << "discarded " << numDiscarded << std::endl;
	EndReport(reportOut);
	reportOut.close();
	return 0;
}


template <typename T_Spectrum>
ssize_t SolidifySequence(DNASequence &read, T_Spectrum &spectrum, int tupleSize,
										 IntMatrix &votes, IntVector &solid,
										 ssize_t minVotes, SolidStats &stats, ssize_t numSearch, ssize_t DoDeletion, ssize_t DoInsertion) {
	
	ssize_t s;
	std::deque<ssize_t> history;
	//UNUSED// ssize_t changeMade;
	ssize_t startPos, fixPos;
	fixPos = -1;
	ssize_t iter = 0;
	ssize_t numFixed = 0;
	do {
		//		std::cout << "iter: " << iter << std::endl;

		if (fixPos > 0)
			startPos = fixPos;
		else 
			startPos = 0;

		for (s = 1; s <= numSearch; s++) {
			InitVotingMatrix(read, votes);
			InitSolidVector(read, solid);
			VoteSequence(read, spectrum, tupleSize, startPos, 
									 votes, solid, s, DoDeletion, DoInsertion, history);
			//			PrintMatrix(votes, std::cout, 2);
			++numFixed;
			if (FixSequence(read, spectrum, votes, solid, tupleSize, minVotes, stats, s, fixPos)) 
				return numFixed;
		}
		//		std::cout << "fp: " << fixPos << std::endl;
		++iter;
	} while (fixPos > 0);
	return 0;
}
 
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

template <typename T_Spectrum>
ssize_t CheckSolid(DNASequence &seq, T_Spectrum &spectrum, int tupleSize) {
	ssize_t p;
	typename T_Spectrum::TupleType tuple;
	for (p = 0; p < seq.length - tupleSize +1; p++ ) {
		tuple.assign((char*) &seq.seq[p]);
		if (spectrum.FindTuple(tuple) == -1) {
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

template <typename T_Spectrum>
void VoteSequence(DNASequence &seq, T_Spectrum &spectrum, int tupleSize,
									ssize_t startPos,
									IntMatrix &votes, IntVector &solid, 
									ssize_t numSearch,
									ssize_t checkInsertions, ssize_t checkDeletions,
									std::deque<ssize_t> &history) {

	DNASequence dnaseq;
	dnaseq.seq    = seq.seq;
	dnaseq.length = seq.length;
	dnaseq._ascii = 1;
	ssize_t p;
	/*
		std::cout << "fixing seq "  << seq.namestr << " ";
		seq.PrintSeq(std::cout); std::cout << std::endl;
	*/
	typename T_Spectrum::TupleType tempTuple;
	for (p = startPos; p < seq.length - tupleSize + 1; p++ ) {
		//		std::cout << "pos: " << p << std::endl;
		tempTuple.assign((char*) &seq.seq[p]);
		if (tempTuple.Valid()) {
			if (spectrum.FindTuple(tempTuple) != -1) {

				solid[p] = 1;
				/*				std::cout << "solid: " << p << std::endl;*/
			}
			else {
				// Cast votes for mutations
				ssize_t vp;
				//UNUSED// ssize_t nucIndex;
				unsigned char mutNuc;
				unsigned char un;
				ssize_t mut;
				//UNUSED// ssize_t histPos, histMut;
				for (vp = 1; vp < tupleSize - 1; vp++) {
					un = toupper(tempTuple[vp]);
					mutNuc = nextNuc[seq.seq[p + vp]];
					char start = seq.seq[p + vp];
					MutateTuple((char*) seq.seq, p + vp, mutNuc);

					for (mut = 0; mut < 3; mut++ ) {
						tempTuple.assign((char*) &seq.seq[p]);
						if (spectrum.FindTuple(tempTuple) != -1) {
							VoteHistory(votes, history);
							votes[vp + p][unmasked_nuc_index[mutNuc]]++;
							/*	
							std::cout << vp + p << " " 
												<< votes[vp + p][unmasked_nuc_index[mutNuc]] << std::endl;
							*/
						}
						else {
							if (numSearch > 1) {
								history.push_back(vp+p);
								history.push_back(mutNuc);

								VoteSequence(seq, spectrum, tupleSize, p + vp + 1, votes, solid, numSearch-1,
														 checkInsertions, checkDeletions, history);
								history.pop_back();
								history.pop_back();
							}
						}
						mutNuc = nextNuc[mutNuc];
						MutateTuple((char*) seq.seq, p + vp, mutNuc);
					}
					// Put this guy back
					//					MutateTuple((char*) seq.seq, p + vp, mutNuc);
					assert(seq.seq[p + vp] == start);
				}
				/*
				if (checkDeletions) {
					// Cast votes for deletions.  Don't delete at the boundaries
					// since this will just move to the next tuple
					if (p < seq.length - tupleSize) {
						Tuple delTup;
						// fill the vacated spot with the character after this tuple
						unsigned char fillChar = seq.seq[p + tupleSize];
						for (vp = 2; vp < tupleSize; vp++) {
							delTup = tuple;
							DeleteTuple((char*) delTup.c_str(), delTup.size(), vp, fillChar);
							if (FindKmer(delTup, spectrum) != -1) {
								VoteHistory(votes, history);
								votes[vp + p][4]++;
							}
							else {
								if (numSearch > 1) {
									history.push_back(vp + p);
									history.push_back(4);
									VoteSequence(seq, spectrum, tupleSize, votes, solid, numSearch-1,
															 checkInsertions, checkDeletions, history);
									history.pop_back();
									history.pop_back();
								}
							}
						}
					}
				}
				*/
				/*
					Fix this later.
					Use the faster nuc lookup instead of nuc->number->nuc.
				if (checkInsertions) {
					for (vp = 1; vp < tupleSize - 1; vp++ ) {
						un = toupper(tuple.c_str()[vp]);
						Tuple insTup;
						ssize_t ins;
						nucIndex = StartNuc(un);
						for (ins = 0; ins < 3; ins++) {
							insTup = tuple;
							InsertTuple((char*) insTup.c_str(), insTup.size(), vp, nuc_char[nucIndex]);
							if (FindKmer(insTup, spectrum) != -1) {
								VoteHistory(votes, history);
								votes[vp + p][5 + nucIndex]++;
							}
							else {
								if (numSearch > 1) {
									history.push_back(vp + p);
									history.push_back(5 + nucIndex);
									VoteSequence(seq, spectrum, tupleSize, votes, solid, numSearch-1,
															 checkInsertions, checkDeletions, history);

									history.pop_back();
									history.pop_back();
								}
							}
							nucIndex = NextNuc(nucIndex);
						}
					}
				}
				*/
			}
		}
	}
}

template <typename T_Spectrum>
ssize_t FixSequence(DNASequence &seq, 
								T_Spectrum &spectrum,
								IntMatrix &votes, IntVector &alreadySolid,
								int tupleSize, ssize_t voteThreshold, SolidStats &stats, ssize_t numSearch,
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
			if (votes[p][m] > voteThreshold) 
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
			//			std::cout << "orig seq" << std::endl;
			//			seq.PrintSeq(std::cout);
			//			std::cout << std::endl;
			if (maxPos.size() > 1 and maxPos[0] != maxPos[1]) {
				//				std::cout << "FIXING A TIE, fixing pos: " << maxPos[0] << std::endl;
			}
			unsigned char prev, cur;
			for (s = 0; s < numSearch and s < maxPos.size(); s++) {
				mod = maxMod[s];
				pos = maxPos[s];
				fixPos = pos;
				if (mod < 4) {
					prev = seq.seq[pos];
					cur = nuc_char[mod];
					seq.seq[pos] = nuc_char[mod];
					stats.nMut++;
				}
				else if (mod == 4) {
					ssize_t i;
					for (i = pos; i < seq.length; i++ ) {
						seq.seq[i] = seq.seq[i+1];
					}
					seq.length--;
					stats.nDel++;
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
					stats.nIns++;
				}
			}
			ssize_t solidRes = CheckSolid(seq, spectrum, tupleSize);
			/*
			std::cout << seq.namestr<< " found a fix " << solidRes << " " << prev << " " << cur 
								<< " " << numSearch << " " << fixPos << std::endl;
			*/
			//			seq.PrintSeq(std::cout);
			//			std::cout << std::endl;
			return solidRes;
		} 
		else {
			//			std::cout << "fix has " << numTies << " ties / " << numSearch << std::endl;
			stats.nTies++;
			return 0;
		}
	}
	else {
		//		std::cout << "none above threshold. " << std::endl;
		stats.nNotSolid++;
		return 0;
	}
}

template <typename T_Spectrum>
void TrimSequence(DNASequence &seq, T_Spectrum &spectrum,
									int tupleSize, ssize_t &seqStart, ssize_t &seqEnd) {
	ssize_t s0, e0;

	ssize_t i;
	seqStart = 0;
	s0 = seqStart;
	typename T_Spectrum::TupleType tempTuple;
	for (i = 0; i < seq.length - tupleSize + 1; i++ ){ 
		tempTuple.assign((char*) &seq.seq[i]);
		if (spectrum.FindTuple(tempTuple) != -1) {
			break;
		}
		// Not solid yet, advance
		seqStart++;
	}
	seqEnd = seq.length;
	e0 = seqEnd;
	for (i = seqStart + 1; i < seq.length - tupleSize; i++ ) {
		tempTuple.assign((char*) &seq.seq[i]);
		if (spectrum.FindTuple(tempTuple) == -1) {
			break;
		}
	}
	if (i == seq.length - tupleSize) 
		// The sequence is not trimmed.
		seqEnd = seq.length;
	else 
		// The sequence is trimmed. Trim end is the index of the first
		// 'bad' nucleotide. Since seqStart is the index of the first
		// 'good' nucleotide, seqEnd - seqStart is the length of the
		// untrimmed seq.  In other words, it's half open 0 based
		// indexing.
		seqEnd = i + tupleSize-1;

}
	
template <typename T_Spectrum>
ssize_t SolidSubsequence(DNASequence &seq, T_Spectrum &spectrum,
										 int tupleSize, ssize_t &seqStart, ssize_t &seqEnd) {
	ssize_t i;
	ssize_t solidSubsequence = 1;
	typename T_Spectrum::TupleType tempTuple;
	for (i = seqStart; i < seqEnd - tupleSize + 1; i++) {
		tempTuple.assign((char*) &seq.seq[i]);
		if (spectrum.FindTuple(tempTuple) == -1) {
			solidSubsequence = 0;
			break;
		}
	}
	return solidSubsequence;
}
		

	
void PrintMap(DNASequence &seq, ssize_t start, ssize_t end, std::ofstream &out) {
	out << ">" << seq.namestr << std::endl;
	out << start << " " << end << std::endl;
}
