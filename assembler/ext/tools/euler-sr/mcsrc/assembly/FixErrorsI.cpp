/***************************************************************************
 * Title:          FixErrorsI.cpp 
 * Author:         Mark Chaisson, Glenn Tesler
 * Created:        2007
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
using namespace std;
#ifdef _OPENMP
#include <omp.h>
#endif

#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "utils.h"
#include "hash/HashUtils.h"
#include "IntegralTupleStatic.h"
#include "SAP.h"
#include "Voting.h"
#include "RuleList.h"

#include <vector>
#include <iostream>
#include <ext/hash_map>
#include <map>

int SEdge::tupleSize;

void PrintUsage() {
	cout << "usage: fixErrors readsFile spectrumFile tupleSize outputFile [options] " << endl;
	cout << "options: " << endl;
	cout << " -minMult m    (20)    The minimum multiplicity k-mer to consider solid " << endl
			 << " -misMatch N   (1)     Cost to mutate a nucleotide " << endl
			 << " -minVotes v   (2)     Require at least 'v' separate votes to fix any position."<< endl
			 << "                       A vote is cast for a position p, nucleotide n, if a" << endl
			 << "                       change at (p,n) makes a tuple t change from below m" << endl
			 << "                       to above." << endl
			 << " -edgeLimit N  (3)     Only find fixes if they are more than N away from"
			 << "                       the edge." << endl << endl
			 << endl
			 << " -span N               Minimum span needed for solid region" << endl
			 << " -maxTrim   N  (7)     The maximum amount to clip off the ends of reads. " << endl
			 << endl
			 << "Indels: " << endl
			 << " -correctIndels        Try to correct indels" << endl
			 << " -noCorrectIndels      Don't try to correct indels" << endl
			 << " -titleRules file      Read rules from 'file' for determining which type" << endl
			 << "                       of error correction to apply based on read name format." << endl
			 << "                       Type=0 allows indels.  Type=1 doesn't. " << endl
			 << endl
			 << " -maxGap g     (7)     The maximum gap to consider (insertions or deletions). " << endl
			 << "                         Note: increasing maxGap may substantially increase" << endl
			 << "                         the run time. " << endl
			 << " -gapExtend N  (1)     Cost to extend a gap " << endl
			 << endl
			 << " -startScore score (3) Start fixing reads with 'score' max score.  This "<< endl
			 << "                       helps fixing reads that only have 1 or 2 errors." << endl << endl
			 << " -stepScore  step (2)  Fix reads starting at 'startScore', and increase the " << endl
			 << "                       search space by 'step' mutations until" << endl
			 << "                       'maxScore' is reached." << endl << endl
			 << " -maxScore  (infinity) Maximum score to permit extension of graph " << endl
			 << endl
			 << "Other: " << endl
			 << " -replaceN             Replace all unfixed 'N's with random nucleotides." << endl 
			 << " -numJobs N (1)        Start N jobs for fixing errors. " << endl
			 << " -discardFile 'file'   Name of file to output unfixable reads." << endl
			 << "                       Not specifying this keeps discards in original file." << endl
			 << " -readFixFile filename Print the numbers of fixes (indel, mutation) for" <<endl
			 << "                       each read to 'filename'." << endl
			 << endl;

	// TODO: Parameter discrepancies due to older versions of code:
	// "-lockFile" listed but not implemented;
	// -printMap, -maxMods: parsed but not implemented;
	// -trim: parsed.  If we defaulted doTrim=0 then -trim would be useful
	//                 (in setting doTrim=1), but doTrim=1 by default
	//                 so -trim is effectively a nop.
	

		//			 << " -lockFile file        Wait on lock 'file' to read the spectrum." << endl
		//			 << "                       This is helpful when many processes are reading the."<<endl 
		//			 << "                       same file at the same time off an NFS." << endl << endl

	//			 << " -readsPerJob M (0)    Give M jobs to each thread at at time, so that R " << endl
	//			 << "                       reads will be fixed in ceil(R/(M*N)) processes" << endl
	//			 << "                       The default value of 0 gives all reads to one job." << endl;
	
}

ssize_t GetReadType(DNASequence &read, RuleList &rules, ssize_t useNameRules, ssize_t defaultType) {
	ssize_t readType;
	// 
	// Look to see if we can deduce read type from the FASTA 
	// name.
	//
	ssize_t readRule;
	if (useNameRules) {
		if (!GetReadRule(rules, read.namestr, readRule)) {
			readType = defaultType;
		}
		else {
			if (readRule == -1)
				readType = defaultType;
			else
				readType = rules[readRule].type;
		}
	}
	else {
		// 
		// No way to deduce read type from the name, deduce from 
		// length
		if (read.length < 50) {
			readType = 1;
		}
		else {
			readType = 0;
		}
	}
	return readType;
}

ssize_t FixReadDP(DNASequence &read, 	CountedIntegralTupleDict &spectrum, FixParams &params, Stats &totalStats, IntMatrix &scoreMat) {
	Stats stats;
	//
	// Fix the read, but try to do so with a low score at first, since
	// that takes much less time.
	//
	//UNUSED// int fracMaxTrim = (int) (read.length * 0.2);
	ssize_t readWasModified = 0;			
	ssize_t readFixed = 0;
	ssize_t score;
	for (score = params.startScore; score <= params.maxScore and !readFixed; score+= params.stepScore ) {
		params.scoreThreshold = score;
		stats.Reset();
		readWasModified = 0;
		readFixed = SolidifyRead(read, spectrum, scoreMat, params, stats, readWasModified);
	} 
	totalStats += stats;
	return readFixed;
}

int main(int argc, char* argv[]) {
  if (argc < 4) {
    PrintUsage();
    exit(0);
  }
  int tupleSize;
  int argi = 1;
  string readsFile        = argv[argi++];
  string spectrumFileName = argv[argi++];
  tupleSize = atoi(argv[argi++]);
  string outputFileName   = argv[argi++];
  string discardFileName  = "";
  string readFixFile      = "";
	string titleRulesFile   = "";
  FixParams params;
	ssize_t maxMods;
	int fixIndels = -1; // -1: not specified; 0: allow indels; 1: don't
  params.gapOpen = 4;
  params.gapExtend = 1;
  params.maxGap = 7;
  params.span = 2;
  params.maxTrim = 7;
  params.edgeLimit = 3;
  params.scoreThreshold = 100000;
	//	IntegralTuple::tupleSize = tupleSize;
	IntegralTuple::SetTupleSize(tupleSize);
  ssize_t minMult;
  minMult = 20;
	maxMods     = 99999;
  params.misMatch = 1;
	ssize_t replaceN = 0;
  params.startScore = 3;
  params.stepScore  = 2;
	ssize_t useNameRules = 0;
	ssize_t minVotes = 2;
	ssize_t doTrim = 1;
	ssize_t maxTrim = 0;
	ssize_t printMap = 0;
	ssize_t numJobs = 1;
	//	ssize_t readsPerJob = 0;
  while (argi < argc) {
    if (strcmp(argv[argi], "-startScore") == 0) {
      params.startScore = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-stepScore") == 0) {
      params.stepScore = atoi(argv[++argi]);
    }
		else if (strcmp(argv[argi], "-printMap") == 0) {
			printMap = 1;
			//			mapFileName = argv[++argi];
		}
		else if (strcmp(argv[argi], "-maxMods") == 0) {
			maxMods = atoi(argv[++argi]);
		}
    else if (strcmp(argv[argi], "-maxTrim") == 0 ) {
      ++argi;
      params.maxTrim = atoi(argv[argi]);
			maxTrim = params.maxTrim;
			doTrim = 1;
    }
    else if (strcmp(argv[argi], "-minMult") == 0 ) {
      ++argi;
      minMult = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-span") == 0 ) {
      ++argi;
      params.span = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-maxGap") == 0) {
      ++argi;
      params.maxGap = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-gapExtend") == 0 ) {
      ++argi;
      params.gapExtend = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-misMatch") == 0 ) {
      ++argi;
      params.misMatch = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-maxScore") == 0 ) {
      ++argi;
      params.scoreThreshold = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-discardFile") == 0 ) {
      ++argi;
      discardFileName = argv[argi];
    }
    else if (strcmp(argv[argi], "-edgeLimit") == 0 ) {
      ++argi;
      params.edgeLimit = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-readFixFile") == 0 ) {
      ++argi;
      readFixFile = argv[argi];
    }
		else if (strcmp(argv[argi], "-titleRules") == 0 ) {
			titleRulesFile = argv[++argi];
		}
    else if (strcmp(argv[argi], "-minVotes") == 0 ) {
      ++argi;
      minVotes = atoi(argv[argi]);
    }
		else if (strcmp(argv[argi], "-trim") == 0) {
			doTrim = 1;
		}
		else if (strcmp(argv[argi], "-replaceN") == 0) {
			replaceN = 1;
		}
		else if (strcmp(argv[argi], "-numJobs") == 0) {
			numJobs = atoi(argv[++argi]);
		}
		//		else if (strcmp(argv[argi], "-readsPerJob") == 0) {
		//			readsPerJob = atoi(argv[++argi]);
		//		}
		else if (strcmp(argv[argi], "-correctIndels") == 0) {
			fixIndels = 0;
		}
		else if (strcmp(argv[argi], "-noCorrectIndels") == 0) {
			fixIndels = 1;
		}
    else {
      PrintUsage();
      cout << "Bad option: " << argv[argi] << endl;
			exit(1);
    }
    ++argi;
  }
	string reportFileName = FormReportName(readsFile);
	std::ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);

	// Determine number of jobs
#ifdef _OPENMP
	if (numJobs > 1) {
//		int maxJobs = omp_get_max_threads();
		int maxJobs = 1;
		if (numJobs > maxJobs) {
			cout << "Reducing from " << numJobs << " jobs to " << maxJobs << std::endl;
			numJobs = maxJobs;
		} else if (numJobs < maxJobs) {
//			omp_set_num_threads(numJobs);
		}
		cout << "Using " << numJobs << " threads" << std::endl;
	}
#else
	// not parallel mode, so only one job
	if (numJobs != 1) {
		cout << "'-numJobs " << numJobs << "' is only available if the code is compiled via 'gmake parallel'" << std::endl;
		numJobs = 1;
	}
#endif


	//

  ssize_t i, j;
  IntMatrix ScoreMat;
  CreateMatrix(ScoreMat, 5, 5);
  for (i = 0; i < 5; i++ ) {
    for (j = 0; j < 5; j++ ) {
      if (i != j) {
				ScoreMat[i][j] = params.misMatch;
      }
      else {
				ScoreMat[i][j] = 0;
      }
    }
  }
  for (i = 0; i < 4; i++) {
    ScoreMat[i][4] = params.misMatch;
    ScoreMat[4][i] = params.misMatch;
  }
  ofstream seqOut, discardOut, readFixOut;
  openck(outputFileName, seqOut, ios::out, report);

  if (discardFileName != "") 
    openck(discardFileName, discardOut, ios::out, report);
  if (readFixFile != "") 
    openck(readFixFile, readFixOut, ios::out, report);

	RuleList rules;

	if (titleRulesFile != "") {
		ParseRuleFile(titleRulesFile, rules, report);
		useNameRules = 1;
	}

	cout << "reading tuples..." <<endl;


	//UNUSED// ssize_t spectrumSize;

	

		/****
				 Old code that does not use the dictionary.
					vector<CountedIntegralTuple> spectrum;
					ReadMultBoundedBinaryTupleList(spectrumFileName, minMult, spectrum); 
		****/


	// Remove tupless less than minMult from the list.
	SEdge::tupleSize = tupleSize;
	//	IntegralTuple::tupleSize = tupleSize;
	IntegralTuple::SetTupleSize(tupleSize);

  //UNUSED// ssize_t s;
  Stats totalStats;
	std::vector<Stats> threadStats;


  // Now try to fix all reads.
  ssize_t r = 0;
  ssize_t numFixed = 0;
  vector<char> fixed;
	vector<char> modified;
	//UNUSED// ssize_t fixIteration = 0;

  params.maxScore = params.scoreThreshold;
  //UNUSED// ssize_t score;
  //UNUSED// ssize_t step;
  ssize_t nFixed; 
	ifstream readsIn;
	params.maxScore -= (params.maxScore - params.startScore) % params.stepScore;

	nFixed = 0; 
	openck(readsFile, readsIn, ios::in, report);
	r = 0;
	cout << "**************************************************";
	totalStats.PrintHeader(cout);

	//UNUSED// ssize_t noGaps;
	//UNUSED// ssize_t yesGaps;
	//UNUSED// ssize_t errorModel;
	ssize_t readType;
	ssize_t defaultType = 0;
	//UNUSED// ssize_t readRule;
	ssize_t job;
	std::vector<DNASequenceList> readBufferList;
	BufferedSeqReader<10000> bufferedSeqReader;
	bufferedSeqReader.seqIn = &readsIn;
	std::vector<std::vector<char> > readStatusList;
	readStatusList.resize(numJobs);
	readBufferList.resize(numJobs);
	threadStats.resize(numJobs);

	CountedIntegralTupleDict spectrum;
	spectrum.InitFromFile(spectrumFileName, minMult, report);


	while(readsIn and numJobs > 0) {
		
		//
		// Step 1. Load the sequence buffers from a file.
		//
		ssize_t bufferSize;

		//		cout << "charging buffers." << endl;
		//UNUSED// ssize_t blah;

		for (job = 0; job < numJobs; job++ ){ 
			bufferSize = bufferedSeqReader.Recharge(readBufferList[job]);
			if (bufferSize == 0) {
				numJobs = job;
				break;
			}
			readBufferList[job].resize(bufferSize);
			readStatusList[job].resize(bufferSize);
			std::fill(readStatusList[job].begin(), readStatusList[job].end(), 0);
		}
		//
		// Step 2. Fix all sequences in a sequence buffer.  
		//
		//		cout << "fixing reads." << endl;
#pragma omp parallel for
		for (job = 0; job < numJobs; job++ ) {

			DNASequence read;

			// ******** BEGIN PARALLEL CODE.  **********
			Stats stats;
			DNASequenceList seqList = readBufferList[job];
			DNASequenceList fixedList;
			std::vector<char> readStatus = readStatusList[job];
			ssize_t readFixed = 0;

			fixedList.resize(seqList.size());
			// Input:
			//  Not shared:
			//     seqList   An array of DNA Sequences.  Each sequence is a read and a title.
			//               Both the reads and the titles are of variable length.
			//               The number of reads per list is not always fixed, since the last
			//               call my have a truncated length if the number of reads is not divisible 
			//               by the read buffer length (the likely case).
			// 
			// 
			//  Shared:
			//     spectrum  A long (up to many GB) list of integer/count pairs as well as a dictionary
			//               index for faster lookup (this speeds up the method by >2X).
			//               
			//        rules  A list of regular expression rules used to distinguish what type
			//               of machine was used to sequence a read, since the algorithm used
			//               to fix errors changes depending on this.  Usually, it is safe to guess
			//               the read type based on length (short=Illumina, long=454 or Sanger, both of
			//               which use the same error correction).
			// 
			//       params  A structure containing all of the parameters for the method.
			// 
			//
			// Output:  (shared, of unknown size, and only updated by threads, they are added back
			//           to the master.)
			//
			//   fixedList   A list of reads of the same number of reads as seqList. The total size of the
			//               data structure may be different.
			//   
			//   readStatus  A list of length fixedList.size() that indicates whether or not the error
			//               correction method was successful on the data.
			//  threadStats  A structure containing a record of the total number of changes made
			//               to the list of reads.
			//
			//
			
			// 
			// The following variables are 

			ssize_t ri;
			IntMatrix votes;
			IntVector solid;

			for (ri = 0; ri < seqList.size(); ri++ ) {

				DNASequence read = seqList[ri];
				//SeqReader::GetSeq(readsIn, read, SeqReader::noConvert)) {


				if (fixIndels == -1) {
					readType = GetReadType(read, rules, useNameRules, defaultType);
				} else {
					readType = fixIndels;
				}
				
				//
				// Now fix the read based on which type it is.
				// 
				//  Type 0 allows indels.
				//  Type 1 does not.
				//
				readFixed = 0;
				if (readType == 0) {
					readFixed = FixReadDP(read, spectrum, params, threadStats[job], ScoreMat);
				}
				else if (readType == 1) {
					// 
					// Fix without indels, but do a more thorough search for a solid position.
					//
					ssize_t numChanges = 0;
					stats.Reset();
					//					cout << "solidifying " << read.namestr << endl;
					if ((numChanges = 
							 SolidifySequence(read, // private
																spectrum, // public
																tupleSize, // public  (not updated)
																minMult,   // public  (not updated)
																votes,     // private  
																solid, // private
																minVotes,  // public (not updated) 
																stats, // private
																1, 0, 0, 0 // obviously not changing.
																)) != 0) {
						readFixed = 1;
						//						cout << "stats: " << stats << endl;
						threadStats[job] += stats;
					}
				}
				if (discardFileName == "")
					readFixed = 1;
				
				if (replaceN) {
					for (i = 0; i < read.length; i++) {
						if (unmasked_nuc_index[read.seq[i]] >= 4)
							read.seq[i] = RandomNuc();
					}
				}
				readStatus[ri] = (char) readFixed;
				fixedList[ri]  = read;
			}

			// The fixed reads need to be copied back from the buffer in the 
			// child process back to the master so that they may be collected
			// and flushed to a file.

			readBufferList[job] = fixedList;
			readStatusList[job] = readStatus;
		}

		// 
		// ******* END PARALLEL ******
		//

		
		//
		// Step 3. collect the results, and flush them out to a file.
		//
		//		cout << "flushing reads. " << endl;
		ssize_t globReadFixed;
		for (job = 0; job < numJobs; job++) {
			ssize_t rji;
			totalStats += threadStats[job];

			for (rji = 0; rji < readBufferList[job].size(); rji++ ){ 
				r++;
				globReadFixed = readStatusList[job][rji];
				// 
				// Print some progress dots.
				//
				if (r % 1000 == 999) {
					cout << ".";
					cout.flush();
				}
				if (r % 50000 == 49999 )
					cout << " " << totalStats ;
				
				if (globReadFixed ) {
					readBufferList[job][rji].PrintlnSeq(seqOut);
					numFixed++;
					nFixed++;
				}
				else {
					if (discardFileName != "") {
						readBufferList[job][rji].PrintSeq(discardOut);
						discardOut << endl;
					}
				}
			}
		}
	}
	
	cout << endl;
	readsIn.close();
	readsIn.clear();
		
  cout << "------------------------------------------------------------" << endl;
  cout << "fixed: " << numFixed << endl;
  cout << "# mutations:  " << totalStats.numMut << endl;
  cout << "# insertions: " << totalStats.numIns << endl;
  cout << "# deletions:  " << totalStats.numDel << endl;
  cout << "Sequences were not fixed because of: " << endl;
  cout << "  " << totalStats.numNotSolid << " sequences did not have any solid tuples." << endl;
  cout << "  " << totalStats.numMultiplePaths << " sequences had multiple paths." << endl;
  cout << "  " << totalStats.numNoPathFound << " sequences had no valid support paths "<< endl;
	cout << "  " << totalStats.numErrorAtEnd << " could not be fixed due to an error within " << params.edgeLimit << " of an end "
						<< endl;

	//	IntegralTuple::DumpProcCounts(); // DBG

	EndReport(report);
	report.close();
}

