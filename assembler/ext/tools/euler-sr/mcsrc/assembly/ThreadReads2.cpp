/***************************************************************************
 * Title:          ThreadReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DeBruijnGraph.h"
#include "IntervalGraph.h"
#include "SeqReader.h"
#include "DNASequence.h"
#include "ReadMap.h"
#include "ParseTitle.h"
#include "ThreadUtils.h"
#include "align/alignutils.h"
#include <map>

#include "IntegralTupleStatic.h"

#ifdef USE_CSA_
#include "bbbwt/BBBWTQuery.h"
#endif

using namespace std;

void PrintUsage() {
	cout << "usage: threadReads graphBase reads graphOut" << endl;
	cout << " -genome genome   Use 'genome' to validate threads against corrected sequence."
			 << endl
			 << " -maxBranches b(20) Only attempt to thread reads when there are fewer than b" << endl
			 << "                    branches from a vertex within 'maxThreadLength'." << endl
			 << " -minThreadGap g(2) Only consider a thread valid if the gap between " << endl
			 << "                    first best and second best is greater than 'g'" << endl
			 << " -maxThreadLength l(100)  Only consider threads when the thread length" << endl
			 << "                    is less than l nucleotides." << endl
			 << " -minThreadEnd e(4) Only consider a thread valid if it extends at least" << endl
			 << "                    'e' nucleotides after a branch." << endl
			 << " -catastrophicFraction cf (0.5) Consider the end of a read catastrophic if" << endl
			 << "                    more than cf*length(suffix) errors are in it. " << endl
			 << "                    'e' nucleotides into an edge." << endl
			 << " -minCatastropheLength len (3) Only look for catastrophes in suffixes more than 'len'" << endl
			 << " -onlyMultiEdge     Only print reads that are threaded through multiple edges." << endl
			 << " -allowGaps         Allow gaps in alignments." << endl;
	
	cout << " -maxScore s        Only consider a thread valid if it has less than 's'" << endl
			 << "                    errors." << endl;
	cout << " -verbose           Print why thread-paths are rejected. " << endl;
}

int main(int argc, char* argv[]) {

  if (argc < 3) {
    PrintUsage();
    exit(1);
  }
	int argi = 1;
  string graphBase         = argv[argi++];
	string readsFile         = argv[argi++];
	string graphOut          = argv[argi++];
	string genomeFileName = "";

	ssize_t minThreadGap    = 2;
	ssize_t maxThreadLength = 4;
	ssize_t minThreadEnd    = 4;
	ssize_t minTailLength   = 4;
	ssize_t maxScore        = 4;
	ssize_t printMultiEdge  = 0;
	double catastrophicFraction = 0.5;
	ssize_t minCatastropheLength = 3;
	ssize_t verbose = 0;
	ssize_t allowGaps = 0;
	map<ssize_t,ssize_t> pathLengthHist;
#ifdef USE_CSA_
	string csaFileName;
#endif
	while (argi < argc) {
		if (strcmp(argv[argi], "-genome") == 0) {
			genomeFileName = argv[++argi];
		}
#ifdef USE_CSA_
		else if (strcmp(argv[argi], "-csa") == 0) {
			csaFileName = argv[++argi];
		}
#endif
		else if (strcmp(argv[argi], "-minThreadGap") == 0) {
			minThreadGap = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-maxThreadLength") == 0){ 
			maxThreadLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minThreadEnd") == 0) {
			minThreadEnd = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-maxScore") == 0) {
			maxScore = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-catastrophicFraction") == 0) {
			catastrophicFraction = atof(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minCatastropheLength") == 0) {
			minCatastropheLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minTailLength") == 0) {
			minTailLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-verbose") == 0) {
			verbose = 1;
		}
		else if (strcmp(argv[argi], "-onlyMultiEdge") == 0) {
			printMultiEdge = 1;
		}
		else if (strcmp(argv[argi], "-allowGaps") == 0) {
			allowGaps = 1;
		}
		else {
			PrintUsage();
			cout << "bad option " << argv[argi] << endl;
			exit(1);
		}
		argi++;
	}
	IntervalGraph graph;
	//graph.vertexSize = vertexSize;
	std::string reportFileName = FormReportName(graphBase);
	ofstream reportOut;
	openck(reportFileName, reportOut, std::ios::app, reportOut);
	BeginReport(argc, argv, reportOut);
	int vertexSize;
	ReadIntervalGraph(graphBase, graph, vertexSize, reportOut);

	cout << "done reading graph." << endl;

	//	ofstream threadsOut;
	//openck(threadedReadsFile, threadsOut, ios::out, reportOut);

#ifdef USE_CSA_ 
	BBBWT csa;
	if (csaFileName != "") {
		BW::Read(csaFileName, csa);
	}
#endif
	//UNUSED// ssize_t r;
	string readName, origName;

	DNASequence genome;
	SeqReader::MaskRepeats();
	if (genomeFileName != "" ){
		SeqReader::GetSeq(genomeFileName, genome, SeqReader::noConvert);
	}

	//UNUSED// ssize_t edgeIndex;
	ssize_t readIndex = 0;
	char *edgeSequence;
	ssize_t edgeSequencLength;
	//UNUSED// ssize_t intv;
	vector<ssize_t> bestScores;
	ssize_t numBestScores = 2;
	bestScores.resize(2);

	ssize_t rejectedGood, rejectedBad, acceptedGood, acceptedBad;
	ssize_t rejectGap, rejectLength, rejectEnd;
	rejectGap = rejectLength = rejectEnd = 0;
	rejectedGood = rejectedBad = acceptedGood = acceptedBad = 0;
	ssize_t numThreaded = 0;
	ssize_t numPrinted  = 0;
	ssize_t totalThreadLength = 0;
	ssize_t numPrintedPaths = 0;
	ssize_t numCatastrophe = 0;
	ssize_t numMaxScore    = 0;
	ssize_t numTooLong     = 0;
	ssize_t numTooSimilar  = 0;
	ssize_t numTooEarly    = 0;
	IntMatrix scoreMat;
	//UNUSED// ssize_t origReadIndex;
	//UNUSED// ssize_t fixedReadIndex;
	map<ssize_t, ssize_t> threadLengthHist;
	
	ifstream readsIn;
	openck(readsFile, readsIn, ios::in, reportOut);

	DNASequence readPrefix;
	//UNUSED// ssize_t curOrigReadIndex = 0;
	// fetch first original read.
	ssize_t readPrefixIndex = -1;

	InitScoreMat(scoreMat, -1, 2);

	std::vector<std::vector<char*> > vertexBranches;
	std::vector<std::vector<ssize_t> > vertexBranchLengths;
	vertexBranches.resize(graph.vertices.size());
	vertexBranchLengths.resize(graph.vertices.size());
	ssize_t v;
	ssize_t maxBranches = 24;
	for (v = 0; v < graph.vertices.size(); v++ ){ 
		if (graph.vertices[v].OutDegree() > 1) {
			cout << "Collecting branches for: " << v << endl;
			std::vector<char*> branchSequences;
			std::vector<ssize_t> branchSeqLengths;
			branchSequences.resize(maxBranches);
			branchSeqLengths.resize(maxBranches);
			ssize_t curSeqIndex = 0;
			if (CollectBranchSequences(graph, 120, v, maxBranches,
																 branchSequences, branchSeqLengths,
																 "", curSeqIndex, 0) == 0)
			cout << endl << endl << endl;
		}
	}
	exit(0);


	DNASequence read;
	while(SeqReader::GetSeq(readsIn, read, SeqReader::noConvert)) {

		ssize_t edgeIndex, intv;
		// Don't thread reads that don't exist in the graph.

		ssize_t pathIndex = readIndex * 2;

		//
		// Don't thread this read if it has no map to the graph.
		//
		if (graph.pathLengths[pathIndex] == 0)
			continue;

		// Find the first spot where this read maps to the graph.
		edgeIndex = graph.paths[pathIndex][0].edge;
		intv      = graph.paths[pathIndex][0].index;

		ssize_t pathSeqLength = graph.GetPathLength(pathIndex);

		// The full-length read is represented here.
		if (pathSeqLength == read.length)
			continue;

		ssize_t pathLength = graph.pathLengths[pathIndex];
		ssize_t lastEdge = graph.paths[pathIndex][pathLength-1].edge;
		ssize_t lastIntv = graph.paths[pathIndex][pathLength-1].index;


		//
		// Check if this read cannot thread past the end of the last edge
		// on the path, so don't bother threading it.
		if ((*graph.edges[lastEdge].intervals)[lastIntv].edgePos + (read.length - pathLength) < 
				graph.edges[lastEdge].length)
			continue;
		

		// quick access to the edge
		edgeSequence = (char *) graph.edges[edgeIndex].seq.seq;
		edgeSequencLength = graph.edges[edgeIndex].seq.length;
		
		// Now, to try and thread reads through the graph,
		// find out what reads are mapped to this edge.
		// These are stored in the interval list.

		// Try and thread this read through the graph.

		// Try threading the read.
		// We need to pass in the graph (we will be searching hte graph)
		//           an index of the edge the read is on
		//           The interval that says where a read maps.
		//           The actual sequence of the read.
		//				readIndex = (*graph.edges[edgeIndex].intervals)[intv].read;
				
		// Only thread reads in the forward direction to save half the time.
		ssize_t prefixStart, prefixEnd;

		prefixStart = -1;
		prefixEnd   = -1;
		//				if(readIndex % 2 == 0) {
		// Find the coordinates of the original read in the graph.
		//					readIndex /= 2;
		ssize_t goodThread = 1;
		ssize_t origReadIndex = -1;
					
		//UNUSED// ssize_t mp;

		const char* readTail;
		ssize_t readTailLength;

		
		ssize_t firstEdge = graph.paths[pathIndex][0].edge;
		ssize_t firstIntv = graph.paths[pathIndex][0].index;
		//UNUSED// ssize_t firstVertex = graph.edges[firstEdge].src;
		//UNUSED// int firstIntvLength = (*graph.edges[firstEdge].intervals)[firstIntv].length;
		//UNUSED// ssize_t lastVertex = graph.edges[lastEdge].src;
		//UNUSED// int lastIntvLength = (*graph.edges[lastEdge].intervals)[lastIntv].length;

		readTail = (char*) &read.seq[prefixStart];
		readTailLength   = read.length - prefixStart;
					
		char *fixedCatTail;
		ssize_t fixedCatTailLength;

					
		ssize_t minThreadScore = -1;
		string threadedSeq;
		ThreadPath minThreadPath;
		DNASequence perfectRead;
		// Thread the mapped read.
					
		// 
		// If there is a reference genome, do some comparisons.
		//

		if (genome.length > 0) {
			ssize_t origReadPos;
			ssize_t strand;

			if (ParseKeyword<ssize_t>(read.namestr, string("pos"), origReadPos)) {
				perfectRead.Reset(read.length);
				if (ParseKeyword<ssize_t>(read.namestr, string("strand"), strand)) {
					ssize_t p;

					if (strand == 0) {
						for (p = 0; p < read.length ; p++) {
							perfectRead.seq[p] = genome.seq[origReadPos + p];
						}
					}
					else {
						for (p = 0; p < read.length; p++) {
							perfectRead.seq[perfectRead.length - p - 1] = comp_ascii[genome.seq[origReadPos + p]];
						}
					}
				}
			}
		}


		// Thread the read.
					
		ssize_t bestScore;
		for (bestScore = 0; bestScore < numBestScores; bestScore++) 
			bestScores[bestScore] = -1;
					
		//UNUSED// double maxErrorRate = 0.15;
					
					
		// Concatenate the fixed prefix, and the unfixed suffix,
		// then thread through the graph.
						
		//Make sure we are concatenating something.
		assert(readPrefix.length + readTailLength > 0);
		perfectRead.namestr = readPrefix.namestr;


		minThreadScore = Thread(graph, //lastEdge, lastIntv,
														firstEdge, firstIntv,
														readTail, readTailLength,
														minThreadPath, maxScore, maxThreadLength, perfectRead, 
														bestScores, numBestScores, scoreMat, origReadIndex, allowGaps);

		fixedCatTail = NULL;
		fixedCatTailLength = 0;

					
		if (minThreadScore != -1) {
			ThreadToSeq(graph, minThreadPath, threadedSeq);
		}
					
		// 
		// Run some heuristics to filter suspicious threads.
		//

		// Don't trust threads that are too long.
		ssize_t threadEnd = -1;
		ssize_t gap = -1;
		ssize_t threadPathSize = -1;
		ssize_t maxScoreExceeded = 0;
		ssize_t endsTooEarly     = 0;
		ssize_t tooSimilar       = 0;
		ssize_t noMinThreadFound = 0;
		ssize_t pathTooLong      = 0;
		//UNUSED// ssize_t isCatastrophe    = 0;
					
		if (minThreadScore != -1) {

			// 
			// Don't trust paths that branch too much.
			//
			threadPathSize = minThreadPath.size();
			if ((ssize_t) minThreadPath.size() - (ssize_t) graph.pathLengths[readPrefixIndex] >= maxThreadLength) {
				cout << "path too long " << (ssize_t) minThreadPath.size() - (ssize_t) graph.pathLengths[readPrefixIndex] << " " << maxThreadLength  << endl;
				goodThread  = 0;
				pathTooLong = minThreadPath.size();
			}
					 
			//
			// Don't trust threads that end too early.
			//
			if (minThreadPath.size() > 1) {
				//UNUSED// ssize_t lastPathIntv = minThreadPath.size() - 1;	
				ThreadPath::iterator pathIt = minThreadPath.end();
				pathIt--;
				ssize_t vertex;
				int vertexSize;
				ssize_t lastEdge;

				lastEdge = pathIt->edge;
				vertex = graph.edges[lastEdge].src;
				vertexSize = graph.vertices[vertex].vertexSize;
				threadEnd = pathIt->length - vertexSize;
				if (pathIt->length - vertexSize < minThreadEnd) {
					//								cout << "Rejecting thread ends too early. " << pathIt->length - vertexSize <<  endl;
					//								reads[readIndex].PrintlnSeq(cout);
					cout << "early end: " << pathIt->length - vertexSize << endl;
					goodThread = 0;
					endsTooEarly = 1;
					numTooEarly++;
				}
			}
						
			//
			// Don't trust threads that have too close to a
			// score to the next best path.
			//
			assert(bestScores[0] != -1); // by entering heree minThreadScore != -1, nor this
			gap = bestScores[1] - bestScores[0];
			if (numBestScores > 1 and
					bestScores[0] != -1 and
					bestScores[1] != -1 and
					bestScores[1] - bestScores[0] < minThreadGap) {
				cout << "Rejecting thread too similar. " 
						 << bestScores[0] << " " << bestScores[1] << endl;
															
							
				goodThread = 0;
				tooSimilar = bestScores[1] - bestScores[0];
				numTooSimilar++;
			}

			//
			// Find and truncate catastrophic reads.
			//
			//UNUSED// ssize_t threadLength = read.length - prefixEnd;

			//UNUSED// ssize_t threadPos;
			//UNUSED// ssize_t nErrors = 0;
			//UNUSED// ssize_t catastrophe = 0;
			//UNUSED// ssize_t suffixLength = 0;
			//UNUSED// ssize_t lenDiff = read.length - threadedSeq.size();
			//
			// Don't trust threads that have too many errors at the end.
			//
						
			if (minThreadScore >= maxScore and minThreadPath.size() > 1) {
				cout << "too many errors. " << minThreadScore <<  " " << readPrefix.namestr << endl;
				goodThread = 0;
				maxScoreExceeded = minThreadScore;
				numMaxScore++;
			}
		}
		else {
			cout << "no path found." << endl;
			goodThread = 0;
			if (origReadIndex >= 0) {
				//							assert(origReadIndex < originalReads.size());
				if (prefixEnd != read.length)
					noMinThreadFound = 1;
			}
		}

		// Print all sorts of diagnostics 
		if (genome.length > 0) {
			// 
			// Print a comparison with the genome if one is available.
			//CCL: Quited.

			string originalReadStr((const char*) read.seq, 
														 read.length);
			string readPrefixStr((const char*) readPrefix.seq,
													 readPrefix.length);
			string threadedSeqCopy(threadedSeq);
			if (goodThread == 0) {
			}
#ifdef USE_CSA_
			ssize_t low, high;
			DNASequence threadedSeqDNA;
			threadedSeqDNA.seq = (unsigned char*) threadedSeq.c_str();
			threadedSeqDNA.length = threadedSeq.size();
			BW::Query(threadedSeqDNA, csa, low, high);
			if (high - low > 0) {
				//							cout << " with exact match" << endl;
				if (goodThread) 
					acceptedGood++;
				else {
					rejectedGood++;
				}
			}
			else {
				//							cout << "no exact match." << endl;
				if (goodThread) {
					acceptedBad++;
					cout << "s:" << minThreadScore << " p: " 
							 << minThreadPath.size() << " e: " << threadEnd 
							 << " g: " << gap << " b: " << bestScores[1] << endl;
								

					ssize_t p;
					ssize_t nomm = 0;
					for (p = 0; p < origRead.length; p++ ) {
						if (originalReadStr[p] == perfectRead.seq[p]) 
							originalReadStr[p] = tolower(originalReadStr[p]);
						else 
							nomm++;
					}
					string perfectReadStr((const char*) perfectRead.seq, 
																perfectRead.length);
					cout << "per: " << perfectReadStr << endl;
					cout << "ori: " << originalReadStr << " " << nomm << endl;
					cout << "tru: ";
					for (p = 0; p < readStart; p++)
						cout << " ";
					for (p = 0; p < readPrefix.length; p++ )
						cout << readPrefix.seq[p];
					cout << endl;
						
					ssize_t t;
					for (t = readStart; t < readStart + threadedSeqCopy.size(); t++) 
						if (threadedSeqCopy[t-readStart] == perfectRead.seq[t])
							threadedSeqCopy[t-readStart] =tolower(threadedSeqCopy[t-readStart]);
					cout << "thr: ";
					for (t = 0; t < readStart ; t++) 
						cout << " ";
					cout << threadedSeqCopy << endl;

				}
				else {
					rejectedBad++;
				}
			}
#endif
		}						
		// End diagnostics.
		if (!goodThread) {
		}
		if (goodThread) {
			if (threadLengthHist.find(minThreadPath.size()) == threadLengthHist.end())
				threadLengthHist[minThreadPath.size()] = 1;
			else
				threadLengthHist[minThreadPath.size()]++;

			// don't bother printing this if it is on one edge
			// and we are only printing threads with multiple edges.
			if (printMultiEdge and minThreadPath.size() == 1)
				continue;

			totalThreadLength += threadedSeq.size();
			DNASequence threadedRead;
			ThreadPath prefixPath;
			graph.PathToThread(readPrefixIndex*2, prefixPath);
			string prefixPathSeq;
			prefixPath.pop_back();
			threadedRead.seq = (unsigned char*) threadedSeq.c_str();
			threadedRead.length = threadedSeq.size();
			threadedRead.namestr = readPrefix.namestr;
			//		threadedRead.PrintlnSeq(threadsOut);
						
			numPrinted++;

			numThreaded++;
			if (threadedRead.length < readPrefix.length ) {
			}
		}
		else {
			numPrintedPaths++;
			numPrinted++;
		}
	}
	// end reading ref sequences
	cout << "threaded: " << numThreaded << " printed paths: " << numPrintedPaths << endl;
	cout << "A total of " << numPrinted << " reads were printed." << endl;
	cout << "A total of " << totalThreadLength << " bases were threaded." << endl;
	if (genome.length > 0) {
		cout << "acceptedGood acceptedBad rejectedGood rejectedBad" << endl;
		cout << acceptedGood << " " << acceptedBad << " " 
				 << rejectedGood << " " << rejectedBad << endl;
	}
	cout << "nCatas: " << numCatastrophe 
			 << " nTooLong: " << numTooLong 
			 << " nMaxScore: " << numMaxScore
			 << " nTooSimilar: " << numTooSimilar 
			 << " nTooEarly: " << numTooEarly << endl;
	map<ssize_t,ssize_t>::iterator histit;
	for (histit = threadLengthHist.begin(); histit != threadLengthHist.end();
			 ++histit) {
		cout << (*histit).first << " " << (*histit).second << endl;
	}

	EndReport(reportOut);
	reportOut.close();

	return 0;
}





