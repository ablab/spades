/***************************************************************************
 * Title:          ThreadUtils.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "ThreadUtils.h"
#include "align/alignutils.h"
#include "AlignmentPrinter.h"


using namespace std;

void AlignReadToBranches(DNASequence &read,
												 ssize_t startPos,
												 std::vector<char*> &branchSequences,
												 std::vector<ssize_t> &branchLengths,
												 IntMatrix &scores,
												 IntMatrix &paths,
												 IntMatrix &scoreMat,
												 ssize_t band) {
	DNASequence readTail, branchSeq;
	ssize_t b;
	vector<ssize_t> branchScores;
	branchScores.resize(branchSequences.size());
	readTail.seq = &read.seq[startPos];
	readTail.length = read.length - startPos;
	
	// 
	// Align the read tail to the possible branches.
	//
	for (b = 0; b < branchSequences.size(); b++ ){
		branchSeq.seq    = (unsigned char*) &branchSequences[b];
		branchSeq.length = branchLengths[b];
		ssize_t *blank = NULL;
		branchScores[b] = BandedAlign(readTail, branchSeq, -1, 1, 1, band, 
																	blank, scores, paths, scoreMat, NULL);
	}
	
	//
	// Look for a spread between opt path and other paths.
	//
	ssize_t bestScore = 99999999, bestIndex = -1;
	ssize_t secondBest= 99999999, secondIndex = -1;
	for (b = 0; b < branchSequences.size(); b++ ) {
		if (bestScore > branchScores[b]) {
			secondBest  = bestScore;
			secondIndex = bestIndex;
			bestScore   = branchScores[b];
		}
	}
}
												


ssize_t CollectBranchSequences(IntervalGraph &graph, ssize_t branchLength, ssize_t vertex, ssize_t maxBranches,
													 std::vector<char*> &branchSequences, 
													 std::vector<ssize_t> &branchLengths, 
													 string curSeq, ssize_t &curSeqIndex, ssize_t depth) {
	ssize_t outEdge, outEdgeIndex;
	ssize_t dest;
	for (outEdgeIndex = graph.vertices[vertex].FirstOut();
			 outEdgeIndex != graph.vertices[vertex].EndOut();
			 outEdgeIndex = graph.vertices[vertex].NextOut(outEdgeIndex)) {
		if (curSeqIndex >= maxBranches)
			return 0;
		outEdge = graph.vertices[vertex].out[outEdgeIndex];
		dest = graph.edges[outEdge].dest;
		// add the sequence of the out edge.
		string branchSeq = curSeq;
		if (branchLength - curSeq.size() < graph.edges[outEdge].length) {
			branchSeq.append((const char*) graph.edges[outEdge].seq.seq, branchLength - curSeq.size());
			branchSequences[curSeqIndex] = new char[branchSeq.size()];
			branchSeq.copy(branchSequences[curSeqIndex], branchSeq.size());
			branchLengths[curSeqIndex] = branchSeq.size();
			++curSeqIndex;
		}
		else {
			branchSeq.append((const char*) graph.edges[outEdge].seq.seq, 
											 graph.edges[outEdge].length - 
											 graph.vertices[dest].vertexSize);
			if (CollectBranchSequences(graph, branchLength, dest, maxBranches,
																 branchSequences, branchLengths, branchSeq, curSeqIndex, depth + 1) == 0) {
				return 0;
			}
		}
	}
	// Not too many branches were found.
	return 1;
}

ssize_t Thread(IntervalGraph &graph,ssize_t edgeIndex, ssize_t intvIndex,
					 const char* read, ssize_t readLength, ThreadPath &minThreadPath, 
					 ssize_t maxScore, ssize_t maxDepth,
					 DNASequence &perfectRead, 
					 std::vector<ssize_t> &scoreList, ssize_t scoreListLength,
					 IntMatrix &scoreMat, ssize_t readIndex, ssize_t allowGaps) {
	std::cout << "Beginning thread for " << edgeIndex << " " << intvIndex << std::endl;
	std::vector<std::string> pathSequences;
	ThreadPath curThreadPath;
	ssize_t curThreadScore, minThreadScore;
	minThreadScore = -1;
	curThreadScore = 0;
	ssize_t foundThread = 0;

	DNASequence readSeq;
	readSeq.seq = (unsigned char*) read;
	readSeq.length = readLength;
	char *threadSeq = new char[readSeq.length];
	map<ssize_t,ssize_t> nTraverse;
	ThreadRecursivelyBandedAlign(readSeq, graph, edgeIndex, 
															 (*graph.edges[edgeIndex].intervals)[intvIndex].edgePos,
															 threadSeq, 0,
															 curThreadPath, curThreadScore,
															 minThreadPath, minThreadScore, 0, maxDepth, nTraverse, "");
															 
	if (minThreadPath.size() > 1) {
		ThreadPath::iterator thIt;
		std::cout << readIndex << "  min thread path: " << std::endl;
		for (thIt = minThreadPath.begin(); thIt != minThreadPath.end(); ++thIt) {
			std::cout << (*thIt).edge << " ";
		}
		std::cout << std::endl;
	}
	if (foundThread == 0) {
		//		std::cout << "no thread found. " << std::endl;
	}
	return minThreadScore;
}

void ThreadToSeq(IntervalGraph &graph,
								 ThreadPath &path,
								 std::string &sequence) {
	sequence = "";
	ThreadPath::iterator threadIt, threadEnd;
	ssize_t edge, edgeEnd;
	//		std::cout << "tl: ";
	for (threadIt = path.begin(); threadIt != path.end(); ++threadIt) {
		edge = (*threadIt).edge;
		edgeEnd = (*threadIt).pos + (*threadIt).length;
		if (edgeEnd > graph.edges[edge].length)
			edgeEnd = graph.edges[edge].length;

		ssize_t edgePos;
		//				std::cout << (*threadIt).length << " ";
		for (edgePos = (*threadIt).pos; edgePos < edgeEnd; ++edgePos) {
			sequence += graph.edges[edge].seq.seq[edgePos];
		}
	}
	//	std::cout << std::endl;
}


ssize_t ThreadRecursivelyBandedAlign(DNASequence &read,
																 IntervalGraph &graph,
																 ssize_t edgeIndex, ssize_t edgePos,
																 char* curPathSeq, ssize_t curPathSeqLength,
																 ThreadPath &curThreadPath, ssize_t curThreadScore,
																 ThreadPath &minThreadPath, ssize_t minThreadScore,
																 ssize_t curDepth, ssize_t maxDepth, map<ssize_t,ssize_t> &nTraverse, string spacing) {
	
	static ssize_t *locations = NULL;
	static ssize_t locationsLength = 0;
	DNASequence pathSeq;
	ssize_t dest = graph.edges[edgeIndex].dest;
	if (nTraverse.find(dest) == nTraverse.end()) {
		nTraverse[dest] = 1;
	}
	else {
		nTraverse[dest]++;
		if (nTraverse[dest] > 2) 
			return 0;
	}

	
	if (locations == NULL) {
		locations = new ssize_t[read.length + 1];
		locationsLength = read.length + 1;
	}
	else if (locationsLength < read.length + 1) {
		delete[] locations;
		locations = new ssize_t[read.length + 1];
		locationsLength = read.length + 1;
	}
	
	if (curPathSeqLength + graph.edges[edgeIndex].length - edgePos >= read.length) {

		// The sequence fits on this edge, align it.
		ssize_t edgePathLength = read.length - curPathSeqLength;
		assert(edgePathLength <= graph.edges[edgeIndex].length - edgePos);
		memcpy((char*) &curPathSeq[curPathSeqLength], 
					 (char*) &graph.edges[edgeIndex].seq.seq[edgePos], edgePathLength);
		
		
		pathSeq.seq = (unsigned char*) curPathSeq;
		pathSeq.length = read.length;
		ssize_t score = (ssize_t) BandedAlign(read, pathSeq, -2, 1, 3, 10, locations);

		ssize_t nMatches = CalcNumMatches(read, pathSeq, locations);
		
		assert(pathSeq.length > 0);
		if (pathSeq.length > 15 and 
				(double(nMatches) / double(pathSeq.length) < 0.70)) {
			return 0;
		}
		if (curThreadPath.size() > 1) {

		cout << spacing << "for " << curThreadPath.size() << " edges." << endl;
		ThreadPath::iterator pathIt;
		for (pathIt = curThreadPath.begin(); pathIt != curThreadPath.end(); ++pathIt) {
			cout << (*pathIt).edge << " " ;
		}
		cout << endl;

		cout << "path seq len: " << pathSeq.length << " nmatches: " << nMatches << endl;
		cout << spacing << "SCORE: " << score << endl;
		read.PrintlnSeq(std::cout);
		pathSeq.PrintlnSeq(std::cout);
		PrintAlignment(read, pathSeq, 0, 0, locations, read.length, std::cout);
		IntMatrix scoreMat;
		score = (ssize_t) AffineAlign(read, pathSeq, -2, 1, 3, 3, 1, locations, scoreMat);
		cout << spacing << "SCORE (Affine): " << score << endl;
		read.PrintlnSeq(std::cout);
		pathSeq.PrintlnSeq(std::cout);
		PrintAlignment(read, pathSeq, 0, 0, locations, read.length, std::cout);
		
		score = (ssize_t) Align(read, pathSeq, -2, 1, 1, locations, scoreMat);
		cout << spacing << "SCORE (needle): " << score << endl;
		read.PrintlnSeq(std::cout);
		pathSeq.PrintlnSeq(std::cout);
		PrintAlignment(read, pathSeq, 0, 0, locations, read.length, std::cout);
		}
		return 1;
	}
	else {
		if (curDepth >= maxDepth)
			return 0;

		// Add part of this edge to the curPathSeq
		ssize_t destVertex = graph.edges[edgeIndex].dest;
		ssize_t destLength = graph.vertices[destVertex].vertexSize;
		ssize_t pathSegmentLength = graph.edges[edgeIndex].length - edgePos - destLength;

		// Don't try fitting reads that are stuffed onto edges for now.
		if (pathSegmentLength <= 0)
			return 0;
		if (pathSegmentLength > 0) {
			assert(pathSegmentLength > 0);
			memcpy((char*) &(curPathSeq[curPathSeqLength]), 
						 (char*) &(graph.edges[edgeIndex].seq.seq[edgePos]),
						 pathSegmentLength);
			
			pathSeq.seq = (unsigned char*) curPathSeq;
			pathSeq.length = curPathSeqLength + pathSegmentLength;
			DNASequence readPrefix;
			readPrefix.seq = read.seq;
			readPrefix.length = pathSeq.length;
			
			ssize_t score = (ssize_t) BandedAlign(readPrefix, pathSeq, -2, 1, 3, 10, 
																		locations);
			ssize_t nMatches = CalcNumMatches(readPrefix, pathSeq, locations);
			
			if (pathSeq.length > 15 and 
					(double(nMatches) / double(pathSeq.length) < 0.70)) {
				cout << spacing << "this branch probably sucks!!!" << endl;
				return 0;
			}
		
			cout << spacing << curThreadPath.size()
					 << " score of prefix of length: " << pathSeq.length << " " << score << endl;
		}
		curPathSeqLength += pathSegmentLength;
		curThreadPath.push_back(ThreadPathInterval(edgeIndex, pathSegmentLength, edgePos));
		ssize_t outEdge, outEdgeIndex;
		std::cout << spacing << "searchign " << destVertex 
							<< " through " << graph.vertices[destVertex].OutDegree() << " out edges." << endl;
		for (outEdgeIndex = graph.vertices[destVertex].FirstOut();
				 outEdgeIndex != graph.vertices[destVertex].EndOut();
				 outEdgeIndex = graph.vertices[destVertex].NextOut(outEdgeIndex)) {
			outEdge = graph.vertices[destVertex].out[outEdgeIndex];

			ThreadRecursivelyBandedAlign(read, graph, outEdge, 0,
																	 curPathSeq, curPathSeqLength,
																	 curThreadPath, curThreadScore,
																	 minThreadPath, minThreadScore, curDepth + 1, maxDepth, 
																	 nTraverse, spacing + " ");
		}
		curThreadPath.pop_back();
	}
	return 1; // TODO: check
}


ssize_t ThreadRecursively(IntervalGraph &graph,
											ssize_t edgeIndex, ssize_t edgePos,
											const char*curSeq,ssize_t curSeqLength,
											ThreadPath &curThreadPath, ssize_t curThreadScore,
											ThreadPath &minThreadPath, ssize_t &minThreadScore, 
											ssize_t maxScore, ssize_t maxDepth,
											DNASequence &perfectRead,
											std::vector<ssize_t> &scoreList, ssize_t scoreListLength, 
											ssize_t depth, IntMatrix &scoreMat, ssize_t printNext, ssize_t readIndex, ssize_t allowGaps) {

	// Thread this sequence into the path;
	ssize_t seqPos;
	const char *edgeSeq = (const char *) graph.edges[edgeIndex].seq.seq;
	ssize_t edgeLength = graph.edges[edgeIndex].length;
	if (printNext) {
		std::cout << "at the next iteration:" << std::endl;
		std::string edgeStr;
		edgeStr.append(edgeSeq, std::min(edgeLength, curSeqLength));
		std::string curSeqStr;
		curSeqStr.append(curSeq, curSeqLength);
		std::cout << edgeStr << std::endl;
		std::cout << curSeqStr << std::endl;
		printNext = 0;
	}
	ssize_t destVertex;
	destVertex =graph.edges[edgeIndex].dest;
	// Advance forward in the sequence
	ssize_t destVertexLength = graph.vertices[destVertex].vertexSize;
	ssize_t threadLength;
	
	//
	// If this sequence may be threaded on the entire edge, do that and finish.
	// 
	if (edgePos + curSeqLength <= edgeLength) {
		threadLength = curSeqLength;
	}
	else {
		threadLength = edgeLength - destVertexLength - edgePos;
	}
	if (threadLength <= 0)
		return 0;

	DNASequence threadAlignSeq, edgeAlignSeq;

	// we want to thread the sequence all the way to the end of the edge
	ssize_t threadAlignLength;
	if (threadLength < curSeqLength)
		threadAlignLength = threadLength + destVertexLength;
	else 
		threadAlignLength = threadLength;
	
	double alignScore;
	ssize_t *alignLocations = NULL;
	// Try local alignment if the path doesn't fit on the whole edge
	// 
	if (allowGaps and (curThreadPath.size() > 0 or threadAlignLength != threadLength)) {
		threadAlignSeq.seq = (unsigned char *) &curSeq[0];
		threadAlignSeq.length = threadAlignLength;
		edgeAlignSeq.seq = (unsigned char*) &edgeSeq[edgePos];
		edgeAlignSeq.length = threadAlignLength;
	
		// use alignment that only has affine gaps.
		//		alignScore =  AffineAlign(threadAlignSeq, edgeAlignSeq, -1, 1, 1000, 15, 3, alignLocations, scoreMat);
		alignScore = LocalAlign(threadAlignSeq, edgeAlignSeq, -1, 1, 5, alignLocations, scoreMat);
		
	}

	ssize_t curHammingNumMismatch = 0;

	// Store the number of mismatches and gaps in the current alignment
	// that are for the same distance that is matched on the hamming
	// distance alignment.
	ssize_t curAlignNumMismatch = 0;
	ssize_t curAlignNumGap = 0;
	for (seqPos = 0; seqPos < threadLength; seqPos++) {
		//		std::cout << edgeSeq[seqPos + edgePos];
		if (numeric_nuc_index[(unsigned char) curSeq[seqPos]] != numeric_nuc_index[(unsigned char) edgeSeq[seqPos + edgePos]])
			curHammingNumMismatch++;
	}

	ssize_t readPosGapOffset = 0;
	ssize_t edgePosGapOffset = 0;
	ssize_t useLocalAlignment = 0; 
	ssize_t numReadGaps, numEdgeGaps;
	numReadGaps = numEdgeGaps = 0;
	if (allowGaps and (curThreadPath.size() > 0 or threadAlignLength != threadLength) ) {
		ssize_t loc;
		ssize_t nMismatch = 0;
		ssize_t nMatch = 0;
		ssize_t nGap = 0;
		if (alignLocations[0] != -1) {
			nGap += alignLocations[0];
			if (alignLocations[0] > threadLength)
				curAlignNumGap += threadLength;
			else
				curAlignNumGap += alignLocations[0];
		}
		for (loc = 0 ; loc < threadAlignSeq.length; loc++ ) {
			if (alignLocations[loc] != -1) {
				// There is a gap in the alignment
				if (loc > 0 and alignLocations[loc] - alignLocations[loc-1] != 1) {
					if (loc < threadLength) {
						curAlignNumGap += alignLocations[loc] - alignLocations[loc-1];
					}
					nGap += alignLocations[loc] - alignLocations[loc-1] - 1;
				}
				if (toupper(threadAlignSeq.seq[loc]) != 
						toupper(edgeAlignSeq.seq[alignLocations[loc]])) {
					// Count this against the
					// alignment score if this fits on the thread (since we 
					// aligned more than the thread).
					if (loc < threadLength) {
						curAlignNumMismatch++;
					}
					nMismatch++;
				}
				else {
					nMatch++;
				}
			}
			else {
				if (loc < threadLength) {
					curAlignNumGap++;
				}
				nGap++;
			}
		}
		// Look to see if it is worthwhile to use this local alignment.
		// If the alignment is of suitable quality, use it.

		//		if (double(nMatch) / (curAlignNumMismatch + nGap) > 0.7) {


		// Check to see if a local alignment was used to fit the sequence.
		// If it was, it is possible the alignment has gaps, and we will 
		// need to offset either the read or the edge sequence 
		// by an appropriate number of gaps when searching the next edge.

		if (curAlignNumGap  + curAlignNumMismatch < curHammingNumMismatch) {
			useLocalAlignment = 1;
		// Compute statistics about the read gaps.
		// Find the first gap.
		
			for (loc = 0; loc < threadAlignSeq.length; loc++ ) 
				if (alignLocations[loc]!= -1)
					break;
			ssize_t readStart = loc;
			ssize_t edgeStart = alignLocations[loc];
			for (loc = threadAlignSeq.length-1; loc > readStart; loc--) {
				if (alignLocations[loc] != -1)
					break;
			}
			//UNUSED// ssize_t readEnd = loc;
			ssize_t edgeEnd = alignLocations[loc];
			for (++loc; loc < threadAlignSeq.length; ++loc) {
				alignLocations[loc] = ++edgeEnd;
			}
			numReadGaps = numEdgeGaps = 0;
			for (loc = readStart; loc < threadLength; loc++ ) {
				if (alignLocations[loc] == -1) 
					numReadGaps++;
				if (loc > readStart and alignLocations[loc] - alignLocations[loc-1] != 1)
					numEdgeGaps++;
			}

			// Now some tricky logic to decide if there is a gap.
			
			// Check to see if the local alignment starts after the beginning
			// of the read.
			if (readStart > 0 or edgeStart > 0) {

				// Check two conditions to see if there was 
				// a gap at the beginning of the read.
				if (readStart - edgeStart > 0 and 
						readStart - edgeStart <= threadLength) {
					readPosGapOffset = readStart - edgeStart;
				}
				else if (edgeStart - readStart > 0 and
								 edgeStart - readStart <= threadLength ) {
					edgePosGapOffset = edgeStart - readStart;
				}
				else if (readStart == edgeStart) {
					if (readStart < threadLength) {
						curAlignNumMismatch+= readStart;
					}
					else {
						curAlignNumMismatch+= threadLength;
					}
				}
			}

			// If there are no gaps, it is not necessary to use local alignment
			// maybe assert that the number of mismatches should equal the hamming distance.
			// 
			if (numReadGaps + numEdgeGaps + readPosGapOffset + edgePosGapOffset == 0)
				useLocalAlignment = 0;

			// Compute some statistics 
			if (nGap + readStart > 0) {
				/*
				std::cout << "----------------------------------------------------------------------"
									<< std::endl;
				std::cout << "**********************************************************************"
									<< std::endl;
				std::cout << "read: " << readIndex << " edge: " << edgeIndex << std::endl;
				std::cout << "(" << curThreadPath.size() << ")  thread diff: " 
									<< curAlignNumMismatch << " " << curHammingNumMismatch << " "
									<< nMismatch << " " << threadLength << std::endl;
								PrintAlignment( threadAlignSeq, edgeAlignSeq, 0, 0, alignLocations, threadAlignSeq.length, std::cout);
				threadAlignSeq.PrintlnSeq(std::cout);
				edgeAlignSeq.PrintlnSeq(std::cout);
				if (readPosGapOffset > 0)
					std::cout << "gap offset: " << readPosGapOffset << std::endl;
				*/
			}
		}
		delete[] alignLocations;
	}

	// If local alignment was used to compute the
	// distance between two strings, use that score for this
	// part of the thread.  Otherwise, use the hamming 
	// distance.
	if (useLocalAlignment) {
		curThreadScore += curAlignNumMismatch + curAlignNumGap;
		/*		std::cout << "thread score: " << curThreadScore << " " << minThreadScore << " " << curAlignNumMismatch << " " << curAlignNumGap << std::endl;*/
	}
	else {
		curThreadScore += curHammingNumMismatch;
	}

	//
	// Only search if there is a chance the optimal threading may be 
	// found.

	if (curThreadScore < minThreadScore or minThreadScore == -1) {
		std::string threadSeq;
		ssize_t ep;
		for (ep = edgePos ; ep <edgePos + threadLength; ep++ ) {
			threadSeq += graph.edges[edgeIndex].seq.seq[ep];
		}

		curThreadPath.push_back(ThreadPathInterval(edgeIndex, threadLength, edgePos));
		ssize_t res = 0;

		if (seqPos == curSeqLength) {
			//			std::cout << readIndex << " " << curThreadScore << " " << minThreadScore << std::endl;
			// This is now the global min path
			minThreadScore = curThreadScore;
			minThreadPath  = curThreadPath;
			res = 1;
			
			// We'll keep the top "scorListLength" thread scores.
			// place this one on the podium if it is low enough.
			//UNUSED// ssize_t m;
			ssize_t ms, minScoreIndex;
			for (ms = 0; ms < scoreListLength; ms++ ) {
				if (scoreList[ms] == -1 or
						scoreList[ms] >= minThreadScore) {
					break;
				}
			}
			if (ms < scoreListLength) {
				minScoreIndex = ms;
				// Make room for the new next best score
				for (ms = scoreListLength-1; ms > minScoreIndex; ms--) 
					scoreList[ms] = scoreList[ms-1];
				scoreList[minScoreIndex] = minThreadScore;
			}
		} else {
			// The entire read has not been threaded. Keep going.
			ssize_t outEdge,outEdgeIndex;
			if (readPosGapOffset + numReadGaps > 0) {
				std::cout << "Advancing an extra: " << readPosGapOffset + numReadGaps << std::endl;
				std::cout << "Stepping back: " << edgePosGapOffset + numEdgeGaps << std::endl;
				printNext = 1;
			}
			curSeq = &curSeq[seqPos + readPosGapOffset + numReadGaps];
			//			assert(seqPos != 0);
			if (readPosGapOffset + numReadGaps > 0 and
					curThreadPath.size() > 1)  {
				std::cout << " skipping on internal edge!!" << std::endl;
			}
			if (seqPos > 0) {
				curSeqLength -= (seqPos + readPosGapOffset + numReadGaps);
				// If curseqlength was 0 we would have threaded 
				// the entire sequence and wouldn't be here.
				//assert(curSeqLength >= 0);
				if (curSeqLength >= 0) {
					for (outEdgeIndex = graph.vertices[destVertex].FirstOut();
							 outEdgeIndex < graph.vertices[destVertex].EndOut();
							 outEdgeIndex = graph.vertices[destVertex].NextOut(outEdgeIndex)) {
						outEdge = graph.vertices[destVertex].out[outEdgeIndex];
						if (depth < maxDepth - 1) 
							res |= ThreadRecursively(graph, outEdge, edgePosGapOffset + numEdgeGaps, 
																			 curSeq, curSeqLength,
																			 curThreadPath, curThreadScore,
																			 minThreadPath, minThreadScore, maxScore, maxDepth,
																			 perfectRead,
																			 scoreList, scoreListLength, depth+ 1, scoreMat, 
																			 printNext, readIndex, allowGaps);
					}
				}
			}
		}
		curThreadPath.pop_back();
		return res;
	}
	else {
		return 0;
	}
}
											

ssize_t StorePathSequences(IntervalGraph &graph, 
											 ssize_t edgeIndex,ssize_t edgePos, 
											 std::string curSequence,
											 ssize_t searchLength, 
											 std::vector<std::string> &sequences) {
	std::string newSequence;
	assert(searchLength >= 0);
	if (searchLength == 0) return sequences.size(); 
	if (edgePos + searchLength < graph.edges[edgeIndex].length ) {
		newSequence = std::string((const char*) (graph.edges[edgeIndex].seq.seq + edgePos),
															searchLength);
		sequences.push_back(curSequence + newSequence);
	}
	else {
		// Not done searching graph for sequence of length 'searchLength'

		newSequence = std::string((const char*) (graph.edges[edgeIndex].seq.seq + edgePos), 
															graph.edges[edgeIndex].length - edgePos);
		// Append the substring
		curSequence = curSequence + newSequence;
				
		ssize_t outEdgeIndex, outEdge;
		ssize_t dest;
		dest = graph.edges[edgeIndex].dest;
		for (outEdgeIndex = graph.vertices[dest].FirstOut();
				 outEdgeIndex < graph.vertices[dest].EndOut();
				 outEdgeIndex = graph.vertices[dest].NextOut(outEdgeIndex)) {
			outEdge = graph.vertices[dest].out[outEdgeIndex];

			StorePathSequences(graph, outEdge, 
												 graph.vertices[dest].vertexSize, curSequence, 
												 searchLength - (graph.edges[edgeIndex].length - edgePos), 
												 sequences);
		}
	}
	return sequences.size();
}

void AppendReverseComplements(DNASequenceList &sequences) {
  ssize_t seq;
  ssize_t numSeq = sequences.size();
  SimpleSequence simpleSeq;
  simpleSeq.seq = NULL;
  sequences.resize(numSeq*2);
	// First stratify the sequences
	for (seq = numSeq-1; seq > 0; seq--) {
		sequences[seq*2].seq = sequences[seq].seq;
		sequences[seq*2].length = sequences[seq].length;
	}
	// fill in the reverse compliments between
  for (seq = 0; seq < numSeq; seq++ ) {
    MakeRC((char*) sequences[seq*2].seq, sequences[seq*2].length, sequences[seq*2+1].seq);
    sequences[seq*2+1].length = sequences[seq*2].length;
  }
}
