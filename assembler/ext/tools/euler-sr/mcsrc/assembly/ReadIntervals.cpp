/***************************************************************************
 * Title:          ReadIntervals.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "compatibility.h"
#include "ReadIntervals.h"
#include "BufferedSeqReader.h"

int EdgeIntervalList::tupleSize = 0;


			
ssize_t CountEdges(std::string &readIntervalFileName) {
  std::ifstream readIntervalIn;
  openck(readIntervalFileName, readIntervalIn, std::ios::in);
  std::string word;
  ssize_t numEdges = 0;
  while(readIntervalIn) {
    if (! (readIntervalIn >> word)) break;
    if (word == "EDGE")
      ++numEdges;
  }
  readIntervalIn.close();
  return numEdges;
}

void StoreEdgeInterval(ReadPositions &overlaps, _INT_ overlapSize,
											 SimpleSequenceList &edges,
											 SimpleSequence &sequence,
											 EdgeIntervalList &edgeInterval,
											 std::vector<ssize_t> &mult) {
  // use dbEdgeIndex to index into the list of de Bruijn edges that 
  // are generated from the reads.  If an edge is covered by more than
  // one read, it will be present multiple times in the list 'edges'.
  // we need to increment an interval for every edge in the list
  //UNUSED// ssize_t edgeIndex;

  //UNUSED// ssize_t s;
  ssize_t overlapIndex;
  ssize_t nextN;

	//    std::cout << "processing sequence " << s << std::endl;
	nextN = -1;
	ssize_t p;
	if (sequence.length < overlapSize)
		return;
	for (p = std::min((_SSZT_) (overlapSize-1), sequence.length); p >= 0; p--) { // TODO: remove cast
		if (numeric_nuc_index[sequence.seq[p]] >= 4) {
			nextN = p;
			break;
		}
	}
	for (p = 0; p < sequence.length - overlapSize + 1; p++ ) {
		if (numeric_nuc_index[sequence.seq[p + overlapSize - 1]] >= 4) 
			nextN = p + overlapSize - 1;
		if (p > nextN) {
			overlapIndex = LocateTuple(edges, overlaps, overlapSize, (char*) &(sequence.seq[p]));
			if (overlapIndex < 0) {
				DNASequence tmpSeq;
				tmpSeq.seq = &(sequence.seq[p]);
				tmpSeq.length = overlapSize;
				tmpSeq.PrintSeq(std::cout);
				std::cout << std::endl;
				std::cout << "error, should have found the read sequence in an overlap!!" << std::endl;
				exit(1);
			}
			/*      std::cout << "adding interval: " << overlaps[overlapIndex].read << " " 
							<< overlaps[overlapIndex].pos << std::endl;
			*/
			edgeInterval.IncrementEdgeInterval(overlaps[overlapIndex].read, 
																				 overlaps[overlapIndex].pos, 
																				 p,
																				 mult[overlaps[overlapIndex].read]);
		}
	}
}



ssize_t IncrementReadIntervalList(ssize_t edge, ssize_t edgePos, ssize_t read, ssize_t prevPos, ssize_t readPos,
															ReadIntervalList &readIntervals, int tupleSize,
															ssize_t prevEdge) {
	//
	// Input: an edge and a position on the edge of a read interval that may
	//       follow other read intervals.  If the read interval is the adjacent position
	//       on the same edge as the previous read interval, the previous read interval
	//       is simply incremented as length.  If it is a new edge, or a gap on the same
	//       edge, a new interval is created.
	// Output: 
	//       return value: 0 if the read interval length is grown.
	//                     1 if a new read interval si created.
	//

	ssize_t cur = readIntervals.size() - 1;
	ssize_t gap = readPos - prevPos;
	if (cur >= 0 and
			prevEdge == edge and
			readIntervals[cur].edgePos + 
			readIntervals[cur].length - tupleSize == edgePos - gap) {
		readIntervals[cur].length += gap;
		return 0;
	}
	else {
		++cur;
		readIntervals.push_back(ReadInterval(edge, read, readPos, edgePos, tupleSize, cur));
		return 1;
	}
}


void StoreAllReadIntervals(ReadPositions &dbEdges, _INT_ dbEdgeSize, 
													 SimpleSequenceList &edges, 
													 std::string &sequenceListFileName, 
													 std::vector<BareReadIntervalList> &edgeReadIntervals,
													 PathIntervalList &paths,
													 std::vector<ssize_t> &pathLengths, _INT_ skipGapped) {

	edgeReadIntervals.resize(edges.size());

	BufferedSeqReader<1000> seqReader;
	//std::ifstream readsIn;
	//	openck(sequenceListFileName, readsIn, std::ios::in);
	//	SeqReader seqReader(&readsIn);
	seqReader.Init(sequenceListFileName);
	DNASequence read, readRC;
	PathIntervalList path;
	ssize_t readIndex;
	DNASequence* readPair[2], *readPtr;
	readPair[0] = &read;
	readPair[1] = &readRC;
	readIndex = 0;
	ssize_t rn =0;
	
	while (seqReader.GetSeq(read)) {

		if (skipGapped) {
			ssize_t seqIsGapped = 0;
			ssize_t curPos;
			for (curPos = 0; curPos < read.length; curPos++ ){ 
				if (unmasked_nuc_index[read.seq[curPos]] >= 4) {
					seqIsGapped = 1;
					break;
				}
			}
			if (seqIsGapped) {
				continue;
			}
		}

		MakeRC(read, readRC);
		//UNUSED// ssize_t edgeIndex;
		//UNUSED+// ssize_t s;
		ssize_t r ;
		ssize_t overlapIndex;
		ssize_t nextN;
		ReadIntervalList readPath;
		ssize_t pathPos;
		//		std::cout << "read: " << rn << " " << read.namestr << std::endl;
		//read.PrintlnSeq(std::cout);
		//		std::cout << read.namestr << std::endl;
		++rn;
		PrintStatus(rn, 100000);
		for (r = 0; r < 2; r++) {
			//    std::cout << "processing sequence " << s << std::endl;
			readPath.clear();
			readPtr = readPair[r];
			nextN   = -1;
			ssize_t p;
			ssize_t prevEdge = -1;
			ssize_t fullReadMapped = 1;

			if (readPtr->length < dbEdgeSize) {
				fullReadMapped = 0;
			}
			else {
				for (p = std::min((_SSZT_) (dbEdgeSize-1), readPtr->length); p >= 0; p--) { // TODO: remove cast
					if (numeric_nuc_index[read.seq[p]] >= 4) {
						nextN = p;
						fullReadMapped = 0;
						break;
					}
				}
				ssize_t readPos;
				std::vector<ssize_t> edgeList;
				std::vector<ssize_t> readList, posList, prevEdgeList;
				for (readPos = 0; 
						 readPos < readPtr->length - dbEdgeSize + 1 and fullReadMapped; readPos++ ) {
					if (numeric_nuc_index[readPtr->seq[readPos + dbEdgeSize - 1]] >= 4) {
						nextN = readPos + dbEdgeSize - 1;
						fullReadMapped = 0;
					}
					if (readPos > nextN) {
						overlapIndex = LocateTuple(edges, dbEdges, dbEdgeSize, (char*) &(readPtr->seq[readPos]));
						if (overlapIndex < 0) {
							/*
							DNASequence tmpSeq;
							tmpSeq.seq = &(readPtr->seq[readPos]);
							tmpSeq.length = dbEdgeSize;
							tmpSeq.PrintSeq(std::cout);
							std::cout << std::endl;
							std::cout << "error, should have found the read sequence in an overlap!!" << std::endl;
							exit(1);
							*/
							fullReadMapped = 0;
						}
						/*      
										std::cout << "adding interval: " << overlaps[overlapIndex].read << " " 
										<< overlaps[overlapIndex].pos << std::endl;
						*/
						readList.push_back(dbEdges[overlapIndex].read);
						posList.push_back(dbEdges[overlapIndex].pos);
						prevEdgeList.push_back(prevEdge);
															 
						prevEdge = dbEdges[overlapIndex].read;
					}
				}
				
				if (fullReadMapped) {
					ssize_t prevPos = -1;
					for (readPos = 0; readPos < readPtr->length - dbEdgeSize + 1; readPos++ ) {
						if (IncrementReadIntervalList(readList[readPos], posList[readPos], 
																					readIndex, prevPos, readPos, readPath, dbEdgeSize,
																					prevEdgeList[readPos])) {
							edgeList.push_back(readList[readPos]);
						}
					}
				
				

			
					// Now process the read interval list.
					paths.push_back(new PathInterval[readPath.size()]);
					pathLengths.push_back(readPath.size());
					for (pathPos = 0; pathPos < readPath.size(); pathPos++) {
						ssize_t edgeIndex = edgeList[pathPos]; //readPath[pathPos].edge;
						edgeReadIntervals[edgeIndex].push_back(BareReadInterval(readIndex, 
																																		readPath[pathPos].readPos,
																																		readPath[pathPos].edgePos,
																																		readPath[pathPos].length));
						paths[readIndex][pathPos].edge  = edgeIndex;
						paths[readIndex][pathPos].index = readPath[pathPos].readPos;
						/*				int last = edgeReadIntervals[edgeIndex].size()-1;
											assert(edgeReadIntervals[edgeIndex][last].read == readIndex);
						*/
					}
				}
			}
			if (!fullReadMapped) {
				// We couldn't map the full read either due to N's, the sequence wasn't 
				// stored, or the sequence was too short.  Later modify this so that 
				// a warning of some sort is printed when the read cannot be mapped.
				paths.push_back(NULL);
				pathLengths.push_back(0);
			}
			++readIndex;
		}
	}
}

void PathReadPosToIntervalIndex(std::vector<BareReadIntervalList> &edgeReadIntervals,
																std::vector<std::vector<ssize_t> > &edgeIntervalIndices,
																PathIntervalList &paths,
																std::vector<ssize_t> &pathLengths) {
	ssize_t p, pi;
	//UNUSED// ssize_t ei;
	ssize_t e;
	//UNUSED// ssize_t print = 0;
	for (p = 0; p < paths.size(); p++) {
		for (pi = 0; pi < pathLengths[p]; pi++) {
			e = paths[p][pi].edge;

			// Do a binary search for the read index.
			ssize_t start = 0;
			ssize_t end   = edgeIntervalIndices[e].size();
			ssize_t cur   = end/2;
			ssize_t curIndex = edgeIntervalIndices[e][cur];
			while(start+1 < end and edgeReadIntervals[e][curIndex].read != p) {
				if (p < edgeReadIntervals[e][curIndex].read ) {
					end = cur;
				}
				else {
					start = cur;
				}
				cur = (start + end) / 2;
				curIndex = edgeIntervalIndices[e][cur];
			}
			if (edgeReadIntervals[e][curIndex].read != p) {
				std::cout << "ERROR, searched for a read index that should have been found." 
									<< std::endl;
				exit(1);
			}
			// Rewind to the first instance of this read on this edge
			// or the match, if it is here.
			ssize_t back = cur;
			 
			while(// cur isn't the match
						edgeReadIntervals[e][edgeIntervalIndices[e][back]].readPos != paths[p][pi].index and 
						// still is space
						back > 0 and  
						// prev still read
						edgeReadIntervals[e][edgeIntervalIndices[e][back-1]].read == p )
				back--;

			if (edgeReadIntervals[e][edgeIntervalIndices[e][back]].readPos == 
					paths[p][pi].index) {
				paths[p][pi].index = edgeIntervalIndices[e][back];
			}
			else {
				// already searched cu
				cur++;
				while (cur < edgeIntervalIndices[e].size() and
							 edgeReadIntervals[e][edgeIntervalIndices[e][cur]].read == p and
							 edgeReadIntervals[e][edgeIntervalIndices[e][cur]].readPos != paths[p][pi].index) {
					cur++;
				}
				if (cur == edgeIntervalIndices[e].size() or 
						edgeReadIntervals[e][edgeIntervalIndices[e][cur]].read != p) {
					std::cout << "ERROR, searched for a read index that should have been found."
										<< std::endl;
					exit(1);
				}
				paths[p][pi].index = edgeIntervalIndices[e][cur];
			}
		}
	}
}


void StoreEdgeIntervals(ReadPositions &overlaps, _INT_ overlapSize,
												SimpleSequenceList &edges,
												SimpleSequenceList &sequences,
												EdgeIntervalListList &edgeIntervals,
												PathIntervalList &paths,
												PathLengthList &pathLengths,
												std::vector<ssize_t> &mult) {

  edgeIntervals.resize(sequences.size());
  mult.resize(edges.size());
	paths.resize(sequences.size());
	pathLengths.resize(sequences.size());
	std::fill(pathLengths.begin(), pathLengths.end(), 0);
  //UNUSED+// ssize_t  v,p;
  ssize_t e ; 
	ssize_t s;
	for (s = 0; s < sequences.size(); s++ ) {
		StoreEdgeInterval(overlaps, overlapSize, edges, sequences[s], edgeIntervals[s], mult);
		paths[s] = new PathInterval[edgeIntervals[s].edgeIntervals.size()];
		pathLengths[s] = edgeIntervals[s].edgeIntervals.size();
		for (e = 0; e < edgeIntervals[s].edgeIntervals.size(); e++ ) {
			paths[s][e].edge = edgeIntervals[s].edgeIntervals[e].edge;
		}
	}
}


ssize_t CountEdgeIntervals(EdgeIntervalListList &edgeIntervals) {

  ssize_t numEdgeIntervals = 0;
  ssize_t e;
  for (e = 0; e < edgeIntervals.size(); e++ ) {
    numEdgeIntervals += edgeIntervals[e].edgeIntervals.size();
  }
  return numEdgeIntervals;
}

void StoreReadIntervals(EdgeIntervalListList &edgeIntervals,
												ReadIntervalList &readIntervals) {
  ssize_t e, i;
  ssize_t intv;
  intv = 0;
  ssize_t readPos;
  for (e = 0; e < edgeIntervals.size(); e++ ) {
    readPos = 0;
    for (i = 0; i < edgeIntervals[e].edgeIntervals.size(); i++ ) {
      readIntervals[intv]      = edgeIntervals[e].edgeIntervals[i];
      readIntervals[intv].read = e;
			readIntervals[intv].pathPos = i;
      /*
				std::cout << intv << " " << readIntervals[intv].read << " " 
				<< readIntervals[intv].edgePos << " " 
				<< readIntervals[intv].readPos << " "
				<< readIntervals[intv].length << std::endl;
      */
      intv++;
    }
  }
  std::cout << "stored " << intv << " read intervals " << std::endl;
}


void SortReadIntervalsByReadPos(ReadIntervalList &list) {
	CompareReadIntervalsByReadPos comp;
  std::sort(list.begin(), list.end(), comp);
}

void SortReadIntervalIndicesByRead(BareReadIntervalList &list, 
																	std::vector<ssize_t> &listIndices) {
	CompareBareReadIntervalIndicesByRead comp;
	comp.readIntervalPtr = &list;
	std::sort(listIndices.begin(), listIndices.end(), comp);
}

void SortReadIntervals(ReadIntervalList &list) {
  CompareReadIntervals comp;
  std::sort(list.begin(), list.end(), comp);
}

void SortBareReadIntervalsByReadPos(BareReadIntervalList &list) {
  CompareBareReadIntervalsByReadPos comp;
	std::sort(list.begin(), list.end(), comp);
}

void SortBareReadIntervalsByEdgePos(BareReadIntervalList &list) {
  CompareBareReadIntervalsByEdgePos comp;
	std::sort(list.begin(), list.end(), comp);
}

void EdgeToReadIntervals(EdgeIntervalListList &edgeIntervals,
												 ReadIntervalList &readIntervals) {

  ssize_t numEdgeIntervals = CountEdgeIntervals(edgeIntervals);
  readIntervals.resize(numEdgeIntervals);
  StoreReadIntervals(edgeIntervals, readIntervals);

  // sort read intervals according to their positions in the graph
  // then by their positions along each edge
  // then by their read occurrence
  // finally by their occurence in each read (if a read maps to an edge more 
  // than once)
  SortReadIntervals(readIntervals);
}
