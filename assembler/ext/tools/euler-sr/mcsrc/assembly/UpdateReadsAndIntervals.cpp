/***************************************************************************
 * Title:          UpdateReadsAndIntervals.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DeBruijnGraph.h"
#include "IntervalGraph.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "SimpleSequence.h"
#include <string>
#include "IntegralTupleStatic.h"


void PrintUsage() {
  std::cout << "updateReadsAndIntervals graphFile intvFile seqFile newIntvFile newSeqFile" << std::endl;
}

int main(int argc, char* argv[]) {

  std::string graphFileName, intvFileName, readFileName;
  std::string intvOutName, readsOutName;

  int argi = 1;
  if (argc < 5) {
    PrintUsage();
    return 1;
  }
  graphFileName = argv[argi++];
  intvFileName  = argv[argi++];
  readFileName  = argv[argi++];
  intvOutName   = argv[argi++];
  readsOutName  = argv[argi++];
  
  std::cout << "reading interval graph" << std::endl;
  // Get the graph, including the read-intervals
  IntervalGraph graph;
  graph.ReadIntervalGraph(graphFileName, intvFileName);

  // Get the reads.  Chances are there are more reads 
  // than there are read intervals, otherwise the intervals
  // wouldn't need to be updated.
  std::cout << "reading reads " << std::endl;
  SimpleSequenceList seqList;
  ReadSimpleSequences(readFileName, seqList);
  AppendReverseComplements(seqList);

  std::vector<char> representedReads;
  representedReads.resize(seqList.size());
  ssize_t r;
  for (r = 0; r < representedReads.size(); r++ ) {
    representedReads[r] = 0;
  }
  std::cout << "determining what reads are represented " << std::endl;

  // tally which reads are used in this graph
  ssize_t e, i;
  for (e = 0; e < graph.edges.size(); e++ ) {
    for (i = 0; i < graph.edges[e].intervals->size(); i++) {
      representedReads[(*graph.edges[e].intervals)[i].read] = 1;
    }
  }

  std::cout << "outputting reads " << std::endl;
  // Output the reads that are used in the graph
  std::ofstream readsOut;
  openck(readsOutName, readsOut, std::ios::out);
  DNASequence tmpSeq;
  std::string blank;
  for (r = 0; r < seqList.size()/2; r++ ) {
    if (representedReads[r]) {
      blank = "";
      tmpSeq.titlestream->str("");
      tmpSeq.seq     = seqList[r].seq;
      *tmpSeq.titlestream << r;
      tmpSeq.length  = seqList[r].length;
      tmpSeq.PrintSeq(readsOut);
      readsOut << std::endl;
    }
  }
  readsOut.close();

  ssize_t numRepresentedReads = 0;
  for (r = 0; r < representedReads.size(); r++ ) {
    numRepresentedReads += representedReads[r];
  }
  std::vector<ssize_t> readIndices;
  readIndices.resize(numRepresentedReads);
  ssize_t readIndex = 0;
  for (r = 0; r < representedReads.size(); r++ ) {
    if (representedReads[r]) {
      readIndices[readIndex++] = r;
    }
  }

  // Just in case
  assert(readIndex == numRepresentedReads);

  // Update the indices of the read-intervals so that 
  // they correspond to the fewer number of reads
  std::cout << "updating interval indices " << std::endl;
  ssize_t oldReadIndex, newReadIndex;
  for (e = 0; e < graph.edges.size(); e++ ) {
    for (i = 0; i < graph.edges[e].intervals->size(); i++) {
      oldReadIndex = (*graph.edges[e].intervals)[i].read;
      newReadIndex = std::binary_search(readIndices.begin(), readIndices.end(), oldReadIndex);
      (*graph.edges[e].intervals)[i].read = newReadIndex;
    }
  }

  std::cout << "printing read intervals" << std::endl;
  PrintReadIntervals(graph.edges, intvOutName);

  return 0;
}
