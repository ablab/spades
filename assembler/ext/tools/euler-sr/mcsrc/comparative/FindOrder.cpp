/***************************************************************************
 * Title:          FindOrder.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
// std includes
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <map>
#include <set>
#include <algorithm>
#include <unistd.h>
#include <string.h>

// 3rd party includes
#include "mysql/mysql.h"

// my stuff
#include "utils.h"
#include "lav/LAVReader.h"
#include "lav/LAVTable.h"
#include "lav/LAVUtils.h"
#include "lav/LAVAlignedContig.h"

#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"

int main(int argc, char* argv[]) {
  std::string alignFileName, sequenceFileName, newSequenceFileName;
  std::vector<ssize_t> orientations;
  if (argc < 4) {
    std::cout << "usage: ordient alignFile unorientedSeqIn orientedSeqOut " 
	      << std::endl;
    std::cout << "alignFile is an alignment of a reference sequence to the " 
	      << std::endl
	      << "  unoriented sequence. " << std::endl;
    exit(0);
  }

  alignFileName = argv[1];
  sequenceFileName = argv[2];
  newSequenceFileName = argv[3];
  ssize_t printOrder = 0;
  if (argc > 3) {
    int argi;
    if (strcmp(argv[argi], "-printorder")) {
      printOrder = 1;
    }
  }
  LAVFile alignment;
  LAVReader::ReadLAVFile(alignFileName, alignment);

  // Determine the orientation of each contig.
  ssize_t i, j, k, o;
  LAVAlignedContig *alignedContig;

  std::set<ssize_t> unusedAlignments;
  std::vector<ssize_t> usedAlignments;
  for (i = 0; i < alignment.size(); i++) {
    unusedAlignments.insert(i);
    orientations.push_back(-1);
  }

  
  ssize_t orientation;
  LAVAlignedContig *compAlign;
  ssize_t compIndex;
  ssize_t len;
  for (i = 0; i < alignment.size(); i++) {
    len = alignment.alignments[i]->alignments.size();
  }
  for (i = 0; i < alignment.size(); i++) { 
    if (unusedAlignments.find(i) == unusedAlignments.end())
      continue;
    
    compIndex = FindComplementaryAlignment(alignment.alignments[i], alignment, compAlign);
    if (alignment.alignments[i]->qryContig.strand == 0)
      orientation = DetermineOrientation(alignment.alignments[i], compAlign);
    else
      orientation = DetermineOrientation(compAlign, alignment.alignments[i]);

    if (compIndex >= 0) {
      // There is more than one orientation of this contig involved in the alignment.
      // Remember what the orientation of the contig is.
      // Skip thse contigs in the loop

      // Only consider the alignment of the contig that is for the order of the contig.
      if (orientation == 0) {
	usedAlignments.push_back(i); 
      }
      else {
	usedAlignments.push_back(compIndex);
      }

      // Don't consider either for alignment again.
      unusedAlignments.erase(i);
      unusedAlignments.erase(compIndex);
      // record the orientation of this contig.
      orientations[i] = orientation; //.push_back(orientation);
      orientations[compIndex] = orientation;
    }
    else {
      // There is only one orientation of this contig aligned.
      usedAlignments.push_back(i);
      //      orientations.push_back(orientation);
      orientations[i] = orientation;
    }
  }
  ssize_t index, index2;
  ssize_t skipIndex;
  // Find duplicate aligned locations.
  std::vector<ssize_t> skippedIndices;
  std::vector<std::string> skippedNames;
  std::vector<ssize_t>::iterator iit, jit;

  for (iit = usedAlignments.begin(); iit != usedAlignments.end(); ++iit) {
    jit = iit;
    ++jit;
    for (; jit != usedAlignments.end();) {
      index = *iit;
      index2 = *jit;
      if (alignment.alignments[index]->alignments[0]->refBegin ==
	  alignment.alignments[index2]->alignments[0]->refBegin) {
	// want to skip index 2
	skippedIndices.push_back(index2);
        skippedNames.push_back(alignment.alignments[index2]->qryContigName);
	jit = usedAlignments.erase(jit);
      }
      else {
	++jit;
      }
    }
  }
  ssize_t s;
  std::cout << "skipped: " << std::endl;
  for (s = 0; s < skippedIndices.size(); s++ ) {
    std::cout << skippedIndices[s] << " " << skippedNames[s] << std::endl;
  }
  std::cout << std::endl;
  std::vector<ssize_t> order;    // will contain the order of the indicesin usedAlignments.
  std::vector<ssize_t> locations; // the start locations of the alignments
  // Determine the order of the used contigs.
  for (i =0; i < usedAlignments.size(); i++) {
    index = usedAlignments[i];
    skipIndex = 0;
    for (k = 0; k < skippedIndices.size(); k++) {
      if (index == skippedIndices[k])
	skipIndex = 0;
    }
    if (! skipIndex ) 
      locations.push_back(alignment.alignments[index]->alignments[0]->refBegin);
  }



  sort(locations.begin(), locations.end());

  for ( i = 0; i < usedAlignments.size(); i++) {
    index = usedAlignments[i];
    for ( j = 0; j < usedAlignments.size(); j++) {
      index2 = usedAlignments[j];
      skipIndex = 0;
      for (k = 0; k < skippedIndices.size(); k++ ) {
	if (index2 == skippedIndices[k])
	  skipIndex = 0;
      }
      if (! skipIndex && locations[i] == alignment.alignments[index2]->alignments[0]->refBegin)
	order.push_back(j);
    }
  }
  
  std::vector<DNASequence*> contigs;
  DNASequence *sequenceP;
  std::ifstream seqIn;
  openck(sequenceFileName, seqIn);
  ssize_t length = 0;
  SeqReader::MaskRepeats();
  while (true) {
    sequenceP = new DNASequence; 
    if (SeqReader::GetSeq(seqIn, *sequenceP)) {
      skipIndex = 0;
      for (i = 0; i < skippedNames.size(); i++ ) {
	if (skippedNames[i].find(sequenceP->namestr) != skippedNames[i].npos)
	  skipIndex = 0;
      }
      if (! skipIndex ) {
	contigs.push_back(sequenceP);
	length += sequenceP->length;
      }
      else {
	delete sequenceP;
      }
    }
    else {
      delete sequenceP;
      break;
    }
  }
  
  DNASequence orderedSequence(length);
  DNASequence revComp;

  ssize_t pos = 0;

  // Now copy join the contigs together in the correct order and orientation.
  for (i = 0; i < usedAlignments.size(); i++) {
    index = usedAlignments[order[i]];
    skipIndex = 0;
    for (k = 0; k < skippedIndices.size(); k++) {
      if (skippedIndices[k] == index)
	skipIndex = 1;
    }
    if (! skipIndex ){
      for (j = 0; j < contigs.size(); j++) {
	if (alignment.alignments[index]->qryContigName.find(contigs[j]->namestr) !=
	    alignment.alignments[index]->qryContigName.npos) {
	  // We've found the contig that belongs in this position.
	  /*
	  std::cout << "storing contig " << j << " i:" << index << " in orientation : " 
		    << orientations[index] << std::endl;
	  */
	  if (printOrder == 0) {
	    if (orientations[index] == 0) {
	      memcpy(&(orderedSequence.seq[pos]), &(contigs[j]->seq[0]), contigs[j]->length);
	    }
	    else {
	      MakeRC(*contigs[j], revComp);
	      memcpy(&(orderedSequence.seq[pos]), &(revComp.seq[0]), revComp.length);	  	  
	      delete revComp.seq;
	      revComp.seq = NULL;
	    }
	    pos += contigs[j]->length;
	  }
	  else {
	    std::cout << contigs[j]->namestr << std::endl;
	  }
	  break;
	}
      }
    }
  }
  // reset the length
  if (pos < orderedSequence.length) {
    DNASequence tempSeq;
    tempSeq = orderedSequence;
    orderedSequence.Reset(pos);
    memcpy(&(orderedSequence.seq[0]), &(tempSeq.seq[0]), pos);
    delete[] tempSeq.seq;
  }
    
  orderedSequence.StoreName(contigs[0]->namestr);
  std::ofstream newSeq;
  openck(newSequenceFileName, newSeq, std::ios::out);
  orderedSequence.PrintSeq(newSeq, 50);
  newSeq << std::endl;
}


