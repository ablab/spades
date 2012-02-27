/***************************************************************************
 * Title:          GFFFile.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "GFFFile.h"


ssize_t GFFFile::ParseGFFFile(std::istream &in) {
  GFFEntry *entry;
  entry = new GFFEntry;
  while (entry->ParsePALSLine(in)) {
    entries.push_back(entry);
    entry = new GFFEntry;
  }
  // created one too many entries
  delete entry;

	return 1; // success
}


void GFFFile::ToLAVFile(LAVFile &lavFile, 
		       std::string &refName, ssize_t refLength, 
		       std::string &qryName, ssize_t qryLength) {
  // For now only process one sequence, in + and - strands.
  lavFile.blastzOpts = LAVFile::standardBlastzOpts;
  LAVAlignedContig *forwardContig, *reverseContig;
  // Process the + strand
  forwardContig = new LAVAlignedContig;
  reverseContig = new LAVAlignedContig;
  lavFile.alignments.push_back(forwardContig);
  lavFile.alignments.push_back(reverseContig);
  if (entries.size() > 0) {
    forwardContig->refContig.sequenceName = "\"" + entries[0]->refSeq + "\"";
    forwardContig->refContig.strand = 0;
    forwardContig->refContig.start = 1;
    forwardContig->refContig.end = 1 + refLength;
    forwardContig->refContigName = "\">" + refName + "\"";
    forwardContig->qryContig.sequenceName = "\"" + entries[0]->qrySeq + "\"";
    forwardContig->qryContig.strand = 0;
    forwardContig->qryContig.end = 1 + qryLength;
    forwardContig->qryContigName = "\">" + qryName + "\"";

    reverseContig->refContig.sequenceName = "\"" + entries[0]->refSeq + "\"";
    reverseContig->refContig.strand = 0;
    reverseContig->refContig.start = 1;
    reverseContig->refContig.end = 1 + refLength;
    reverseContig->refContigName = "\">" + refName + "\"";
    reverseContig->qryContig.sequenceName = "\"" + entries[0]->qrySeq + "-\"";
    reverseContig->qryContig.strand = 1;
    reverseContig->qryContig.end = 1 + qryLength;
    reverseContig->qryContigName = "\">" + qryName + " (reverse complement)\"";
  }
  ssize_t i;
  for (i = 0; i < entries.size(); i++) {
    if (entries[i]->strand == '+') {
      LAVBlock *block = new LAVBlock;
      block->refBegin  = entries[i]->refStart;
      block->refEnd    = entries[i]->refEnd;
      block->qryBegin  = entries[i]->qryStart;
      block->qryEnd    = entries[i]->qryEnd;
      block->refALBegin.push_back(entries[i]->refStart);
      block->refALEnd.push_back(entries[i]->refEnd);
      block->qryALBegin.push_back(entries[i]->qryStart);
      block->qryALEnd.push_back(entries[i]->qryEnd);
      block->alIdentity.push_back((ssize_t)(1-entries[i]->errorRatio)*100);
      forwardContig->push_back(block);
    }
  }
  for (i = 0; i < entries.size(); i++) {
    if (entries[i]->strand == '-') {
      LAVBlock *block = new LAVBlock;
      block->refBegin  = entries[i]->refStart;
      block->refEnd    = entries[i]->refEnd;
      block->qryBegin  = qryLength - (entries[i]->qryEnd + 1);
      block->qryEnd    = qryLength - (entries[i]->qryStart + 1);
      block->enumeration = 0;
      block->refALBegin.push_back(entries[i]->refStart);
      block->refALEnd.push_back(entries[i]->refEnd);
      block->qryALBegin.push_back(qryLength - (entries[i]->qryEnd + 1));
      block->qryALEnd.push_back(qryLength - (entries[i]->qryStart + 1));
      block->alIdentity.push_back((ssize_t)(1-entries[i]->errorRatio)*100);
      reverseContig->push_back(block);
    }
  }
}
