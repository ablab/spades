/***************************************************************************
 * Title:          Params.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _PARAMS_H_
#define _PARAMS_H_
#include <string>
#include <ios>
#include <ostream>
#include <istream>
#include <fstream>

class T_Params {
public:
  ssize_t hashStep;
  ssize_t hashLength;
  ssize_t wordLength;
  ssize_t doIndel;
  ssize_t neighborDist;
  ssize_t alignStrips;
  ssize_t alignGaps;
  double lisGapPenalty;
  double lisInsertPenalty;
  ssize_t minGapDist, maxGapDist;
  ssize_t distThreshold;
  ssize_t minHashLen;
  ssize_t productSize;
  ssize_t lisThreshold;
  ssize_t lisWindow;
  ssize_t mergeThreshold;
  double ratioThreshold;
  ssize_t doCheckpoint;
  ssize_t loadCheckpoint;
  ssize_t outputTuple;
  ssize_t outputGlue;
  ssize_t maskRepeats;
  std::string checkpointFile;
  T_Params() {
    hashLength     = 20;
    wordLength     = 14;
    doIndel        = 0;
    neighborDist   = 0;
    hashStep       = 2;
    alignStrips    = 0;
    alignGaps      = 0;
    minGapDist     = 3;
    maxGapDist     = 1000000;
    distThreshold  = 0;
    minHashLen     = hashLength;
    productSize    = 50000;
    lisWindow      = 0;
    lisThreshold   = 0;
    mergeThreshold = 1;
    ratioThreshold = 0.1;
    doCheckpoint   = 0;
    loadCheckpoint = 0;
    checkpointFile = "";
    lisGapPenalty  = 0.1;
    lisInsertPenalty = 0.01;
    outputTuple    = 0;
    outputGlue     = 0;
    maskRepeats    = 0;
  }
  T_Params &operator=(const T_Params &p) {
		if (this != &p) {
			hashLength    = p.hashLength;   
			wordLength    = p.wordLength;  
			doIndel       = p.doIndel;      
			neighborDist  = p.neighborDist;
			hashStep      = p.hashStep;     
			alignStrips   = p.alignStrips;  
			alignGaps     = p.alignGaps;    
			minGapDist    = p.minGapDist;   
			maxGapDist    = p.maxGapDist;   
			distThreshold = p.distThreshold;
			minHashLen    = p.minHashLen;
			productSize   = p.productSize;
			lisWindow     = p.lisWindow;
			lisThreshold  = p.lisThreshold;
			mergeThreshold= p.mergeThreshold;
			ratioThreshold= p.ratioThreshold;
			doCheckpoint  = p.doCheckpoint;
			loadCheckpoint= p.loadCheckpoint;
			checkpointFile= p.checkpointFile;
			lisGapPenalty = p.lisGapPenalty;
			lisInsertPenalty = p.lisInsertPenalty;
			outputTuple    = p.outputTuple;
			outputGlue     = p.outputGlue;
			maskRepeats    = p.maskRepeats;
		}
    return *this;
  }
  friend std::istream & operator>> (std::istream &in, T_Params &p) {
    in.read((char*) & p.hashStep, (std::streamsize)(sizeof(ssize_t)));
    in.read((char*) & p.hashLength, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.wordLength, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.doIndel, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.neighborDist, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.alignStrips, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.alignGaps, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.minGapDist, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.maxGapDist, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.distThreshold, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.minHashLen, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.productSize, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.lisThreshold, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.lisWindow, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.mergeThreshold, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.ratioThreshold, (std::streamsize)(sizeof( double)));
    in.read((char*) & p.doCheckpoint, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.outputGlue, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.outputTuple, (std::streamsize)(sizeof( ssize_t)));
    in.read((char*) & p.maskRepeats, (std::streamsize)(sizeof( ssize_t)));
    return in.read((char*) & p.loadCheckpoint, (std::streamsize)(sizeof( ssize_t)));
  }

  friend std::ostream &operator<<(std::ostream &out, const T_Params &p) {

    out.write((char*) & p.hashStep, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.hashLength, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.wordLength, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.doIndel, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.neighborDist, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.alignStrips, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.alignGaps, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.minGapDist, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.maxGapDist, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.distThreshold, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.minHashLen, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.productSize, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.lisThreshold, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.lisWindow, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.mergeThreshold, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.ratioThreshold, (std::streamsize)(sizeof( double)));
    out.write((char*) & p.doCheckpoint, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.outputGlue, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.outputTuple, (std::streamsize)(sizeof( ssize_t)));
    out.write((char*) & p.maskRepeats, (std::streamsize)(sizeof( ssize_t)));
    return out.write((char*) & p.loadCheckpoint, (std::streamsize)(sizeof( ssize_t)));
  }
  
  void AlignStrips() {
    alignStrips = 1;
    alignGaps   = 0;
  }
  void AlignGaps() {
    alignStrips = 0;
    alignGaps  = 1;
  }
};


class T_ParallelParams {
public:
  ssize_t numFragments;
  ssize_t numChildren;
  ssize_t mytid;
  ssize_t *childIds;
  T_ParallelParams() {
    numFragments = 1; 
    numChildren  = 0;
    mytid = 0;
    childIds = NULL;
  }
};

#endif
