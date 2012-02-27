/***************************************************************************
 * Title:          FragmentUtils.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef FRAGMENT_UTILS_H_
#define FRAGMENT_UTILS_H_

#include <vector>
#include <string>
#include <set>
#include "lav/LAVBlock.h"
#include "DNASequence.h"
#include "SeqReader.h"

class Locus {
public:
  Locus() {
    seqid = "";
    start = -1;
    end   = -1;
  }
  DNASequence locusSeq;
  ssize_t start, end;
  std::string seqid;
};

typedef std::vector<Locus*> Loci;
typedef std::set<std::string> SequenceIds;

ssize_t FindSurroundingRegion(ssize_t start, ssize_t end,
			  std::vector<LAVBlock*> &blocks,
			  ssize_t &regionStart, ssize_t &regionEnd,
			  ssize_t seqLength);


void PrintTempSeq(DNASequence &sequence,
		  ssize_t sbjctStart, ssize_t sbjctEnd,
		  std::string &sbjctName);


ssize_t  ReadLoci(std::string &componentFileName, Loci &components);
ssize_t  ReadLociSequences(Loci &loci, std::string &seqDir);
void PrintLoci(Loci &loci, std::string &outName);

#endif
