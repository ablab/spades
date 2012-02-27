/***************************************************************************
 * Title:          LAVAlignedContig.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _LAV_ALIGNED_CONTIG_
#define _LAV_ALIGNED_CONTIG_

#include <vector>
#include <string>
#include <iostream>
#include <ostream>

#include "LAVSequence.h"
#include "LAVBlock.h"

class LAVAlignedContig {
public:
  LAVSequence refContig, qryContig;
  std::string refContigName, qryContigName;
  std::vector<LAVBlock*> alignments;
  ssize_t size() { return alignments.size(); }
  ssize_t push_back(LAVBlock* block) {
    alignments.push_back(block);
    return alignments.size();
  }
  std::ostream& PrintAlignment(std::ostream &out) {
    out << "s {" << std::endl;
    refContig.PrintSequence(out);
    qryContig.PrintSequence(out);
    out << "}" << std::endl;
    out << "h {" << std::endl;
    out << "   \"" << refContigName << "\"" << std::endl;  
    out << "   \"" << qryContigName << "\"" << std::endl;
    out << "}" << std::endl;
    ssize_t i;
    for (i = 0; i < alignments.size(); i++) {
      alignments[i]->PrintBlock(out);
    }
    return out;
  }
  ~LAVAlignedContig() {
    ssize_t i;
    for (i = 0; i < alignments.size(); i++) {
      delete alignments[i];
    }
  }
};

#endif
