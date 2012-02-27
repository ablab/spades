/***************************************************************************
 * Title:          LAVFile.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/18/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _LAVFILE_H
#define _LAVFILE_H

#include <string>
#include <vector>
#include <map>

#include "LAVBlock.h"
#include "LAVAlignedContig.h"
#include "mctypes.h"


class LAVFile {
public:
  static const char* standardBlastzOpts;
  std::string blastzOpts;
  std::vector <LAVAlignedContig*> alignments;
  ssize_t size() { return alignments.size(); }
  ~LAVFile() {
    ssize_t i;
    for (i =0; i < alignments.size(); i++) {
      delete alignments[i];
    }
  }
  void ParseBlastzOpts(std::string options,
											 IntMatrix &scoreMat,
											 std::map<char, ssize_t> &options2);
};

#endif
