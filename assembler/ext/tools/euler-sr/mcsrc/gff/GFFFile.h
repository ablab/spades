/***************************************************************************
 * Title:          GFFFile.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef GFFFILE_H_
#define GFFFILE_H_

#include <vector>
#include <istream>
#include <fstream>
#include <string>

#include "GFFEntry.h"
#include "lav/LAVFile.h"
#include "utils.h"

class GFFFile {
public:
  std::vector<GFFEntry*> entries;
  ssize_t ParseGFFFile(std::istream &in);
  ssize_t ParseGFFFile(std::string &gffFileName) {
    std::ifstream in;
    openck(gffFileName, in, std::ios::in);
    return ParseGFFFile(in);
  }
  void ToLAVFile(LAVFile &lavFile, 
		std::string &refName, ssize_t refLength, 
		std::string &qryName, ssize_t qryLength);
};

#endif
