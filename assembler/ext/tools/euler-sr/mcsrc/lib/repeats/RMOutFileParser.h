/***************************************************************************
 * Title:          RMOutFileParser.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef RM_OUT_FILE_PARSER
#define RM_OUT_FILE_PARSER

#include "RMEntry.h"
#include <fstream>
class RMOutFileParser {
public:
  static ssize_t ParseRMLine(std::istream &in, RMEntry &repeat);
  static ssize_t ParseRMFile(std::istream &in, RMEntries &repeats);
};

#endif
