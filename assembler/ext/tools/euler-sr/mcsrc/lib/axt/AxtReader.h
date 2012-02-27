/***************************************************************************
 * Title:          AxtReader.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _AXT_READER
#define _AXT_READER

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

#include "AxtEntry.h"

class AxtReader {
  // class to read axt files.  Store axt files in lav format.
  static ssize_t GetGapLocations(std::ifstream &axtIn, std::vector<ssize_t> &gapLocations);
  static ssize_t GetTitle(std::ifstream &axtIn, AxtEntry &axtEntry);

  static ssize_t GetAlignmentCoordinates(std::ifstream &axtIn);

  static ssize_t GetAxtEntry(std::ifstream &axtIn,
			 AxtEntry *&axtEntry);
  static ssize_t IsNuc(char c);
public:
  static void ReadAxtFile(std::string &axtFileName, 
			  AxtEntries &axtEntries);
};

#endif
