/***************************************************************************
 * Title:          GapFunctions.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "GapFunctions.h"

ssize_t MuchGreater(ssize_t a, ssize_t b, double ratio) {
  if (b == 0) 
    return 1;
  if (double(a)/double(b) > ratio)
    return 1;
  else
    return 0;
}

ssize_t GetSpeciesNameFromFile(std::string fileName, std::string &speciesName) {

  ssize_t pos = fileName.find(".EN");
  if (pos != std::string::npos) {
    speciesName = fileName.substr(0,pos);
    return 1;
  }
  else {
    speciesName = "";
    return 0;
  }
}

ssize_t GetRegionNameFromFile(std::string fileName, std::string &regionName) {
  ssize_t pos = fileName.find("EN");
  if (pos != std::string::npos) {
    regionName = fileName.substr(pos,pos+6);
    return 1;
  }
  else {
    regionName = "";
    return 0;
  }
}
