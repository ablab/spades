/***************************************************************************
 * Title:          AxtEntry.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _AXT_ENTRY
#define _AXT_ENTRY

#include <vector>
#include <string>
#include <ostream>

class AxtEntry {
public:
  ssize_t number;
  std::string refTitle;
  std::string qryTitle; 
  ssize_t refStart;
  ssize_t refEnd;
  ssize_t qryStart; 
  ssize_t qryEnd; 
  ssize_t qryStrand;
  ssize_t score;
  std::vector<ssize_t> refGapLocations;
  std::vector<ssize_t> qryGapLocations;
  void Print(std::ostream &out) {
    out << number << " " << refTitle << " " << refStart << " " << refEnd << " "
	<< qryTitle << " " << qryStart << " " << qryEnd << " " << qryStrand << " "
	<< score << " " << std::endl;
    ssize_t i;
    for (i = 0; i < refGapLocations.size(); i+=2) 
      out << refGapLocations[i] << " " << refGapLocations[i+1] << std::endl;
    for (i = 0; i < qryGapLocations.size(); i+=2) 
      out << qryGapLocations[i] << " " << qryGapLocations[i+1] << std::endl;

  }
};

typedef std::vector<AxtEntry*> AxtEntries;
#endif
