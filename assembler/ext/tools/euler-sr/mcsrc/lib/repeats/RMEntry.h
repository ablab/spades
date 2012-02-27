/***************************************************************************
 * Title:          RMEntry.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef RM_ENTRY_H_
#define RM_ENTRY_H_

#include <string>
#include <vector>
// Sample line:
//4203  17.6  0.9  2.8  NT_004321   11055   11639 (1149409) C  L1MB3          LINE/L1                (0) 6183   5616      1

class RMEntry {
public:
  ssize_t swScore;
  double pctDiv;
  double pctDel;
  double pctIns;
  std::string qrySeq;
  ssize_t seqBegin;
  ssize_t seqEnd;
  ssize_t left;
  ssize_t strand;
  std::string repName;
  std::string repClass;
  std::string repStart;
  std::string repEnd;
  std::string repLeft;
  ssize_t id;
};
  
typedef std::vector<RMEntry*> RMEntries;

#endif
