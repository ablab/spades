/***************************************************************************
 * Title:          RMOutFileParser.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/18/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "RMOutFileParser.h"
#include <stdlib.h>
#include <string>

ssize_t RMOutFileParser::ParseRMFile(std::istream &in, RMEntries &repeats) {
  std::string temp;
  // discard the first 3 lines
  std::getline(in, temp);
  std::getline(in, temp);
  std::getline(in, temp);
  if (! in.good()) return 0;
  RMEntry *entry;
  entry = new RMEntry;
  while (ParseRMLine(in, *entry))
    repeats.push_back(entry);
  // created one two many
  delete entry;
  return 1;
}

ssize_t RMOutFileParser::ParseRMLine(std::istream &in, RMEntry &repeat) {
  ssize_t status;
  std::string posStr, leftStr;
  if (!in.good()) status = 0;
  if (!(in >> repeat.swScore) &&
      (in >> repeat.pctDiv) && 
      (in >> repeat.pctDel) &&
      (in >> repeat.pctIns) &&
      (in >> repeat.qrySeq) &&
      (in >> repeat.seqBegin) && 
      (in >> repeat.seqEnd) && 
      (in >> repeat.left) &&
      (in >> repeat.strand) &&
      (in >> repeat.repName) &&
      (in >> repeat.repClass) &&
      (in >> posStr) && 
      (in >> repeat.repEnd) &&
      (in >> leftStr) &&
      (in >> repeat.id))  return 0;
  
  // now parse the (potentially not just a number positions;
  if (posStr[0] == '(') {
    posStr.erase(posStr.begin());
    posStr.erase(posStr.begin() + posStr.size() - 1);
  }
  repeat.repStart = atoi(posStr.c_str());

  if (leftStr[0] == '(') {
    leftStr.erase(leftStr.begin());
    leftStr.erase(leftStr.begin() + leftStr.size() - 1);
  }
  repeat.repLeft = atoi(leftStr.c_str());

  return 1;
}

