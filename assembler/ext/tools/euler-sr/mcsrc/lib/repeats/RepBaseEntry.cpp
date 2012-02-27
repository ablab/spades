/***************************************************************************
 * Title:          RepBaseEntry.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "RepBaseEntry.h"
#include "SeqReader.h"
#include <sstream>

ssize_t RepBaseEntry::ParseRepBaseEntry(std::ifstream &in) {
  ssize_t result;
  result = SeqReader::GetSeq(in, reference);
  if (! result ) return 0;

  // now parse the name and title from the reference sequence name
  std::stringstream strm;
  ssize_t space;
  for (space = 0; space < reference.namestr.size(); space++) 
    if (reference.namestr[space] == '\t') break;

  if (space < reference.namestr.size()) {
    name = std::string(reference.namestr.c_str(), 
		       reference.namestr.c_str() + space);
    repClass = std::string(reference.namestr.c_str() + space + 1,
			   reference.namestr.c_str() + reference.namestr.size());
  }
  return result;
}
