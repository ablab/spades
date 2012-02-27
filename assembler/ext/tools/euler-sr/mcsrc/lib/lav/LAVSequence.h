/***************************************************************************
 * Title:          LAVSequence.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _LAV_SEQUENCE
#define _LAV_SEQUENCE

#include <ostream>
#include <iostream>

class LAVSequence {
public:
  std::string sequenceName;
  ssize_t start, end, strand, contig;
  LAVSequence& operator=(const LAVSequence &s) {
		if (this != &s) {
			sequenceName= s.sequenceName;
			start = s.start;
			end   = s.end;
			strand = s.strand;
			contig = s.contig;
		}
    return *this;
  }
  std::ostream & PrintSequence(std::ostream &out) {
    
    if (sequenceName != "") 
      out << "  \"" << sequenceName << "\" ";
    else
      out << "\"\" ";
    //UNUSED// int tmpstrand;
    out << start << " " << end << " " << strand << " " << contig << std::endl;

		return out;
  }
};

#endif
