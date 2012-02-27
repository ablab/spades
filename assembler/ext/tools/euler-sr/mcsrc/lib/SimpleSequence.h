/***************************************************************************
 * Title:          SimpleSequence.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef SIMPLE_SEQUENCE_H_
#define SIMPLE_SEQUENCE_H_

#include <vector>
#include <string>
#include <iostream>

#include "compatibility.h"

class SimpleSequence {
 public:
  unsigned char* seq;
  _SSZT_ length; // TODO: ideally this should be size_t (unsigned), but would need to fix all uses that combine it with signed numbers

  SimpleSequence() {
    seq = NULL;
    length = 0;
  }
  SimpleSequence & operator=(const SimpleSequence &s) {
		if (this != &s) {
			seq    = s.seq;
			length = s.length;
		}
    return *this;
  }
  std::ostream &PrintSeq(std::ostream &out, std::string title,
												 ssize_t lineLength = 50) {
    // simple sequences do not have titles, we have to supply then
    ssize_t i;
    out << ">" << title << std::endl;
    for (i = 0; i < length; i++ ) {
      if (i % lineLength == 0 and i > 0 and i < length-1)
				out << std::endl;
      out << seq[i];
    }
    out << std::endl;
		return out;
  }
};
typedef std::vector<SimpleSequence> SimpleSequenceList;
#endif
