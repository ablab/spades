/***************************************************************************
 * Title:          NewettTree.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "NewettTree.h"
void Advance(std::ifstream &in) {
  char next = in.peek();
  while (in.good() && 
	 (next == ' ' or
	  next == '\n' or
	  next == '\t') and
	 next != ')' and 
	 next != '(') {
    in.get();
    next = in.peek();
  }
}


void ReadData(std::ifstream &in, std::string &data) {

  data = ""; // reset the input
  char next = in.peek();
  while (next != '(' &&
	 next != ')' &&
	 in.good()) {
    data.append(1, next);
    in.get();
    next = in.peek();
  }
}
