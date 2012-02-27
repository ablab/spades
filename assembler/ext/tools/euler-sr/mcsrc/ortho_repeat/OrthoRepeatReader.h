/***************************************************************************
 * Title:          OrthoRepeatReader.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef ORTHO_REPEAT_READER_H_
#define ORTHO_REPEAT_READER_H_

#include <iostream>
#include <fstream>
#include "OrthoRepeat.h"

class OrthoRepeatReader {
public:
  std::ifstream &infile;
  OrthoRepeatReader(std::ifstream &in) :infile(in) { }
  ssize_t GetNextRepeat(OrthoRepeat *&repeat);
};


#endif
