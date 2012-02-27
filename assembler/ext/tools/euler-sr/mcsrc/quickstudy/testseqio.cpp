/***************************************************************************
 * Title:          testseqio.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqReader.h"
#include "utils.h"

#include <ostream>
#include <fstream>

int main(int argc, char* argv[]) {

  std::string inName = argv[1];
  std::string outName = argv[2];

  DNASequence seq;
  SeqReader::GetSeq(inName, seq);
  std::ofstream out;
  openck(outName, out, std::ios::out);
  seq.PrintSeq(out);
  out << std::endl;
  return 1;
}

