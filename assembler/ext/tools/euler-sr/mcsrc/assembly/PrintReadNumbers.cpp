/***************************************************************************
 * Title:          PrintReadNumbers.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SimpleSequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include <fstream>
#include <sstream>
#include "utils.h"
#include "IntegralTupleStatic.h"

int main(int argc, char*argv[]) {
  if (argc != 3) {
    std::cout << "usage: printReadNumbers readsIn readsOut " << std::endl;
    return 1;
  }

  SimpleSequenceList reads;
  std::string readsFileName = argv[1];
  ReadSimpleSequences(readsFileName, reads);
  AppendReverseComplements(reads);
  std::ofstream out;
  openck(argv[2], out, std::ios::out);
  std::stringstream titleStrm;
  std::string title;

  ssize_t i;
  for (i = 0; i < reads.size(); i++ ) {
    title = "";
    titleStrm.str(title);
    titleStrm << i;
    title = titleStrm.str();
    reads[i].PrintSeq(out, title);
  }
  out.close();
  return 0;
}

