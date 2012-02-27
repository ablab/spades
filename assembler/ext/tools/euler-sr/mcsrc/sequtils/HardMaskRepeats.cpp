/***************************************************************************
 * Title:          HardMaskRepeats.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include "SeqReader.h"
#include "SeqUtils.h"
#include "DNASequence.h"
#include "utils.h"

int main(int argc, char *argv[]) {

  if (argc != 3) {
    std::cout << "usage: mskrep infile outfile " << std::endl;
    std::cout << "       hard-masks the repeats (lower-case to N)" << std::endl;
    exit(1);
  }
  std::string inFileName, outFileName;

  inFileName = argv[1];
  outFileName= argv[2];

  std::ifstream in;
  std::ofstream out;
  
  openck(inFileName, in);
  openck(outFileName, out, std::ios::out);

  DNASequence seq, norep;

  SeqReader::MaskRepeats();

  while (SeqReader::GetSeq(in, seq, SeqReader::noConvert)) {
    if (seq.length == 0) 
      break;
    //UNUSED// ssize_t nonMasked = 0;
    //UNUSED+// unsigned char *p;
    unsigned char *c,  *end;
    end = &seq.seq[0] + seq.length;
    std::cout << "converting seq of len: " << seq.length <<std::endl;
    for (c = &seq.seq[0]; c != end; c++) (numeric_nuc_index[*c] >= 0 and numeric_nuc_index[*c] <= 3) or (*c='N');
    seq.PrintSeq(out);
    out << std::endl;
    seq.Reset();
  }
}
  





  
  
    
