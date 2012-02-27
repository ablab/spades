/***************************************************************************
 * Title:          RmvDuplicates.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <string>
#include "bbbwt/BBBWTQuery.h"
#include "DNASequence.h"
#include "SeqUtils.h"
#include "SeqReader.h"
#include "utils.h"

int main(int argc, char* argv[]){
  BBBWT csa;
  std::string seqFileName, csaFileName;
  std::string outFileName;
  ssize_t numRemoved=0;
  
  if (argc<3){
    std::cout <<"usage: rmvDuplicates infile outfile"<<std::endl;
    return 1;
  }
  seqFileName= argv[1];
  outFileName=argv[2];
  DNASequence seq;
  std::ifstream seqIn;
  openck(seqFileName,seqIn, std::ios::in);
  ssize_t low, high;
  std::ofstream out;
  openck(outFileName+".tmp",out, std::ios::out);
  //first pass check for duplicates
  while(SeqReader::GetSeq(seqIn, seq)){
    if (BW::Query(seq, csa, low, high) <=0){
      BW::Store(seq,csa);
      seq.PrintSeq(out);
      out << std::endl;
    }else{
      numRemoved++;
    }
  }
  out.close();
  seqIn.close();
  std::ifstream seqIn2;
  std::ofstream out2;
  //second pass check for redundant reads (substrings)
  openck(outFileName+".tmp",seqIn2, std::ios::in);
  openck(outFileName,out2,std::ios::out);
  while(SeqReader::GetSeq(seqIn2,seq)){
    if (BW::Query(seq,csa, low,high) <=1){
      seq.PrintSeq(out2);
      out2 << std::endl;
    }else{
      numRemoved++;
    }
  }
  std::string oldFile=outFileName+".tmp";
  std::remove(oldFile.c_str());
  out2.close();
  seqIn2.close();
  std::cout<<"Removed "<<numRemoved<<" redundant reads"<<std::endl;
  return 0;
}
