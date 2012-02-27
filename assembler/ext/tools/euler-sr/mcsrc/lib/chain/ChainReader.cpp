/***************************************************************************
 * Title:          ChainReader.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "ChainReader.h"

#include <iostream>
#include <fstream>

#include "utils.h"


ssize_t ChainReader::ReadChainHeader(std::ifstream &in, 
				 Chain &chain) {
  char tStrandChar, qStrandChar;
  std::string word;
  if (in.peek() == 'c') {
    in >> word; // discard 'chain'
    if (! (in >> chain.header.score >> chain.header.tName >> chain.header.tSize 
	   >> tStrandChar >> chain.header.tStart >> chain.header.tEnd 
	   >> chain.header.qName >> chain.header.qSize >> qStrandChar 
	   >> chain.header.qStart >> chain.header.qEnd >> chain.id)) 
      return 0;
    else {
      in.get(); // get the \n
      if (tStrandChar == '+')
	chain.header.tStrand = 0;
      else 
	chain.header.tStrand = 1;
      
      if (qStrandChar == '+')
	chain.header.qStrand = 0;
      else
	chain.header.qStrand = 1;
      return 1;
    }
  }
  return 0;
}

ssize_t ChainReader::ReadAlignment(std::ifstream &in,
			       Chain &chain){

  ssize_t size, dt, dq;
  size = 0;
  dt = 0;
  dq = 0;
  if (!in.good() || in.peek() == '\0')
    return 0;
  // If this call was started on a newline, get that and return.
  if (in.peek() == '\n') {
    in.get();
    return 0;
  }

  // Not a newline, read the size
  if (! (in >> size) ) {
    in.get();
    return -1;
  }

  // Successfully read in size.  If that was the only thing on the line
  // don't try to read dt and dq
  if (in.peek() == '\n') {
    in.get(); // get rid of newline
  }
  else {
    if (! (in >> dt >> dq)) {
      // Handle error reading in the 
      return -1;
    }
    in.get(); // read the newline after dt and dq
  }
  //  std::cout << size << " " << dt << " " << dq << std::endl;
  chain.size.push_back(size);
  chain.dt.push_back(dt);
  chain.dq.push_back(dq);
  return 1;
}

void ChainReader::ReadChainFile(std::string chainFileName, 
			       std::vector<Chain*> &chains) {

  std::ifstream chainIn;
  openck(chainFileName, chainIn);
  
  Chain chain, *newChain;

  while (ReadChainHeader(chainIn, chain)) {
    newChain = new Chain;
    *newChain = chain;
    while (ReadAlignment(chainIn, *newChain) > 0) ;
    chains.push_back(newChain);
  }
}
