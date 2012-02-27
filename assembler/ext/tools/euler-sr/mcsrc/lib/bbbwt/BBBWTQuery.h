/***************************************************************************
 * Title:          BBBWTQuery.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef BBWT_QUERY_H_
#define BBWT_QUERY_H_

#include "compatibility.h"
#include "BBBWT.h"
#include "bitdictionary.h"
#include "BWdictionary.h"
#include "DNASequence.h"

typedef BW::bbbwt_T<BW::dictionary_T<bitdictionary_T<48> > >  BBBWT;

// Define accessor functions to the bbbwt


namespace BW {
  ssize_t Store(SimpleSequence &seq, BBBWT &csa);
  ssize_t Store(SimpleSequence &seq, BBBWT &csa, ssize_t start, ssize_t length);
  ssize_t Query(SimpleSequence &seq, BBBWT &csa, ssize_t &low, ssize_t &high);
  ssize_t Query(SimpleSequence &seq, ssize_t start, ssize_t length, BBBWT &csa, ssize_t &low, ssize_t &high);
  void Write(BBBWT &bbb, std::string &outFileName);
  void Read(std::string &inFileName, BBBWT &bbb);
}; // namespace

#endif
