///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef VECSTRING_H
#define VECSTRING_H

#include "String.h"
#include "feudal/MasterVec.h"
#include <cstddef>

typedef MasterVec<String> vecString;

inline
void Convert(const vec<String> &oldVec, vecString &newVec)
{ // probably wanted to clear newVec, but to maintain complete compatibility, we won't
  newVec.reserve(newVec.size()+oldVec.size());
  for (size_t idx = 0; idx < oldVec.size(); ++idx)
    newVec.push_back(oldVec[idx]); }

#endif
