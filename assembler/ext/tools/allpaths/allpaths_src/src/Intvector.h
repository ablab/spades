///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


// This file defines the typedef "intvector", which stores 
// a vector of ints, and the typedef "vecintvector", which
// stores a vector of intvectors.

#ifndef INTVECTOR_H_
#define INTVECTOR_H_

#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"

typedef SerfVec<int> IntVec;
typedef MasterVec<IntVec> VecIntVec;

typedef SerfVec<unsigned int> UIntVec;
typedef MasterVec<UIntVec> VecUIntVec;

typedef SerfVec<unsigned short> UShortVec;
typedef MasterVec<UShortVec> VecUShortVec;

typedef SerfVec<long> LongVec;
typedef MasterVec<LongVec> VecLongVec;

typedef SerfVec<unsigned long> ULongVec;
typedef MasterVec<ULongVec> VecULongVec;

typedef SerfVec<uint64_t> UInt64Vec;
typedef MasterVec<UInt64Vec> UInt64VecVec;

typedef SerfVec<int64_t> Int64Vec;
typedef MasterVec<Int64Vec> Int64VecVec;


#endif // INTVECTOR_H_
