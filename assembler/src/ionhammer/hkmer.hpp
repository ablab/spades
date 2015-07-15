//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_HKMER_HPP__
#define __HAMMER_HKMER_HPP__

#include "HSeq.hpp"

namespace hammer {

const uint32_t K = 16;
typedef HSeq<K> HKMer;

};

#endif // __HAMMER_HKMER_HPP__
