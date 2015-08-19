#ifndef __HAMMER_HKMER_HPP__
#define __HAMMER_HKMER_HPP__

#include "HSeq.hpp"

namespace hammer {

const uint32_t K = 16;
typedef HSeq<K> HKMer;

};

#endif // __HAMMER_HKMER_HPP__
