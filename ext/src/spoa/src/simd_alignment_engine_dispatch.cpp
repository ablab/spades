// Copyright (c) 2020 Mario Brcic, Robert Vaser

#include "simd_alignment_engine_implementation.hpp"

#if defined(__AVX2__)
  #define ARCH Architecture::kAVX2
#elif defined(__SSE4_1__)
  #define ARCH Architecture::kSSE4_1
#else
  #define ARCH Architecture::kSSE2
#endif

namespace spoa {

template class SimdAlignmentEngine<ARCH>;

}  // namespace spoa
