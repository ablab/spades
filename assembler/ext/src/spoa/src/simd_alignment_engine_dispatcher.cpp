// Copyright (c) 2020 Mario Brcic, Robert Vaser

#include "simd_alignment_engine_implementation.hpp"

#ifdef SPOA_GENERATE_DISPATCH

#include "cpuinfo_x86.h"  // NOLINT

static const cpu_features::X86Features features =
    cpu_features::GetX86Info().features;

#endif

namespace spoa {

#ifndef SPOA_GENERATE_DISPATCH

template class SimdAlignmentEngine<Architecture::kAutomatic>;

#endif

std::unique_ptr<AlignmentEngine> CreateSimdAlignmentEngine(
    AlignmentType type,
    AlignmentSubtype subtype,
    std::int8_t m,
    std::int8_t n,
    std::int8_t g,
    std::int8_t e,
    std::int8_t q,
    std::int8_t c) {
#ifdef SPOA_GENERATE_DISPATCH
  if (features.avx2) {
    return SimdAlignmentEngine<Architecture::kAVX2>::Create(
        type, subtype, m, n, g, e, q, c);
  } else if (features.sse4_1) {
    return SimdAlignmentEngine<Architecture::kSSE4_1>::Create(
        type, subtype, m, n, g, e, q, c);
  } else {
    return SimdAlignmentEngine<Architecture::kSSE2>::Create(
        type, subtype, m, n, g, e, q, c);
  }
#else
  return SimdAlignmentEngine<Architecture::kAutomatic>::Create(
      type, subtype, m, n, g, e, q, c);
#endif
}

}  // namespace spoa
