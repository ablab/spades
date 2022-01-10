// Copyright (c) 2020 Robert Vaser

#ifndef SIMD_ALIGNMENT_ENGINE_IMPLEMENTATION_HPP_
#define SIMD_ALIGNMENT_ENGINE_IMPLEMENTATION_HPP_

#include "simd_alignment_engine.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

extern "C" {
#ifdef SPOA_USE_SIMDE
  #if defined(__AVX2__)
    #include "simde/x86/avx2.h"
  #else
    #include "simde/x86/sse4.1.h"  // SSE4.1 is covered better
  #endif
#elif defined(__AVX2__) || defined(__SSE4_1__)
  #include <immintrin.h>  // AVX2 and lower
#endif
}

#include "spoa/graph.hpp"

namespace spoa {

// Taken from https://gcc.gnu.org/viewcvs/gcc?view=revision&revision=216149
inline void* align(
    std::size_t __align,
    std::size_t __size,
    void*& __ptr,  // NOLINT
    std::size_t& __space) noexcept {  // NOLINT
  const auto __intptr = reinterpret_cast<uintptr_t>(__ptr);
  const auto __aligned = (__intptr - 1u + __align) & -__align;
  const auto __diff = __aligned - __intptr;
  if ((__size + __diff) > __space) {
    return nullptr;
  } else {
    __space -= __diff;
    return __ptr = reinterpret_cast<void*>(__aligned);
  }
}

template<Architecture A, typename T>
T* AllocateAlignedMemory(
    T** storage,
    std::size_t size,
    std::size_t alignment) {
  *storage = new T[size + alignment - 1];
  void* ptr = static_cast<void*>(*storage);
  std::size_t storage_size = (size + alignment - 1) * sizeof(T);
  return static_cast<T*>(align(alignment, size * sizeof(T), ptr, storage_size));
}

template<Architecture A, typename T>
struct InstructionSet;

#if defined(__AVX2__)

constexpr std::uint32_t kRegisterSize = 256;
using __mxxxi = __m256i;

inline __mxxxi _mmxxx_load_si(__mxxxi const* mem_addr) {
  return _mm256_load_si256(mem_addr);
}

inline void _mmxxx_store_si(__mxxxi* mem_addr, const __mxxxi& a) {
  _mm256_store_si256(mem_addr, a);
}

inline __mxxxi _mmxxx_or_si(const __mxxxi& a, const __mxxxi& b) {
  return _mm256_or_si256(a, b);
}

#define _mmxxx_slli_si(a, n) n < 16 ? \
  _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 2, 0)), 16 - n) : \
  _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 2, 0))

#define _mmxxx_srli_si(a, n) \
  _mm256_srli_si256(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(2, 0, 0, 1)), n - 16)  // NOLINT

template<Architecture A>
struct InstructionSet<A, std::int16_t> {
  using type = std::int16_t;
  static constexpr std::uint32_t kNumVar = kRegisterSize / (8 * sizeof(type));
  static constexpr std::uint32_t kLogNumVar = 4;
  static constexpr std::uint32_t kLSS = 2;   // Left Shift Size
  static constexpr std::uint32_t kRSS = 30;  // Right Shift Size
  static inline __mxxxi _mmxxx_add_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm256_add_epi16(a, b);
  }
  static inline __mxxxi _mmxxx_sub_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm256_sub_epi16(a, b);
  }
  static inline __mxxxi _mmxxx_min_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm256_min_epi16(a, b);
  }
  static inline __mxxxi _mmxxx_max_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm256_max_epi16(a, b);
  }
  static inline __mxxxi _mmxxx_set1_epi(type a) {
    return _mm256_set1_epi16(a);
  }
  static inline void _mmxxx_prefix_max(
      __mxxxi& a,  // NOLINT
      const __mxxxi* masks,
      const __mxxxi* penalties) {
    a = _mmxxx_max_epi(a, _mmxxx_or_si(masks[0], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[0]), 2)));  // NOLINT
    a = _mmxxx_max_epi(a, _mmxxx_or_si(masks[1], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[1]), 4)));  // NOLINT
    a = _mmxxx_max_epi(a, _mmxxx_or_si(masks[2], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[2]), 8)));  // NOLINT
    a = _mmxxx_max_epi(a, _mmxxx_or_si(masks[3], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[3]), 16)));  // NOLINT
  }
};

template<Architecture A>
struct InstructionSet<A, std::int32_t> {
  using type = std::int32_t;
  static constexpr std::uint32_t kNumVar = kRegisterSize / (8 * sizeof(type));
  static constexpr std::uint32_t kLogNumVar = 3;
  static constexpr std::uint32_t kLSS = 4;
  static constexpr std::uint32_t kRSS = 28;
  static inline __mxxxi _mmxxx_add_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm256_add_epi32(a, b);
  }
  static inline __mxxxi _mmxxx_sub_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm256_sub_epi32(a, b);
  }
  static inline __mxxxi _mmxxx_min_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm256_min_epi32(a, b);
  }
  static inline __mxxxi _mmxxx_max_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm256_max_epi32(a, b);
  }
  static inline __mxxxi _mmxxx_set1_epi(type a) {
    return _mm256_set1_epi32(a);
  }
  static inline void _mmxxx_prefix_max(
      __mxxxi& a,  // NOLINT
      const __mxxxi* masks,
      const __mxxxi* penalties) {
    a = _mmxxx_max_epi(a, _mmxxx_or_si(masks[0], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[0]), 4)));  // NOLINT
    a = _mmxxx_max_epi(a, _mmxxx_or_si(masks[1], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[1]), 8)));  // NOLINT
    a = _mmxxx_max_epi(a, _mmxxx_or_si(masks[2], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[2]), 16)));  // NOLINT
    }
};

#elif defined(__SSE4_1__) || defined(SPOA_USE_SIMDE)

constexpr std::uint32_t kRegisterSize = 128;
using __mxxxi = __m128i;

inline __mxxxi _mmxxx_load_si(__mxxxi const* mem_addr) {
  return _mm_load_si128(mem_addr);
}

inline void _mmxxx_store_si(__mxxxi* mem_addr, const __mxxxi& a) {
  _mm_store_si128(mem_addr, a);
}

inline __mxxxi _mmxxx_or_si(const __mxxxi& a, const __mxxxi& b) {
  return _mm_or_si128(a, b);
}

#define _mmxxx_slli_si(a, n) \
    _mm_slli_si128(a, n)

#define _mmxxx_srli_si(a, n) \
    _mm_srli_si128(a, n)

template<Architecture A>
struct InstructionSet<A, std::int16_t> {
  using type = std::int16_t;
  static constexpr std::uint32_t kNumVar = kRegisterSize / (8 * sizeof(type));
  static constexpr std::uint32_t kLogNumVar = 3;
  static constexpr std::uint32_t kLSS = 2;
  static constexpr std::uint32_t kRSS = 14;
  static inline __mxxxi _mmxxx_add_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm_add_epi16(a, b);
  }
  static inline __mxxxi _mmxxx_sub_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm_sub_epi16(a, b);
  }
  static inline __mxxxi _mmxxx_min_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm_min_epi16(a, b);
  }
  static inline __mxxxi _mmxxx_max_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm_max_epi16(a, b);
  }
  static inline __mxxxi _mmxxx_set1_epi(type a) {
    return _mm_set1_epi16(a);
  }
  static inline void _mmxxx_prefix_max(
      __mxxxi& a,  // NOLINT
      const __mxxxi* masks,
      const __mxxxi* penalties) {
    a = _mmxxx_max_epi(a, _mmxxx_or_si(masks[0], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[0]), 2)));  // NOLINT
    a = _mmxxx_max_epi(a, _mmxxx_or_si(masks[1], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[1]), 4)));  // NOLINT
    a = _mmxxx_max_epi(a, _mmxxx_or_si(masks[2], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[2]), 8)));  // NOLINT
  }
};

template<Architecture A>
struct InstructionSet<A, std::int32_t> {
  using type = std::int32_t;
  static constexpr std::uint32_t kNumVar = kRegisterSize / (8 * sizeof(type));
  static constexpr std::uint32_t kLogNumVar = 2;
  static constexpr std::uint32_t kLSS = 4;
  static constexpr std::uint32_t kRSS = 12;
  static inline __mxxxi _mmxxx_add_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm_add_epi32(a, b);
  }
  static inline __mxxxi _mmxxx_sub_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm_sub_epi32(a, b);
  }
  static inline __mxxxi _mmxxx_min_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm_min_epi32(a, b);
  }
  static inline __mxxxi _mmxxx_max_epi(const __mxxxi& a, const __mxxxi& b) {
    return _mm_max_epi32(a, b);
  }
  static inline __mxxxi _mmxxx_set1_epi(type a) {
    return _mm_set1_epi32(a);
  }
  static inline void _mmxxx_prefix_max(
      __mxxxi& a,  // NOLINT
      const __mxxxi* masks,
      const __mxxxi* penalties) {
    a = _mmxxx_max_epi(a, _mmxxx_or_si(masks[0], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[0]), 4)));  // NOLINT
    a = _mmxxx_max_epi(a, _mmxxx_or_si(masks[1], _mmxxx_slli_si(_mmxxx_add_epi(a, penalties[1]), 8)));  // NOLINT
  }
};

#endif

#if defined(__AVX2__) || defined(__SSE4_1__) || defined(SPOA_USE_SIMDE)

template<Architecture A, typename T>
void _mmxxx_print(const __mxxxi& a) {
  __attribute__((aligned(kRegisterSize / 8))) typename T::type unpacked[T::kNumVar];  // NOLINT
  _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), a);

  for (std::uint32_t i = 0; i < T::kNumVar; i++) {
    std::cout << unpacked[i] << " ";
  }
}

template<Architecture A, typename T>
typename T::type _mmxxx_max_value(const __mxxxi& a) {
  typename T::type max_score = 0;
  __attribute__((aligned(kRegisterSize / 8))) typename T::type unpacked[T::kNumVar]; // NOLINT
  _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), a);

  for (std::uint32_t i = 0; i < T::kNumVar; i++) {
    max_score = std::max(max_score, unpacked[i]);
  }
  return max_score;
}

template<Architecture A, typename T>
typename T::type _mmxxx_value_at(const __mxxxi& a, std::uint32_t i) {
  __attribute__((aligned(kRegisterSize / 8))) typename T::type unpacked[T::kNumVar]; // NOLINT
  _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), a);

  return unpacked[i];
}

template<Architecture A, typename T>
std::int32_t _mmxxx_index_of(
    const __mxxxi* row,
    std::uint32_t row_width,
    typename T::type value) {
  for (std::uint32_t i = 0; i < row_width; ++i) {
    __attribute__((aligned(kRegisterSize / 8))) typename T::type unpacked[T::kNumVar]; // NOLINT
    _mmxxx_store_si(reinterpret_cast<__mxxxi*>(unpacked), row[i]);

    for (std::uint32_t j = 0; j < T::kNumVar; j++) {
      if (unpacked[j] == value) {
        return i * T::kNumVar + j;
      }
    }
  }
  return -1;
}

#endif

template<Architecture A>
std::unique_ptr<AlignmentEngine> SimdAlignmentEngine<A>::Create(
    AlignmentType type,
    AlignmentSubtype subtype,
    std::int8_t m,
    std::int8_t n,
    std::int8_t g,
    std::int8_t e,
    std::int8_t q,
    std::int8_t c) {
#if defined(__AVX2__) || defined(__SSE4_1__) || defined(SPOA_USE_SIMDE)
  return std::unique_ptr<AlignmentEngine>(
      new SimdAlignmentEngine<A>(type, subtype, m, n, g, e, q, c));
#else
  (void) type;
  (void) subtype;
  (void) m;
  (void) n;
  (void) g;
  (void) e;
  (void) q;
  (void) c;
  return nullptr;
#endif
}

template<Architecture A>
struct SimdAlignmentEngine<A>::Implementation {
#if defined(__AVX2__) || defined(__SSE4_1__) || defined(SPOA_USE_SIMDE)
  std::vector<std::uint32_t> node_id_to_rank;

  std::unique_ptr<__mxxxi[]> sequence_profile_storage;
  std::uint64_t sequence_profile_size;
  __mxxxi* sequence_profile;

  std::vector<std::int32_t> first_column;
  std::unique_ptr<__mxxxi[]> M_storage;
  std::uint64_t M_size;
  __mxxxi* H;
  __mxxxi* F;
  __mxxxi* E;
  __mxxxi* O;
  __mxxxi* Q;

  std::unique_ptr<__mxxxi[]> masks_storage;
  std::uint32_t masks_size;
  __mxxxi* masks;

  std::unique_ptr<__mxxxi[]> penalties_storage;
  std::uint32_t penalties_size;
  __mxxxi* penalties;

  Implementation()
      : node_id_to_rank(),
        sequence_profile_storage(nullptr),
        sequence_profile_size(0),
        sequence_profile(nullptr),
        first_column(),
        M_storage(nullptr),
        M_size(0),
        H(nullptr),
        F(nullptr),
        E(nullptr),
        O(nullptr),
        Q(nullptr),
        masks_storage(nullptr),
        masks_size(0),
        masks(nullptr),
        penalties_storage(nullptr),
        penalties_size(0),
        penalties(nullptr) {
  }
#endif
};

template<Architecture A>
SimdAlignmentEngine<A>::SimdAlignmentEngine(
    AlignmentType type,
    AlignmentSubtype subtype,
    std::int8_t m,
    std::int8_t n,
    std::int8_t g,
    std::int8_t e,
    std::int8_t q,
    std::int8_t c)
    : AlignmentEngine(type, subtype, m, n, g, e, q, c),
      pimpl_(new Implementation()) {
}

template<Architecture A>
void SimdAlignmentEngine<A>::Prealloc(
    std::uint32_t max_sequence_len,
    std::uint8_t alphabet_size) {
  if (max_sequence_len > std::numeric_limits<int32_t>::max()) {
    throw std::invalid_argument(
        "[spoa::SimdAlignmentEngine::Prealloc] error: too large sequence!");
  }

#if defined(__AVX2__) || defined(__SSE4_1__) || defined(SPOA_USE_SIMDE)

  std::int64_t worst_case_score = WorstCaseAlignmentScore(
      static_cast<std::int64_t>(max_sequence_len) + 8,
      static_cast<std::int64_t>(max_sequence_len) * alphabet_size);

  if (worst_case_score < std::numeric_limits<std::int32_t>::min() + 1024) {
    return;
  } else if (worst_case_score < std::numeric_limits<std::int16_t>::min() + 1024) {  // NOLINT
    try {
      Realloc(
          (max_sequence_len / InstructionSet<A, std::int32_t>::kNumVar) + 1,
          static_cast<std::uint64_t>(max_sequence_len) * alphabet_size,
          alphabet_size);
    } catch (std::bad_alloc& ba) {
      throw std::invalid_argument(
          "[spoa::SimdAlignmentEngine::Prealloc] error: insufficient memory!");
    }
  } else {
    try {
      Realloc(
          (max_sequence_len / InstructionSet<A, std::int16_t>::kNumVar) + 1,
          static_cast<std::uint64_t>(max_sequence_len) * alphabet_size,
          alphabet_size);
    } catch (std::bad_alloc& ba) {
      throw std::invalid_argument(
          "[spoa::SimdAlignmentEngine::Prealloc] error: insufficient memory!");
    }
  }

#endif
  (void) alphabet_size;
}

template<Architecture A>
void SimdAlignmentEngine<A>::Realloc(
    std::uint64_t matrix_width,
    std::uint64_t matrix_height,
    std::uint8_t num_codes) {
#if defined(__AVX2__) || defined(__SSE4_1__) || defined(SPOA_USE_SIMDE)
  if (pimpl_->node_id_to_rank.size() < matrix_height - 1) {
    pimpl_->node_id_to_rank.resize(matrix_height - 1, 0);
  }
  if (pimpl_->sequence_profile_size < num_codes * matrix_width) {
    __mxxxi* storage = nullptr;
    pimpl_->sequence_profile_size = num_codes * matrix_width;
    pimpl_->sequence_profile = AllocateAlignedMemory<A>(
        &storage,
        pimpl_->sequence_profile_size,
        kRegisterSize / 8);
    pimpl_->sequence_profile_storage.reset();
    pimpl_->sequence_profile_storage = std::unique_ptr<__mxxxi[]>(storage);
  }
  if (subtype_ == AlignmentSubtype::kLinear) {
    if (pimpl_->first_column.size() < matrix_height) {
      pimpl_->first_column.resize(matrix_height, 0);
    }
    if (pimpl_->M_size < matrix_height * matrix_width) {
      __mxxxi* storage = nullptr;
      pimpl_->M_size = matrix_height * matrix_width;
      pimpl_->H = AllocateAlignedMemory<A>(
        &storage,
        pimpl_->M_size,
        kRegisterSize / 8);
      pimpl_->M_storage.reset();
      pimpl_->M_storage = std::unique_ptr<__mxxxi[]>(storage);
    }
  } else if (subtype_ == AlignmentSubtype::kAffine) {
    if (pimpl_->first_column.size() < 2 * matrix_height) {
      pimpl_->first_column.resize(2 * matrix_height, 0);
    }
    if (pimpl_->M_size < 3 * matrix_height * matrix_width) {
      __mxxxi* storage = nullptr;
      pimpl_->M_size = 3 * matrix_height * matrix_width;
      pimpl_->H = AllocateAlignedMemory<A>(
        &storage,
        pimpl_->M_size,
        kRegisterSize / 8);
      pimpl_->F = pimpl_->H + matrix_height * matrix_width;
      pimpl_->E = pimpl_->F + matrix_height * matrix_width;
      pimpl_->M_storage.reset();
      pimpl_->M_storage = std::unique_ptr<__mxxxi[]>(storage);
    }
  } else if (subtype_ == AlignmentSubtype::kConvex) {
    if (pimpl_->first_column.size() < 3 * matrix_height) {
      pimpl_->first_column.resize(3 * matrix_height, 0);
    }
    if (pimpl_->M_size < 5 * matrix_height * matrix_width) {
      __mxxxi* storage = nullptr;
      pimpl_->M_size = 5 * matrix_height * matrix_width;
      pimpl_->H = AllocateAlignedMemory<A>(
          &storage,
          pimpl_->M_size,
          kRegisterSize / 8);
      pimpl_->F = pimpl_->H + matrix_height * matrix_width;
      pimpl_->E = pimpl_->F + matrix_height * matrix_width;
      pimpl_->O = pimpl_->E + matrix_height * matrix_width;
      pimpl_->Q = pimpl_->O + matrix_height * matrix_width;
      pimpl_->M_storage.reset();
      pimpl_->M_storage = std::unique_ptr<__mxxxi[]>(storage);
    }
  }
  if (pimpl_->masks_size < InstructionSet<A, std::int16_t>::kLogNumVar + 1) {
    __mxxxi* storage = nullptr;
    pimpl_->masks_size = InstructionSet<A, std::int16_t>::kLogNumVar + 1;
    pimpl_->masks = AllocateAlignedMemory<A>(
        &storage,
        pimpl_->masks_size,
        kRegisterSize / 8);
    pimpl_->masks_storage.reset();
    pimpl_->masks_storage = std::unique_ptr<__mxxxi[]>(storage);
  }
  if (pimpl_->penalties_size < 2 * InstructionSet<A, std::int16_t>::kLogNumVar) {  // NOLINT
    __mxxxi* storage = nullptr;
    pimpl_->penalties_size = 2 * InstructionSet<A, std::int16_t>::kLogNumVar;
    pimpl_->penalties = AllocateAlignedMemory<A>(
        &storage,
        pimpl_->penalties_size,
        kRegisterSize / 8);
    pimpl_->penalties_storage.reset();
    pimpl_->penalties_storage = std::unique_ptr<__mxxxi[]>(storage);
  }
#endif
  (void) matrix_width;
  (void) matrix_height;
  (void) num_codes;
}

template<Architecture A> template<typename T>
void SimdAlignmentEngine<A>::Initialize(
    const char* sequence,
    const Graph& graph,
    std::uint64_t normal_matrix_width,
    std::uint64_t matrix_width,
    std::uint64_t matrix_height) noexcept {
#if defined(__AVX2__) || defined(__SSE4_1__) || defined(SPOA_USE_SIMDE)
  std::int32_t padding_penatly = -1 * std::max(
      std::max(abs(m_), abs(n_)),
      std::max(abs(g_), abs(q_)));

  __attribute__((aligned(kRegisterSize / 8))) typename T::type unpacked[T::kNumVar] = {};  // NOLINT

  for (std::uint32_t i = 0; i < graph.num_codes(); ++i) {
    char c = graph.decoder(i);
    for (std::uint32_t j = 0; j < matrix_width; ++j) {
      for (std::uint32_t k = 0; k < T::kNumVar; ++k) {
        unpacked[k] = (j * T::kNumVar + k) < normal_matrix_width ?
            (c == sequence[j * T::kNumVar + k] ? m_ : n_) : padding_penatly;
      }
      pimpl_->sequence_profile[i * matrix_width + j] =
          _mmxxx_load_si(reinterpret_cast<const __mxxxi*>(unpacked));
    }
  }

  const auto& rank_to_node = graph.rank_to_node();
  for (std::uint32_t i = 0; i < rank_to_node.size(); ++i) {
    pimpl_->node_id_to_rank[rank_to_node[i]->id] = i;
  }

  typename T::type kNegativeInfinity =
      std::numeric_limits<typename T::type>::min() + 1024;

  __mxxxi negative_infinities = T::_mmxxx_set1_epi(kNegativeInfinity);
  __mxxxi zeroes = T::_mmxxx_set1_epi(0);

  // initialize secondary matrices
  switch (subtype_) {
    case AlignmentSubtype::kConvex:
      for (std::uint32_t j = 0; j < matrix_width; ++j) {
        pimpl_->O[j] = negative_infinities;
        pimpl_->Q[j] = T::_mmxxx_set1_epi(q_ + j * T::kNumVar * c_);

        __mxxxi c = T::_mmxxx_set1_epi(c_);
        for (std::uint32_t k = 1; k < T::kNumVar; ++k) {
          c = _mmxxx_slli_si(c, T::kLSS);
          pimpl_->Q[j] = T::_mmxxx_add_epi(pimpl_->Q[j], c);
        }
      }
      pimpl_->first_column[2 * matrix_height] = 0;
      for (std::uint32_t i = 1; i < matrix_height; ++i) {
        const auto& edges = rank_to_node[i - 1]->inedges;
        std::int32_t penalty = edges.empty() ? q_ - c_ : kNegativeInfinity;
        for (const auto& it : edges) {
          std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;
          penalty = std::max(penalty, pimpl_->first_column[2 * matrix_height + pred_i]);  // NOLINT
        }
        pimpl_->first_column[2 * matrix_height + i] = penalty + c_;
      }
      // fall through
    case AlignmentSubtype::kAffine:
      for (std::uint32_t j = 0; j < matrix_width; ++j) {
        pimpl_->F[j] = negative_infinities;
        pimpl_->E[j] = T::_mmxxx_set1_epi(g_ + j * T::kNumVar * e_);

        __mxxxi e = T::_mmxxx_set1_epi(e_);
        for (std::uint32_t k = 1; k < T::kNumVar; ++k) {
          e = _mmxxx_slli_si(e, T::kLSS);
          pimpl_->E[j] = T::_mmxxx_add_epi(pimpl_->E[j], e);
        }
      }
      pimpl_->first_column[matrix_height] = 0;
      for (std::uint32_t i = 1; i < matrix_height; ++i) {
        const auto& edges = rank_to_node[i - 1]->inedges;
        std::int32_t penalty = edges.empty() ? g_ - e_ : kNegativeInfinity;
        for (const auto& it : edges) {
          std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;
          penalty = std::max(penalty, pimpl_->first_column[matrix_height + pred_i]);  // NOLINT
        }
        pimpl_->first_column[matrix_height + i] = penalty + e_;
      }
      // fall through
    case AlignmentSubtype::kLinear:
      break;
    default:
      break;
  }

  // initialize primary matrix
  switch (type_) {
    case AlignmentType::kSW:
      for (std::uint32_t j = 0; j < matrix_width; ++j) {
        pimpl_->H[j] = zeroes;
      }
      for (std::uint32_t i = 0; i < matrix_height; ++i) {
        pimpl_->first_column[i] = 0;
      }
      break;
    case AlignmentType::kNW:
      switch (subtype_) {
        case AlignmentSubtype::kConvex:
          for (std::uint32_t i = 0; i < matrix_height; ++i) {
            pimpl_->first_column[i] = std::max(
                pimpl_->first_column[matrix_height + i],
                pimpl_->first_column[2 * matrix_height + i]);
          }
          for (std::uint32_t j = 0; j < matrix_width; ++j) {
            pimpl_->H[j] = T::_mmxxx_max_epi(pimpl_->E[j], pimpl_->Q[j]);
          }
          break;
        case AlignmentSubtype::kAffine:
          for (std::uint32_t i = 0; i < matrix_height; ++i) {
            pimpl_->first_column[i] = pimpl_->first_column[matrix_height + i];
          }
          for (std::uint32_t j = 0; j < matrix_width; ++j) {
            pimpl_->H[j] = pimpl_->E[j];
          }
          break;
        case AlignmentSubtype::kLinear:
          pimpl_->first_column[0] = 0;
          for (std::uint32_t i = 1; i < matrix_height; ++i) {
            const auto& edges = rank_to_node[i - 1]->inedges;
            std::int32_t penalty = edges.empty() ? 0 : kNegativeInfinity;
            for (const auto& it : edges) {
              std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;
              penalty = std::max(penalty, pimpl_->first_column[pred_i]);
            }
            pimpl_->first_column[i] = penalty + g_;
          }
          for (std::uint32_t j = 0; j < matrix_width; ++j) {
            pimpl_->H[j] = T::_mmxxx_set1_epi(g_ + j * T::kNumVar * g_);
            __mxxxi g = T::_mmxxx_set1_epi(g_);

            for (std::uint32_t k = 1; k < T::kNumVar; ++k) {
              g = _mmxxx_slli_si(g, T::kLSS);
              pimpl_->H[j] = T::_mmxxx_add_epi(pimpl_->H[j], g);
            }
          }
        default:
          break;
      }
      break;
    case AlignmentType::kOV:
      switch (subtype_) {
        case AlignmentSubtype::kConvex:
          for (std::uint32_t j = 0; j < matrix_width; ++j) {
            pimpl_->H[j] = T::_mmxxx_max_epi(pimpl_->E[j],
            pimpl_->Q[j]);
          }
          break;
        case AlignmentSubtype::kAffine:
          for (std::uint32_t j = 0; j < matrix_width; ++j) {
            pimpl_->H[j] = pimpl_->E[j];
          }
          break;
        case AlignmentSubtype::kLinear:
          for (std::uint32_t j = 0; j < matrix_width; ++j) {
            pimpl_->H[j] = T::_mmxxx_set1_epi(g_ + j * T::kNumVar * g_);
            __mxxxi g = T::_mmxxx_set1_epi(g_);

            for (std::uint32_t k = 1; k < T::kNumVar; ++k) {
              g = _mmxxx_slli_si(g, T::kLSS);
              pimpl_->H[j] = T::_mmxxx_add_epi(pimpl_->H[j], g);
            }
          }
          break;
        default:
          break;
      }
      for (std::uint32_t i = 0; i < matrix_height; ++i) {
        pimpl_->first_column[i] = 0;
      }
      break;
    default:
      break;
  }
#endif
  (void) sequence;
  (void) graph;
  (void) normal_matrix_width;
  (void) matrix_width;
  (void) matrix_height;
}

template<Architecture A>
Alignment SimdAlignmentEngine<A>::Align(
    const char* sequence, std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score) {
  if (sequence_len > std::numeric_limits<int32_t>::max()) {
    throw std::invalid_argument(
        "[spoa::SimdAlignmentEngine::Align] error: too large sequence!");
  }

  if (graph.nodes().empty() || sequence_len == 0) {
    return Alignment();
  }

#if defined(__AVX2__) || defined(__SSE4_1__) || defined(SPOA_USE_SIMDE)

  std::int64_t worst_case_score = WorstCaseAlignmentScore(
      sequence_len + 8,
      graph.nodes().size());

  if (worst_case_score < std::numeric_limits<std::int32_t>::min() + 1024) {
    throw std::invalid_argument(
        "[spoa::SimdAlignmentEngine::Align] error: possible overflow!");
  } else if (worst_case_score < std::numeric_limits<std::int16_t>::min() + 1024) {  // NOLINT
    try {
      Realloc(
          std::ceil(static_cast<double>(sequence_len) / InstructionSet<A, std::int32_t>::kNumVar),  // NOLINT
          graph.nodes().size() + 1,
          graph.num_codes());
    } catch (std::bad_alloc& ba) {
      throw std::invalid_argument(
          "[spoa::SimdAlignmentEngine::Align] error: insufficient memory!");
    }
    Initialize<InstructionSet<A, std::int32_t>>(
        sequence,
        graph,
        sequence_len,
        std::ceil(static_cast<double>(sequence_len) / InstructionSet<A, std::int32_t>::kNumVar),  // NOLINT
        graph.nodes().size() + 1);

    if (subtype_ == AlignmentSubtype::kLinear) {
      return Linear<InstructionSet<A, std::int32_t>>(sequence_len, graph, score);  // NOLINT
    } else if (subtype_ == AlignmentSubtype::kAffine) {
      return Affine<InstructionSet<A, std::int32_t>>(sequence_len, graph, score);  // NOLINT
    } else if (subtype_ == AlignmentSubtype::kConvex) {
      return Convex<InstructionSet<A, std::int32_t>>(sequence_len, graph, score);  // NOLINT
    }
  } else {
    try {
      Realloc(
          std::ceil(static_cast<double>(sequence_len) / InstructionSet<A, std::int16_t>::kNumVar),  // NOLINT
          graph.nodes().size() + 1,
          graph.num_codes());
    } catch (std::bad_alloc& ba) {
      throw std::invalid_argument(
          "[spoa::SimdAlignmentEngine::Align] error: insufficient memory!");
    }
    Initialize<InstructionSet<A, std::int16_t>>(
        sequence,
        graph,
        sequence_len,
        std::ceil(static_cast<double>(sequence_len) / InstructionSet<A, std::int16_t>::kNumVar),  // NOLINT
        graph.nodes().size() + 1);

    if (subtype_ == AlignmentSubtype::kLinear) {
      return Linear<InstructionSet<A, std::int16_t>>(sequence_len, graph, score);  // NOLINT
    } else if (subtype_ == AlignmentSubtype::kAffine) {
      return Affine<InstructionSet<A, std::int16_t>>(sequence_len, graph, score);  // NOLINT
    } else if (subtype_ == AlignmentSubtype::kConvex) {
      return Convex<InstructionSet<A, std::int16_t>>(sequence_len, graph, score);  // NOLINT
    }
  }

#endif
  (void) sequence;
  (void) score;
  return Alignment();
}

template<Architecture A> template <typename T>
Alignment SimdAlignmentEngine<A>::Linear(
    std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score) noexcept {
#if defined(__AVX2__) || defined(__SSE4_1__) || defined(SPOA_USE_SIMDE)
  std::uint64_t normal_matrix_width = sequence_len;
  std::uint64_t matrix_width =
      std::ceil(static_cast<double>(sequence_len) / T::kNumVar);
  const auto& rank_to_node = graph.rank_to_node();

  typename T::type kNegativeInfinity =
      std::numeric_limits<typename T::type>::min() + 1024;

  __attribute__((aligned(kRegisterSize / 8))) typename T::type unpacked[T::kNumVar] = {0};  // NOLINT

  for (std::uint32_t i = 0, j = 0; i < T::kNumVar && j < T::kLogNumVar; ++i) {
    unpacked[i] = kNegativeInfinity;
    if ((i & (i + 1)) == 0) {
      pimpl_->masks[j++] = _mmxxx_load_si(
        reinterpret_cast<const __mxxxi*>(unpacked));
    }
  }
  pimpl_->masks[T::kLogNumVar] = _mmxxx_slli_si(
      T::_mmxxx_set1_epi(kNegativeInfinity),
      T::kLSS);

  pimpl_->penalties[0] = T::_mmxxx_set1_epi(g_);
  for (std::uint32_t i = 1; i < T::kLogNumVar; ++i) {
    pimpl_->penalties[i] = T::_mmxxx_add_epi(
        pimpl_->penalties[i - 1],
        pimpl_->penalties[i - 1]);
  }

  typename T::type max_score = type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;  // NOLINT
  std::int32_t max_i = -1;
  std::int32_t max_j = -1;
  std::uint32_t last_column_id = (normal_matrix_width - 1) % T::kNumVar;
  __mxxxi zeroes = T::_mmxxx_set1_epi(0);
  __mxxxi g = T::_mmxxx_set1_epi(g_);

  // alignment
  for (const auto& it : rank_to_node) {
    __mxxxi* char_profile = &(pimpl_->sequence_profile[it->code * matrix_width]);  // NOLINT

    std::uint32_t i = pimpl_->node_id_to_rank[it->id] + 1;
    std::uint32_t pred_i = it->inedges.empty() ?
        0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

    __mxxxi* H_row = &(pimpl_->H[i * matrix_width]);
    __mxxxi* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

    __mxxxi x = _mmxxx_srli_si(
        T::_mmxxx_set1_epi(pimpl_->first_column[pred_i]),
        T::kRSS);

    for (std::uint64_t j = 0; j < matrix_width; ++j) {
      // get diagonal
      __mxxxi t1 = _mmxxx_srli_si(H_pred_row[j], T::kRSS);
      H_row[j] = _mmxxx_or_si(
          _mmxxx_slli_si(H_pred_row[j], T::kLSS),
          x);
      x = t1;

      // update M
      H_row[j] = T::_mmxxx_max_epi(
          T::_mmxxx_add_epi(H_row[j], char_profile[j]),
          T::_mmxxx_add_epi(H_pred_row[j], g));
    }
    // check other predecessors
    for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
      pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

      H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

      x = _mmxxx_srli_si(
          T::_mmxxx_set1_epi(pimpl_->first_column[pred_i]),
          T::kRSS);

      for (std::uint64_t j = 0; j < matrix_width; ++j) {
        // get diagonal
        __mxxxi t1 = _mmxxx_srli_si(H_pred_row[j], T::kRSS);
        __mxxxi m = _mmxxx_or_si(
            _mmxxx_slli_si(H_pred_row[j], T::kLSS),
            x);
        x = t1;

        // updage M
        H_row[j] = T::_mmxxx_max_epi(
            H_row[j],
            T::_mmxxx_max_epi(
                T::_mmxxx_add_epi(m, char_profile[j]),
                T::_mmxxx_add_epi(H_pred_row[j], g)));
      }
    }

    __mxxxi score = T::_mmxxx_set1_epi(kNegativeInfinity);
    x = _mmxxx_srli_si(
        T::_mmxxx_add_epi(
            T::_mmxxx_set1_epi(pimpl_->first_column[i]),
            g),
        T::kRSS);

    for (std::uint64_t j = 0; j < matrix_width; ++j) {
      // add last element of previous vector into this one
      H_row[j] = T::_mmxxx_max_epi(
          H_row[j],
          _mmxxx_or_si(x, pimpl_->masks[T::kLogNumVar]));

      T::_mmxxx_prefix_max(H_row[j], pimpl_->masks, pimpl_->penalties);

      x = _mmxxx_srli_si(
          T::_mmxxx_add_epi(H_row[j], g),
          T::kRSS);

      if (type_ == AlignmentType::kSW) {
        H_row[j] = T::_mmxxx_max_epi(H_row[j], zeroes);
      }
      score = T::_mmxxx_max_epi(score, H_row[j]);
    }

    if (type_ == AlignmentType::kSW) {
      std::int32_t max_row_score = _mmxxx_max_value<A, T>(score);
      if (max_score < max_row_score) {
        max_score = max_row_score;
        max_i = i;
      }
    } else if (type_ == AlignmentType::kOV) {
      if (it->outedges.empty()) {
        std::int32_t max_row_score = _mmxxx_max_value<A, T>(score);
        if (max_score < max_row_score) {
          max_score = max_row_score;
          max_i = i;
        }
      }
    } else if (type_ == AlignmentType::kNW) {
      if (it->outedges.empty()) {
        std::int32_t max_row_score = _mmxxx_value_at<A, T>(
            H_row[matrix_width - 1],
            last_column_id);
        if (max_score < max_row_score) {
          max_score = max_row_score;
          max_i = i;
        }
      }
    }
  }

  if (max_i == -1 && max_j == -1) {
    return Alignment();
  }
  if (score) {
    *score = max_score;
  }

  if (type_ == AlignmentType::kSW) {
    max_j = _mmxxx_index_of<A, T>(
        &(pimpl_->H[max_i * matrix_width]),
        matrix_width,
        max_score);
  } else if (type_ == AlignmentType::kOV) {
    if (rank_to_node[max_i - 1]->outedges.empty()) {
      max_j = _mmxxx_index_of<A, T>(
          &(pimpl_->H[max_i * matrix_width]),
          matrix_width,
          max_score);
    } else {
      max_j = normal_matrix_width - 1;
    }
  } else if (type_ == AlignmentType::kNW) {
    max_j = normal_matrix_width - 1;
  }

  // backtrack
  std::uint32_t max_num_predecessors = 1;
  for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(max_i); ++i) {
    max_num_predecessors = std::max(
        max_num_predecessors,
        static_cast<std::uint32_t>(rank_to_node[i]->inedges.size()));
  }

  typename T::type* backtrack_storage = nullptr;
  typename T::type* H = AllocateAlignedMemory<A>(
      &backtrack_storage,
      3 * T::kNumVar + 2 * T::kNumVar * max_num_predecessors,
      kRegisterSize / 8);
  typename T::type* H_pred = H + T::kNumVar;
  typename T::type* H_diag_pred = H_pred + T::kNumVar * max_num_predecessors;
  typename T::type* H_left_pred = H_diag_pred + T::kNumVar * max_num_predecessors;  // NOLINT
  typename T::type* profile = H_left_pred + T::kNumVar;

  std::vector<std::uint32_t> predecessors;

  std::int32_t i = max_i;
  std::int32_t j = max_j;
  std::int32_t prev_i = 0, prev_j = 0;

  std::uint32_t j_div = j / T::kNumVar;
  std::uint32_t j_mod = j % T::kNumVar;

  bool load_next_segment = true;

  Alignment alignment;

  do {
    // check stop condition
    if (j == -1 || i == 0) {
      break;
    }

    const auto& it = rank_to_node[i - 1];
    // load everything
    if (load_next_segment) {
      predecessors.clear();

      // load current cells
      _mmxxx_store_si(
          reinterpret_cast<__mxxxi*>(H),
          pimpl_->H[i * matrix_width + j_div]);

      // load predecessors cells
      if (it->inedges.empty()) {
        predecessors.emplace_back(0);
        _mmxxx_store_si(reinterpret_cast<__mxxxi*>(H_pred), pimpl_->H[j_div]);
      } else {
        std::uint32_t store_pos = 0;
        for (const auto& jt : it->inedges) {
          predecessors.emplace_back(pimpl_->node_id_to_rank[jt->tail->id] + 1);
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(&H_pred[store_pos * T::kNumVar]),
              pimpl_->H[predecessors.back() * matrix_width + j_div]);
          ++store_pos;
        }
      }

      // load query profile cells
      _mmxxx_store_si(
          reinterpret_cast<__mxxxi*>(profile),
          pimpl_->sequence_profile[it->code * matrix_width + j_div]);
    }

    // check stop condition
    if (type_ == AlignmentType::kSW && H[j_mod] == 0) {
      break;
    }

    if (j_mod == 0) {
      // border case
      if (j_div > 0) {
        _mmxxx_store_si(
            reinterpret_cast<__mxxxi*>(H_left_pred),
            pimpl_->H[i * matrix_width + j_div - 1]);

        for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(&H_diag_pred[p * T::kNumVar]),
              pimpl_->H[predecessors[p] * matrix_width + (j_div - 1)]);
        }
      } else {
        H_left_pred[T::kNumVar - 1] = pimpl_->first_column[i];

        for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
          H_diag_pred[(p + 1) * T::kNumVar - 1] =
              pimpl_->first_column[predecessors[p]];
        }
      }
    }

    // find best predecessor cell
    bool predecessor_found = false;

    if (i != 0) {
      for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
        if ((j_mod == 0 && H[j_mod] == H_diag_pred[(p + 1) * T::kNumVar - 1] + profile[j_mod]) ||  // NOLINT
            (j_mod != 0 && H[j_mod] == H_pred[p * T::kNumVar + j_mod - 1] + profile[j_mod])) {  // NOLINT
          prev_i = predecessors[p];
          prev_j = j - 1;
          predecessor_found = true;
          break;
        }
      }
    }

    if (!predecessor_found && i != 0) {
      for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
        if (H[j_mod] == H_pred[p * T::kNumVar + j_mod] + g_) {
          prev_i = predecessors[p];
          prev_j = j;
          predecessor_found = true;
          break;
        }
      }
    }

    if (!predecessor_found) {
      if ((j_mod == 0 && H[j_mod] == H_left_pred[T::kNumVar - 1] + g_) ||
          (j_mod != 0 && H[j_mod] == H[j_mod - 1] + g_)) {
        prev_i = i;
        prev_j = j - 1;
        predecessor_found = true;
      }
    }

    alignment.emplace_back(
        i == prev_i ? -1 : rank_to_node[i - 1]->id,
        j == prev_j ? -1 : j);

    // update for next round
    load_next_segment =
        (i == prev_i ? false : true) ||
        (j != prev_j && prev_j % T::kNumVar == T::kNumVar - 1 ? true : false);

    i = prev_i;
    j = prev_j;
    j_div = j / T::kNumVar;
    j_mod = j % T::kNumVar;
  } while (true);

  delete[] backtrack_storage;

  // update alignment for NW (backtrack stops on first row or column)
  if (type_ == AlignmentType::kNW) {
    while (i == 0 && j != -1) {
      alignment.emplace_back(-1, j);
      --j;
    }
    while (i != 0 && j == -1) {
      alignment.emplace_back(rank_to_node[i - 1]->id, -1);

      const auto& it = rank_to_node[i - 1];
      if (it->inedges.empty()) {
          i = 0;
      } else {
        for (const auto& jt : it->inedges) {
          std::uint32_t pred_i = pimpl_->node_id_to_rank[jt->tail->id] + 1;
          if (pimpl_->first_column[i] == pimpl_->first_column[pred_i] + g_) {
            i = pred_i;
            break;
          }
        }
      }
    }
  }

  std::reverse(alignment.begin(), alignment.end());
  return alignment;
#else
  (void) sequence_len;
  (void) graph;
  (void) score;
  return Alignment();
#endif
}

template<Architecture A> template <typename T>
Alignment SimdAlignmentEngine<A>::Affine(
    std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score) noexcept {
#if defined(__AVX2__) || defined(__SSE4_1__) || defined(SPOA_USE_SIMDE)
  std::uint64_t normal_matrix_width = sequence_len;
  std::uint64_t matrix_width =
      std::ceil(static_cast<double>(sequence_len) / T::kNumVar);
  const auto& rank_to_node = graph.rank_to_node();

  typename T::type kNegativeInfinity =
      std::numeric_limits<typename T::type>::min() + 1024;

  typename T::type max_score = type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;  // NOLINT
  std::int32_t max_i = -1;
  std::int32_t max_j = -1;
  std::uint32_t last_column_id = (normal_matrix_width - 1) % T::kNumVar;
  __mxxxi zeroes = T::_mmxxx_set1_epi(0);
  __mxxxi g = T::_mmxxx_set1_epi(g_ - e_);
  __mxxxi e = T::_mmxxx_set1_epi(e_);

  __attribute__((aligned(kRegisterSize / 8))) typename T::type unpacked[T::kNumVar] = {0};  // NOLINT

  for (std::uint32_t i = 0, j = 0; i < T::kNumVar && j < T::kLogNumVar; ++i) {
    unpacked[i] = kNegativeInfinity;
    if ((i & (i + 1)) == 0) {
      pimpl_->masks[j++] = _mmxxx_load_si(
          reinterpret_cast<const __mxxxi*>(unpacked));
    }
  }
  pimpl_->masks[T::kLogNumVar] = _mmxxx_slli_si(
      T::_mmxxx_set1_epi(kNegativeInfinity),
      T::kLSS);

  pimpl_->penalties[0] = T::_mmxxx_set1_epi(e_);
  for (std::uint32_t i = 1; i < T::kLogNumVar; ++i) {
    pimpl_->penalties[i] = T::_mmxxx_add_epi(
        pimpl_->penalties[i - 1],
        pimpl_->penalties[i - 1]);
  }

  // alignment
  for (const auto& it : rank_to_node) {
    __mxxxi* char_profile = &(pimpl_->sequence_profile[it->code * matrix_width]);  // NOLINT

    std::uint32_t i = pimpl_->node_id_to_rank[it->id] + 1;

    __mxxxi* H_row = &(pimpl_->H[i * matrix_width]);
    __mxxxi* F_row = &(pimpl_->F[i * matrix_width]);

    std::uint32_t pred_i = it->inedges.empty() ?
        0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

    __mxxxi* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);
    __mxxxi* F_pred_row = &(pimpl_->F[pred_i * matrix_width]);

    __mxxxi x = _mmxxx_srli_si(
        T::_mmxxx_set1_epi(pimpl_->first_column[pred_i]),
        T::kRSS);

    for (std::uint64_t j = 0; j < matrix_width; ++j) {
      // update F
      F_row[j] = T::_mmxxx_add_epi(
          T::_mmxxx_max_epi(
              T::_mmxxx_add_epi(H_pred_row[j], g),
              F_pred_row[j]),
          e);

      // update H
      H_row[j] = T::_mmxxx_add_epi(
          _mmxxx_or_si(
              _mmxxx_slli_si(H_pred_row[j], T::kLSS),
              x),
          char_profile[j]);
      x = _mmxxx_srli_si(H_pred_row[j], T::kRSS);
    }
    // check other predecessors
    for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
      pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

      H_pred_row = &(pimpl_->H[pred_i * matrix_width]);
      F_pred_row = &(pimpl_->F[pred_i * matrix_width]);

      x = _mmxxx_srli_si(
          T::_mmxxx_set1_epi(pimpl_->first_column[pred_i]),
          T::kRSS);

      for (std::uint64_t j = 0; j < matrix_width; ++j) {
        // update F
        F_row[j] = T::_mmxxx_max_epi(
            F_row[j],
            T::_mmxxx_add_epi(
                T::_mmxxx_max_epi(
                    T::_mmxxx_add_epi(H_pred_row[j], g),
                    F_pred_row[j]),
                e));

        // update H
        H_row[j] = T::_mmxxx_max_epi(
            H_row[j],
            T::_mmxxx_add_epi(
                _mmxxx_or_si(
                    _mmxxx_slli_si(H_pred_row[j], T::kLSS),
                    x),
                char_profile[j]));
        x = _mmxxx_srli_si(H_pred_row[j], T::kRSS);
      }
    }

    __mxxxi* E_row = &(pimpl_->E[i * matrix_width]);
    __mxxxi score = zeroes;
    x = T::_mmxxx_set1_epi(pimpl_->first_column[i]);

    for (std::uint64_t j = 0; j < matrix_width; ++j) {
      H_row[j] = T::_mmxxx_max_epi(H_row[j], F_row[j]);

      E_row[j] = T::_mmxxx_add_epi(
          T::_mmxxx_add_epi(
              _mmxxx_or_si(
                  _mmxxx_slli_si(H_row[j], T::kLSS),
                  _mmxxx_srli_si(x, T::kRSS)),
              g),
          e);

      T::_mmxxx_prefix_max(E_row[j], pimpl_->masks, pimpl_->penalties);

      H_row[j] = T::_mmxxx_max_epi(H_row[j], E_row[j]);
      x = T::_mmxxx_max_epi(
          H_row[j],
          T::_mmxxx_sub_epi(E_row[j], g));

      if (type_ == AlignmentType::kSW) {
        H_row[j] = T::_mmxxx_max_epi(H_row[j], zeroes);
      }
      score = T::_mmxxx_max_epi(score, H_row[j]);
    }

    if (type_ == AlignmentType::kSW) {
      std::int32_t max_row_score = _mmxxx_max_value<A, T>(score);
      if (max_score < max_row_score) {
        max_score = max_row_score;
        max_i = i;
      }
    } else if (type_ == AlignmentType::kOV) {
      if (it->outedges.empty()) {
        std::int32_t max_row_score = _mmxxx_max_value<A, T>(score);
        if (max_score < max_row_score) {
          max_score = max_row_score;
          max_i = i;
        }
      }
    } else if (type_ == AlignmentType::kNW) {
      if (it->outedges.empty()) {
        std::int32_t max_row_score = _mmxxx_value_at<A, T>(
            H_row[matrix_width - 1],
            last_column_id);
        if (max_score < max_row_score) {
          max_score = max_row_score;
          max_i = i;
        }
      }
    }
  }

  if (max_i == -1 && max_j == -1) {
    return Alignment();
  }
  if (score) {
    *score = max_score;
  }

  if (type_ == AlignmentType::kSW) {
    max_j = _mmxxx_index_of<A, T>(
        &(pimpl_->H[max_i * matrix_width]),
        matrix_width, max_score);
  } else if (type_ == AlignmentType::kOV) {
    if (rank_to_node[max_i - 1]->outedges.empty()) {
      max_j = _mmxxx_index_of<A, T>(
          &(pimpl_->H[max_i * matrix_width]),
          matrix_width, max_score);
    } else {
      max_j = normal_matrix_width - 1;
    }
  } else if (type_ == AlignmentType::kNW) {
    max_j = normal_matrix_width - 1;
  }

  // backtrack
  std::uint32_t max_num_predecessors = 1;
  for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(max_i); ++i) {
    max_num_predecessors = std::max(
        max_num_predecessors,
        static_cast<std::uint32_t>(rank_to_node[i]->inedges.size()));
  }

  typename T::type* backtrack_storage = nullptr;
  typename T::type* H = AllocateAlignedMemory<A>(
      &backtrack_storage,
      6 * T::kNumVar + 3 * T::kNumVar * max_num_predecessors,
      kRegisterSize / 8);
  typename T::type* H_pred = H + T::kNumVar;
  typename T::type* H_diag_pred = H_pred + T::kNumVar * max_num_predecessors;
  typename T::type* H_left = H_diag_pred + T::kNumVar * max_num_predecessors;
  typename T::type* F = H_left + T::kNumVar;
  typename T::type* F_pred = F + T::kNumVar;
  typename T::type* E = F_pred + T::kNumVar * max_num_predecessors;
  typename T::type* E_left = E + T::kNumVar;
  typename T::type* profile = E_left + T::kNumVar;

  std::vector<std::uint32_t> predecessors;

  std::int32_t i = max_i;
  std::int32_t j = max_j;
  std::int32_t prev_i = 0, prev_j = 0;

  std::uint32_t j_div = j / T::kNumVar;
  std::uint32_t j_mod = j % T::kNumVar;

  bool load_next_segment = true;

  Alignment alignment;

  do {
    // check stop condition
    if (j == -1 || i == 0) {
      break;
    }

    const auto& it = rank_to_node[i - 1];
    // load everything
    if (load_next_segment) {
      predecessors.clear();

      // load current cells
      _mmxxx_store_si(
          reinterpret_cast<__mxxxi*>(H),
          pimpl_->H[i * matrix_width + j_div]);
      _mmxxx_store_si(
          reinterpret_cast<__mxxxi*>(E),
          pimpl_->E[i * matrix_width + j_div]);

      // load predecessors cells
      if (it->inedges.empty()) {
        predecessors.emplace_back(0);
        _mmxxx_store_si(reinterpret_cast<__mxxxi*>(H_pred), pimpl_->H[j_div]);
        _mmxxx_store_si(reinterpret_cast<__mxxxi*>(F_pred), pimpl_->F[j_div]);
      } else {
        std::uint32_t store_pos = 0;
        for (const auto& jt : it->inedges) {
          predecessors.emplace_back(pimpl_->node_id_to_rank[jt->tail->id] + 1);
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(&H_pred[store_pos * T::kNumVar]),
              pimpl_->H[predecessors.back() * matrix_width + j_div]);
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(&F_pred[store_pos * T::kNumVar]),
              pimpl_->F[predecessors.back() * matrix_width + j_div]);
          ++store_pos;
        }
      }

      // load query profile cells
      _mmxxx_store_si(
          reinterpret_cast<__mxxxi*>(profile),
          pimpl_->sequence_profile[it->code * matrix_width + j_div]);
    }

    // check stop condition
    if (type_ == AlignmentType::kSW && H[j_mod] == 0) {
      break;
    }

    if (j_mod == 0) {
      // border case
      if (j_div > 0) {
        for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(&H_diag_pred[p * T::kNumVar]),
              pimpl_->H[predecessors[p] * matrix_width + (j_div - 1)]);
        }
        _mmxxx_store_si(
            reinterpret_cast<__mxxxi*>(H_left),
            pimpl_->H[i * matrix_width + j_div - 1]);
        _mmxxx_store_si(
            reinterpret_cast<__mxxxi*>(E_left),
            pimpl_->E[i * matrix_width + j_div - 1]);
      } else {
        for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
          H_diag_pred[(p + 1) * T::kNumVar - 1] =
              pimpl_->first_column[predecessors[p]];
        }
        H_left[T::kNumVar - 1] = pimpl_->first_column[i];
        E_left[T::kNumVar - 1] = pimpl_->first_column[i];
      }
    }

    // find best predecessor cell
    bool predecessor_found = false, extend_left = false, extend_up = false;

    if (i != 0) {
      for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
        if ((j_mod == 0 && H[j_mod] == H_diag_pred[(p + 1) * T::kNumVar - 1] + profile[j_mod]) ||  // NOLINT
            (j_mod != 0 && H[j_mod] == H_pred[p * T::kNumVar + j_mod - 1] + profile[j_mod])) {  // NOLINT
          prev_i = predecessors[p];
          prev_j = j - 1;
          predecessor_found = true;
          break;
        }
      }
    }

    if (!predecessor_found && i != 0) {
      for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
        if ((extend_up = H[j_mod] == F_pred[p * T::kNumVar + j_mod] + e_) ||
                         H[j_mod] == H_pred[p * T::kNumVar + j_mod] + g_) {
          prev_i = predecessors[p];
          prev_j = j;
          predecessor_found = true;
          break;
        }
      }
    }

    if (!predecessor_found) {
      if ((j_mod != 0 && ((extend_left = H[j_mod] == E[j_mod - 1] + e_) ||
                                         H[j_mod] == H[j_mod - 1] + g_)) ||
          (j_mod == 0 && ((extend_left = H[j_mod] == E_left[T::kNumVar - 1] + e_ ) ||  // NOLINT
                                         H[j_mod] == H_left[T::kNumVar - 1] + g_))) {  // NOLINT
        prev_i = i;
        prev_j = j - 1;
        predecessor_found = true;
      }
    }

    alignment.emplace_back(
        i == prev_i ? -1 : rank_to_node[i - 1]->id,
        j == prev_j ? -1 : j);

    // update for next round
    load_next_segment =
        (i == prev_i ? false : true) ||
        (j != prev_j && prev_j % T::kNumVar == T::kNumVar - 1 ? true : false);

    i = prev_i;
    j = prev_j;
    j_div = j / T::kNumVar;
    j_mod = j % T::kNumVar;

    if (extend_left) {
      while (true) {
        // load
        if (j_mod == T::kNumVar - 1) {
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(E),
              pimpl_->E[i * matrix_width + j_div]);
        } else if (j_mod == 0) {  // boarder case
          if (j_div > 0) {
            _mmxxx_store_si(
                reinterpret_cast<__mxxxi*>(E_left),
                pimpl_->E[i * matrix_width + j_div - 1]);
          }
        }

        alignment.emplace_back(-1, j);
        --j;
        j_div = j / T::kNumVar;
        j_mod = j % T::kNumVar;
        if ((j == -1) ||
            (j_mod != T::kNumVar - 1 && E[j_mod] + e_ != E[j_mod + 1]) ||
            (j_mod == T::kNumVar - 1 && E_left[j_mod] + e_ != E[0])) {
          break;
        }
      }
      load_next_segment = true;
    } else if (extend_up) {
      while (true) {
        // load
        _mmxxx_store_si(
            reinterpret_cast<__mxxxi*>(F),
            pimpl_->F[i * matrix_width + j_div]);

        prev_i = 0;
        predecessors.clear();
        std::uint32_t store_pos = 0;
        for (const auto& it : rank_to_node[i - 1]->inedges) {
          predecessors.emplace_back(pimpl_->node_id_to_rank[it->tail->id] + 1);
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(&H_pred[store_pos * T::kNumVar]),
              pimpl_->H[predecessors.back() * matrix_width + j_div]);
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(&F_pred[store_pos * T::kNumVar]),
              pimpl_->F[predecessors.back() * matrix_width + j_div]);
          ++store_pos;
        }

        bool stop = false;
        for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
          if ((stop = F[j_mod] == H_pred[p * T::kNumVar + j_mod] + g_) ||
                      F[j_mod] == F_pred[p * T::kNumVar + j_mod] + e_) {
            prev_i = predecessors[p];
            break;
          }
        }

        alignment.emplace_back(rank_to_node[i - 1]->id, -1);
        i = prev_i;

        if (stop || i == 0) {
          break;
        }
      }
    }
  } while (true);

  delete[] backtrack_storage;

  // update alignment for NW (backtrack stops on first row or column)
  if (type_ == AlignmentType::kNW) {
    while (i == 0 && j != -1) {
      alignment.emplace_back(-1, j);
      --j;
    }
    while (i != 0 && j == -1) {
      alignment.emplace_back(rank_to_node[i - 1]->id, -1);

      const auto& it = rank_to_node[i - 1];
      if (it->inedges.empty()) {
        i = 0;
      } else {
        for (const auto& jt : it->inedges) {
          std::uint32_t pred_i = pimpl_->node_id_to_rank[jt->tail->id] + 1;
          if (pimpl_->first_column[i] == pimpl_->first_column[pred_i] + e_) {
            i = pred_i;
            break;
          }
        }
      }
    }
  }

  std::reverse(alignment.begin(), alignment.end());
  return alignment;
#else
  (void) sequence_len;
  (void) graph;
  (void) score;
  return Alignment();
#endif
}

template<Architecture A> template <typename T>
Alignment SimdAlignmentEngine<A>::Convex(
    std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score) noexcept {
#if defined(__AVX2__) || defined(__SSE4_1__) || defined(SPOA_USE_SIMDE)
  std::uint64_t normal_matrix_width = sequence_len;
  std::uint64_t matrix_width =
      std::ceil(static_cast<double>(sequence_len) / T::kNumVar);
  std::uint64_t matrix_height = graph.nodes().size() + 1;
  const auto& rank_to_node = graph.rank_to_node();

  typename T::type kNegativeInfinity =
      std::numeric_limits<typename T::type>::min() + 1024;

  typename T::type max_score = type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;  // NOLINT
  std::int32_t max_i = -1;
  std::int32_t max_j = -1;
  std::uint32_t last_column_id = (normal_matrix_width - 1) % T::kNumVar;
  __mxxxi zeroes = T::_mmxxx_set1_epi(0);
  __mxxxi g = T::_mmxxx_set1_epi(g_ - e_);
  __mxxxi e = T::_mmxxx_set1_epi(e_);
  __mxxxi q = T::_mmxxx_set1_epi(q_ - c_);
  __mxxxi c = T::_mmxxx_set1_epi(c_);

  __attribute__((aligned(kRegisterSize / 8))) typename T::type unpacked[T::kNumVar] = {0};  // NOLINT

  for (std::uint32_t i = 0, j = 0; i < T::kNumVar && j < T::kLogNumVar; ++i) {
    unpacked[i] = kNegativeInfinity;
    if ((i & (i + 1)) == 0) {
      pimpl_->masks[j++] = _mmxxx_load_si(
          reinterpret_cast<const __mxxxi*>(unpacked));
    }
  }
  pimpl_->masks[T::kLogNumVar] = _mmxxx_slli_si(
      T::_mmxxx_set1_epi(kNegativeInfinity),
      T::kLSS);

  pimpl_->penalties[0] = T::_mmxxx_set1_epi(e_);
  for (std::uint32_t i = 1; i < T::kLogNumVar; ++i) {
    pimpl_->penalties[i] = T::_mmxxx_add_epi(
        pimpl_->penalties[i - 1],
        pimpl_->penalties[i - 1]);
  }
  pimpl_->penalties[T::kLogNumVar] = T::_mmxxx_set1_epi(c_);
  for (std::uint32_t i = T::kLogNumVar + 1; i < 2 * T::kLogNumVar; ++i) {
    pimpl_->penalties[i] = T::_mmxxx_add_epi(
        pimpl_->penalties[i - 1],
        pimpl_->penalties[i - 1]);
  }

  // alignment
  for (const auto& it : rank_to_node) {
    __mxxxi* char_profile = &(pimpl_->sequence_profile[it->code * matrix_width]);  // NOLINT

    std::uint32_t i = pimpl_->node_id_to_rank[it->id] + 1;

    __mxxxi* H_row = &(pimpl_->H[i * matrix_width]);
    __mxxxi* F_row = &(pimpl_->F[i * matrix_width]);
    __mxxxi* O_row = &(pimpl_->O[i * matrix_width]);

    std::uint32_t pred_i = it->inedges.empty() ?
        0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

    __mxxxi* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);
    __mxxxi* F_pred_row = &(pimpl_->F[pred_i * matrix_width]);
    __mxxxi* O_pred_row = &(pimpl_->O[pred_i * matrix_width]);

    __mxxxi x = _mmxxx_srli_si(
        T::_mmxxx_set1_epi(pimpl_->first_column[pred_i]),
        T::kRSS);

    for (std::uint64_t j = 0; j < matrix_width; ++j) {
      // update F
      F_row[j] = T::_mmxxx_add_epi(
          T::_mmxxx_max_epi(
              T::_mmxxx_add_epi(H_pred_row[j], g),
              F_pred_row[j]),
          e);

      // update O
      O_row[j] = T::_mmxxx_add_epi(
          T::_mmxxx_max_epi(
              T::_mmxxx_add_epi(H_pred_row[j], q),
              O_pred_row[j]),
          c);

      // update H
      H_row[j] = T::_mmxxx_add_epi(
          _mmxxx_or_si(
              _mmxxx_slli_si(H_pred_row[j], T::kLSS),
              x),
          char_profile[j]);
      x = _mmxxx_srli_si(H_pred_row[j], T::kRSS);
    }
    // check other predecessors
    for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
      pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

      H_pred_row = &(pimpl_->H[pred_i * matrix_width]);
      F_pred_row = &(pimpl_->F[pred_i * matrix_width]);
      O_pred_row = &(pimpl_->O[pred_i * matrix_width]);

      x = _mmxxx_srli_si(
          T::_mmxxx_set1_epi(pimpl_->first_column[pred_i]),
          T::kRSS);

      for (std::uint64_t j = 0; j < matrix_width; ++j) {
        // update F
        F_row[j] = T::_mmxxx_max_epi(
            F_row[j],
            T::_mmxxx_add_epi(
                T::_mmxxx_max_epi(
                    T::_mmxxx_add_epi(H_pred_row[j], g),
                    F_pred_row[j]),
                e));

        // update O
        O_row[j] = T::_mmxxx_max_epi(
            O_row[j],
            T::_mmxxx_add_epi(
                T::_mmxxx_max_epi(
                    T::_mmxxx_add_epi(H_pred_row[j], q),
                    O_pred_row[j]),
                c));

        // update H
        H_row[j] = T::_mmxxx_max_epi(
            H_row[j],
            T::_mmxxx_add_epi(
                _mmxxx_or_si(
                    _mmxxx_slli_si(H_pred_row[j], T::kLSS),
                    x),
                char_profile[j]));

        x = _mmxxx_srli_si(H_pred_row[j], T::kRSS);
      }
    }

    __mxxxi* E_row = &(pimpl_->E[i * matrix_width]);
    __mxxxi* Q_row = &(pimpl_->Q[i * matrix_width]);

    x = T::_mmxxx_set1_epi(pimpl_->first_column[i]);
    __mxxxi y = T::_mmxxx_set1_epi(pimpl_->first_column[i]);

    __mxxxi score = zeroes;

    for (std::uint64_t j = 0; j < matrix_width; ++j) {
      H_row[j] = T::_mmxxx_max_epi(
          H_row[j],
          T::_mmxxx_max_epi(F_row[j], O_row[j]));

      E_row[j] = T::_mmxxx_add_epi(
          T::_mmxxx_add_epi(
              _mmxxx_or_si(
                  _mmxxx_slli_si(H_row[j], T::kLSS),
                  _mmxxx_srli_si(x, T::kRSS)),
              g),
          e);

      T::_mmxxx_prefix_max(E_row[j], pimpl_->masks, pimpl_->penalties);

      Q_row[j] = T::_mmxxx_add_epi(
          T::_mmxxx_add_epi(
              _mmxxx_or_si(
                  _mmxxx_slli_si(H_row[j], T::kLSS),
                  _mmxxx_srli_si(y, T::kRSS)),
              q),
          c);

      T::_mmxxx_prefix_max(Q_row[j], pimpl_->masks, &pimpl_->penalties[T::kLogNumVar]);  // NOLINT

      H_row[j] = T::_mmxxx_max_epi(
          H_row[j],
          T::_mmxxx_max_epi(E_row[j], Q_row[j]));

      x = T::_mmxxx_max_epi(
          H_row[j],
          T::_mmxxx_sub_epi(E_row[j], g));

      y = T::_mmxxx_max_epi(
          H_row[j],
          T::_mmxxx_sub_epi(Q_row[j], q));

      if (type_ == AlignmentType::kSW) {
        H_row[j] = T::_mmxxx_max_epi(H_row[j], zeroes);
      }
      score = T::_mmxxx_max_epi(score, H_row[j]);
    }

    if (type_ == AlignmentType::kSW) {
      std::int32_t max_row_score = _mmxxx_max_value<A, T>(score);
      if (max_score < max_row_score) {
        max_score = max_row_score;
        max_i = i;
      }
    } else if (type_ == AlignmentType::kOV) {
      if (it->outedges.empty()) {
        std::int32_t max_row_score = _mmxxx_max_value<A, T>(score);
        if (max_score < max_row_score) {
          max_score = max_row_score;
          max_i = i;
        }
      }
    } else if (type_ == AlignmentType::kNW) {
      if (it->outedges.empty()) {
        std::int32_t max_row_score = _mmxxx_value_at<A, T>(
            H_row[matrix_width - 1],
            last_column_id);
        if (max_score < max_row_score) {
          max_score = max_row_score;
          max_i = i;
        }
      }
    }
  }

  if (max_i == -1 && max_j == -1) {
    return Alignment();
  }
  if (score) {
    *score = max_score;
  }

  if (type_ == AlignmentType::kSW) {
    max_j = _mmxxx_index_of<A, T>(
        &(pimpl_->H[max_i * matrix_width]),
        matrix_width, max_score);
  } else if (type_ == AlignmentType::kOV) {
    if (rank_to_node[max_i - 1]->outedges.empty()) {
      max_j = _mmxxx_index_of<A, T>(
          &(pimpl_->H[max_i * matrix_width]),
          matrix_width, max_score);
    } else {
      max_j = normal_matrix_width - 1;
    }
  } else if (type_ == AlignmentType::kNW) {
    max_j = normal_matrix_width - 1;
  }

  // backtrack
  std::uint32_t max_num_predecessors = 1;
  for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(max_i); ++i) {
    max_num_predecessors = std::max(
        max_num_predecessors,
        static_cast<std::uint32_t>(rank_to_node[i]->inedges.size()));
  }

  typename T::type* backtrack_storage = nullptr;
  typename T::type* H = AllocateAlignedMemory<A>(
      &backtrack_storage,
      9 * T::kNumVar + 4 * T::kNumVar * max_num_predecessors,
      kRegisterSize / 8);
  typename T::type* H_pred = H + T::kNumVar;
  typename T::type* H_diag_pred = H_pred + T::kNumVar * max_num_predecessors;
  typename T::type* H_left = H_diag_pred + T::kNumVar * max_num_predecessors;
  typename T::type* F = H_left + T::kNumVar;
  typename T::type* F_pred = F + T::kNumVar;
  typename T::type* O = F_pred + T::kNumVar * max_num_predecessors;
  typename T::type* O_pred = O + T::kNumVar;
  typename T::type* E = O_pred + T::kNumVar * max_num_predecessors;
  typename T::type* E_left = E + T::kNumVar;
  typename T::type* Q = E_left + T::kNumVar;
  typename T::type* Q_left = Q + T::kNumVar;
  typename T::type* profile = Q_left + T::kNumVar;

  std::vector<std::uint32_t> predecessors;

  std::int32_t i = max_i;
  std::int32_t j = max_j;
  std::int32_t prev_i = 0, prev_j = 0;

  std::uint32_t j_div = j / T::kNumVar;
  std::uint32_t j_mod = j % T::kNumVar;

  bool load_next_segment = true;

  Alignment alignment;

  do {
    // check stop condition
    if (j == -1 || i == 0) {
      break;
    }

    const auto& it = rank_to_node[i - 1];
    // load everything
    if (load_next_segment) {
      predecessors.clear();

      // load current cells
      _mmxxx_store_si(
          reinterpret_cast<__mxxxi*>(H),
          pimpl_->H[i * matrix_width + j_div]);
      _mmxxx_store_si(
          reinterpret_cast<__mxxxi*>(E),
          pimpl_->E[i * matrix_width + j_div]);
      _mmxxx_store_si(
          reinterpret_cast<__mxxxi*>(Q),
          pimpl_->Q[i * matrix_width + j_div]);

      // load predecessors cells
      if (it->inedges.empty()) {
        predecessors.emplace_back(0);
        _mmxxx_store_si(reinterpret_cast<__mxxxi*>(H_pred), pimpl_->H[j_div]);
        _mmxxx_store_si(reinterpret_cast<__mxxxi*>(F_pred), pimpl_->F[j_div]);
        _mmxxx_store_si(reinterpret_cast<__mxxxi*>(O_pred), pimpl_->O[j_div]);
      } else {
        std::uint32_t store_pos = 0;
        for (const auto& jt : it->inedges) {
          predecessors.emplace_back(pimpl_->node_id_to_rank[jt->tail->id] + 1);
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(&H_pred[store_pos * T::kNumVar]),
              pimpl_->H[predecessors.back() * matrix_width + j_div]);
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(&F_pred[store_pos * T::kNumVar]),
              pimpl_->F[predecessors.back() * matrix_width + j_div]);
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(&O_pred[store_pos * T::kNumVar]),
              pimpl_->O[predecessors.back() * matrix_width + j_div]);
          ++store_pos;
        }
      }

      // load query profile cells
      _mmxxx_store_si(
          reinterpret_cast<__mxxxi*>(profile),
          pimpl_->sequence_profile[it->code * matrix_width + j_div]);
    }

    // check stop condition
    if (type_ == AlignmentType::kSW && H[j_mod] == 0) {
      break;
    }

    if (j_mod == 0) {
      // border case
      if (j_div > 0) {
        for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(&H_diag_pred[p * T::kNumVar]),
              pimpl_->H[predecessors[p] * matrix_width + (j_div - 1)]);
        }
        _mmxxx_store_si(
            reinterpret_cast<__mxxxi*>(H_left),
            pimpl_->H[i * matrix_width + j_div - 1]);
        _mmxxx_store_si(
            reinterpret_cast<__mxxxi*>(E_left),
            pimpl_->E[i * matrix_width + j_div - 1]);
        _mmxxx_store_si(
            reinterpret_cast<__mxxxi*>(Q_left),
            pimpl_->Q[i * matrix_width + j_div - 1]);
      } else {
        for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
          H_diag_pred[(p + 1) * T::kNumVar - 1] = pimpl_->first_column[predecessors[p]];  // NOLINT
        }
        H_left[T::kNumVar - 1] = pimpl_->first_column[i];
        E_left[T::kNumVar - 1] = pimpl_->first_column[i];
        Q_left[T::kNumVar - 1] = pimpl_->first_column[i];
      }
    }

    // find best predecessor cell
    bool predecessor_found = false, extend_left = false, extend_up = false;

    if (i != 0) {
      for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
        if ((j_mod == 0 && H[j_mod] == H_diag_pred[(p + 1) * T::kNumVar - 1] + profile[j_mod]) ||  // NOLINT
            (j_mod != 0 && H[j_mod] == H_pred[p * T::kNumVar + j_mod - 1] + profile[j_mod])) {  // NOLINT
          prev_i = predecessors[p];
          prev_j = j - 1;
          predecessor_found = true;
          break;
        }
      }
    }

    if (!predecessor_found && i != 0) {
      for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
        if ((extend_up = H[j_mod] == F_pred[p * T::kNumVar + j_mod] + e_) ||
                         H[j_mod] == H_pred[p * T::kNumVar + j_mod] + g_  ||
            (extend_up = H[j_mod] == O_pred[p * T::kNumVar + j_mod] + c_) ||
                         H[j_mod] == H_pred[p * T::kNumVar + j_mod] + q_) {
          prev_i = predecessors[p];
          prev_j = j;
          predecessor_found = true;
          break;
        }
      }
    }

    if (!predecessor_found) {
      if ((j_mod != 0 && ((extend_left = H[j_mod] == E[j_mod - 1] + e_)  ||
                                         H[j_mod] == H[j_mod - 1] + g_   ||
                          (extend_left = H[j_mod] == Q[j_mod - 1] + c_)  ||
                                         H[j_mod] == H[j_mod - 1] + q_)) ||
          (j_mod == 0 && ((extend_left = H[j_mod] == E_left[T::kNumVar - 1] + e_) ||  // NOLINT
                                         H[j_mod] == H_left[T::kNumVar - 1] + g_  ||  // NOLINT
                          (extend_left = H[j_mod] == Q_left[T::kNumVar - 1] + c_) ||  // NOLINT
                                         H[j_mod] == H_left[T::kNumVar - 1] + q_))) {  // NOLINT
        prev_i = i;
        prev_j = j - 1;
        predecessor_found = true;
      }
    }

    alignment.emplace_back(
        i == prev_i ? -1 : rank_to_node[i - 1]->id,
        j == prev_j ? -1 : j);

    // update for next round
    load_next_segment =
        (i == prev_i ? false : true) ||
        (j != prev_j && prev_j % T::kNumVar == T::kNumVar - 1 ? true : false);

    i = prev_i;
    j = prev_j;
    j_div = j / T::kNumVar;
    j_mod = j % T::kNumVar;

    if (extend_left) {
      while (true) {
        // load
        if (j_mod == T::kNumVar - 1) {
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(E),
              pimpl_->E[i * matrix_width + j_div]);
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(Q),
              pimpl_->Q[i * matrix_width + j_div]);
        } else if (j_mod == 0) {  // boarder case
          if (j_div > 0) {
            _mmxxx_store_si(
                reinterpret_cast<__mxxxi*>(E_left),
                pimpl_->E[i * matrix_width + j_div - 1]);
            _mmxxx_store_si(
                reinterpret_cast<__mxxxi*>(Q_left),
                pimpl_->Q[i * matrix_width + j_div - 1]);
          }
        }

        alignment.emplace_back(-1, j);
        --j;
        j_div = j / T::kNumVar;
        j_mod = j % T::kNumVar;
        if ((j == -1) ||
            (j_mod != T::kNumVar - 1 &&      E[j_mod] + e_ != E[j_mod + 1]) ||
            (j_mod == T::kNumVar - 1 && E_left[j_mod] + e_ != E[0])         ||
            (j_mod != T::kNumVar - 1 &&      Q[j_mod] + c_ != Q[j_mod + 1]) ||
            (j_mod == T::kNumVar - 1 && Q_left[j_mod] + c_ != Q[0])) {
          break;
        }
      }
      load_next_segment = true;
    } else if (extend_up) {
      while (true) {
        // load
        _mmxxx_store_si(
            reinterpret_cast<__mxxxi*>(F),
            pimpl_->F[i * matrix_width + j_div]);
        _mmxxx_store_si(
            reinterpret_cast<__mxxxi*>(O),
            pimpl_->O[i * matrix_width + j_div]);

        predecessors.clear();
        std::uint32_t store_pos = 0;
        for (const auto& it : rank_to_node[i - 1]->inedges) {
          predecessors.emplace_back(pimpl_->node_id_to_rank[it->tail->id] + 1);
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(&H_pred[store_pos * T::kNumVar]),
              pimpl_->H[predecessors.back() * matrix_width + j_div]);
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(&F_pred[store_pos * T::kNumVar]),
              pimpl_->F[predecessors.back() * matrix_width + j_div]);
          _mmxxx_store_si(
              reinterpret_cast<__mxxxi*>(&O_pred[store_pos * T::kNumVar]),
              pimpl_->O[predecessors.back() * matrix_width + j_div]);
          ++store_pos;
        }

        bool stop = true;
        prev_i = 0;
        for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
          if (F[j_mod] == F_pred[p * T::kNumVar + j_mod] + e_ ||
              O[j_mod] == O_pred[p * T::kNumVar + j_mod] + c_) {
            prev_i = predecessors[p];
            stop = false;
            break;
          }
        }
        if (stop == true) {
          for (std::uint32_t p = 0; p < predecessors.size(); ++p) {
            if (F[j_mod] == H_pred[p * T::kNumVar + j_mod] + g_ ||
                O[j_mod] == H_pred[p * T::kNumVar + j_mod] + q_) {
              prev_i = predecessors[p];
              break;
            }
          }
        }

        alignment.emplace_back(rank_to_node[i - 1]->id, -1);
        i = prev_i;

        if (stop || i == 0) {
          break;
        }
      }
    }
  } while (true);

  delete[] backtrack_storage;

  // update alignment for NW (backtrack stops on first row or column)
  if (type_ == AlignmentType::kNW) {
    while (i == 0 && j != -1) {
      alignment.emplace_back(-1, j);
      --j;
    }
    while (i != 0 && j == -1) {
      alignment.emplace_back(rank_to_node[i - 1]->id, -1);

      const auto& it = rank_to_node[i - 1];
      if (it->inedges.empty()) {
        i = 0;
      } else {
        for (const auto& jt : it->inedges) {
          std::uint32_t pred_i = pimpl_->node_id_to_rank[jt->tail->id] + 1;
          if (pimpl_->first_column[matrix_height + i]     == pimpl_->first_column[matrix_height + pred_i] + e_ ||  // NOLINT
              pimpl_->first_column[2 * matrix_height + i] == pimpl_->first_column[2 * matrix_height + pred_i] + c_ ) {  // NOLINT
            i = pred_i;
            break;
          }
        }
      }
    }
  }

  std::reverse(alignment.begin(), alignment.end());
  return alignment;
#else
  (void) sequence_len;
  (void) graph;
  (void) score;
  return Alignment();
#endif
}

}  // namespace spoa

#endif  // SIMD_ALIGNMENT_ENGINE_IMPLEMENTATION_HPP_
