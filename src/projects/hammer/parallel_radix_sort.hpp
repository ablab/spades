//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

// Copyright 2010, Takuya Akiba
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Takuya Akiba nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef PARALLEL_RADIX_SORT_H_
#define PARALLEL_RADIX_SORT_H_

#include "utils/parallel/openmp_wrapper.h"

#include <stdint.h>
#include <cstring>
#include <cassert>
#include <climits>
#include <algorithm>
#include <utility>

namespace parallel_radix_sort {

namespace internal {
// Size of the software managed buffer
const size_t kOutBufferSize = 32;

// The algorithm is implemented in this internal class
template<typename PlainType, typename UnsignedType, typename Encoder,
         typename ValueManager, int Base>
class ParallelRadixSortInternal {
public:
  ParallelRadixSortInternal();
  ~ParallelRadixSortInternal();

  void Init(size_t max_elems, int max_threads);

  PlainType *Sort(PlainType *data, size_t num_elems, int num_threads,
                  ValueManager *value_manager);

  static void InitAndSort(PlainType *data, size_t num_elems, int num_threads,
                          ValueManager *value_manager);
private:
  size_t max_elems_;
  int max_threads_;

  UnsignedType *tmp_;
  size_t **histo_;
  UnsignedType ***out_buf_;
  size_t **out_buf_n_;

  int num_threads_;
  size_t *pos_bgn_, *pos_end_;
  ValueManager *value_manager_;

  void DeleteAll();

  UnsignedType *SortInternal(UnsignedType *data, size_t num_elems,
                             int num_threads, ValueManager *value_manager);

  // Compute |pos_bgn_| and |pos_end_| (associated ranges for each threads)
  void ComputeRanges(size_t num_elems);

  // First step of each iteration of sorting
  // Compute the histogram of |src| using bits in [b, b + Base)
  void ComputeHistogram(int b, UnsignedType *src);

  // Second step of each iteration of sorting
  // Scatter elements of |src| to |dst| using the histogram
  void Scatter(int b, UnsignedType *src, UnsignedType *dst);
};

template<typename PlainType, typename UnsignedType, typename Encoder,
         typename ValueManager, int Base>
ParallelRadixSortInternal<PlainType, UnsignedType, Encoder, ValueManager, Base>
::ParallelRadixSortInternal()
  : max_elems_(0), max_threads_(0), tmp_(NULL), histo_(NULL),
    out_buf_(NULL), out_buf_n_(NULL), pos_bgn_(NULL), pos_end_(NULL) {
  assert(sizeof(PlainType) == sizeof(UnsignedType));
}

template<typename PlainType, typename UnsignedType, typename Encoder,
         typename ValueManager, int Base>
ParallelRadixSortInternal
<PlainType, UnsignedType, Encoder, ValueManager, Base>
::~ParallelRadixSortInternal() {
  DeleteAll();
}

template<typename PlainType, typename UnsignedType, typename Encoder,
         typename ValueManager, int Base>
void ParallelRadixSortInternal
<PlainType, UnsignedType, Encoder, ValueManager, Base>
::DeleteAll() {
  delete [] tmp_;
  tmp_ = NULL;

  for (int i = 0; i < max_threads_; ++i) delete [] histo_[i];
  delete [] histo_;
  histo_ = NULL;

  for (int i = 0; i < max_threads_; ++i) {
    for (size_t j = 0; j < 1 << Base; ++j) {
      delete [] out_buf_[i][j];
    }
    delete [] out_buf_n_[i];
    delete [] out_buf_[i];
  }
  delete [] out_buf_;
  delete [] out_buf_n_;
  out_buf_ = NULL;
  out_buf_n_ = NULL;

  delete [] pos_bgn_;
  delete [] pos_end_;
  pos_bgn_ = pos_end_ = NULL;

  max_elems_ = 0;
  max_threads_ = 0;
}

template<typename PlainType, typename UnsignedType, typename Encoder,
         typename ValueManager, int Base>
void ParallelRadixSortInternal
<PlainType, UnsignedType, Encoder, ValueManager, Base>
::Init(size_t max_elems, int max_threads) {
  DeleteAll();

  max_elems_ = max_elems;

  if (max_threads == -1) {
    max_threads = omp_get_max_threads();
  }
  assert(max_threads >= 1);
  max_threads_ = max_threads;

  tmp_ = new UnsignedType[max_elems];
  histo_ = new size_t*[max_threads];
  for (int i = 0; i < max_threads; ++i) {
    histo_[i] = new size_t[1 << Base];
  }

  out_buf_ = new UnsignedType**[max_threads];
  out_buf_n_ = new size_t*[max_threads];
  for (int i = 0; i < max_threads; ++i) {
    out_buf_[i] = new UnsignedType*[1 << Base];
    out_buf_n_[i] = new size_t[1 << Base];
    for (size_t j = 0; j < 1 << Base; ++j) {
      out_buf_[i][j] = new UnsignedType[kOutBufferSize];
    }
  }

  pos_bgn_ = new size_t[max_threads];
  pos_end_ = new size_t[max_threads];
}

template<typename PlainType, typename UnsignedType, typename Encoder,
         typename ValueManager, int Base>
PlainType *ParallelRadixSortInternal
<PlainType, UnsignedType, Encoder, ValueManager, Base>
::Sort(PlainType *data, size_t num_elems,
       int num_threads, ValueManager *value_manager) {
  UnsignedType *src = reinterpret_cast<UnsignedType*>(data);
  UnsignedType *res = SortInternal(src, num_elems, num_threads, value_manager);
  return reinterpret_cast<PlainType*>(res);
}

template<typename PlainType, typename UnsignedType, typename Encoder,
         typename ValueManager, int Base>
void ParallelRadixSortInternal
<PlainType, UnsignedType, Encoder, ValueManager, Base>
::InitAndSort(PlainType *data, size_t num_elems,
              int num_threads, ValueManager *value_manager) {
  ParallelRadixSortInternal prs;
  prs.Init(num_elems, num_threads);
  const PlainType *res = prs.Sort(data, num_elems, num_threads, value_manager);
  if (res != data) {
    for (size_t i = 0; i < num_elems; ++i) data[i] = res[i];
  }
}

template<typename PlainType, typename UnsignedType, typename Encoder,
         typename ValueManager, int Base>
UnsignedType *ParallelRadixSortInternal
<PlainType, UnsignedType, Encoder, ValueManager, Base>
::SortInternal(UnsignedType *data, size_t num_elems,
               int num_threads, ValueManager *value_manager) {
  assert(num_elems <= max_elems_);

  if (num_threads == -1) {
    num_threads = omp_get_max_threads();
  }
  assert(1 <= num_threads && num_threads <= max_threads_);
  num_threads_ = num_threads;

  value_manager_ = value_manager;

  // Compute |pos_bgn_| and |pos_end_|
  ComputeRanges(num_elems);

  // Iterate from lower bits to higher bits
  const unsigned bits = CHAR_BIT * sizeof(UnsignedType);
  UnsignedType *src = data, *dst = tmp_;
  for (unsigned b = 0; b < bits; b += Base) {
    ComputeHistogram(b, src);
    Scatter(b, src, dst);

    std::swap(src, dst);
    value_manager->Next();
  }

  return src;
}

template<typename PlainType, typename UnsignedType, typename Encoder,
         typename ValueManager, int Base>
void ParallelRadixSortInternal
<PlainType, UnsignedType, Encoder, ValueManager, Base>
::ComputeRanges(size_t num_elems) {
  pos_bgn_[0] = 0;
  for (int i = 0; i < num_threads_ - 1; ++i) {
    const size_t t = (num_elems - pos_bgn_[i]) / (num_threads_ - i);
    pos_bgn_[i + 1] = pos_end_[i] = pos_bgn_[i] + t;
  }
  pos_end_[num_threads_ - 1] = num_elems;
}

template<typename PlainType, typename UnsignedType, typename Encoder,
         typename ValueManager, int Base>
void ParallelRadixSortInternal
<PlainType, UnsignedType, Encoder, ValueManager, Base>
::ComputeHistogram(int b, UnsignedType *src) {
  // Compute local histogram
  #ifdef _OPENMP
  #pragma omp parallel num_threads(num_threads_)
  #endif
  {
    const int my_id = omp_get_thread_num();
    const size_t my_bgn = pos_bgn_[my_id];
    const size_t my_end = pos_end_[my_id];
    size_t *my_histo = histo_[my_id];

    memset(my_histo, 0, sizeof(size_t) * (1 << Base));
    for (size_t i = my_bgn; i < my_end; ++i) {
      __builtin_prefetch(src + i + 1, 0, 1);
      size_t t = Encoder::extract(src[i], b, Base);
      ++my_histo[t];
    }
  }

  // Compute global histogram
  size_t s = 0;
  for (size_t i = 0; i < 1 << Base; ++i) {
    for (int j = 0; j < num_threads_; ++j) {
      const size_t t = s + histo_[j][i];
      histo_[j][i] = s;
      s = t;
    }
  }
}

template<typename PlainType, typename UnsignedType, typename Encoder,
         typename ValueManager, int Base>
void ParallelRadixSortInternal
<PlainType, UnsignedType, Encoder, ValueManager, Base>
::Scatter(int b, UnsignedType *src, UnsignedType *dst) {
  #ifdef _OPENMP
  #pragma omp parallel num_threads(num_threads_)
  #endif
  {
    const int my_id = omp_get_thread_num();
    const size_t my_bgn = pos_bgn_[my_id];
    const size_t my_end = pos_end_[my_id];
    size_t *my_histo = histo_[my_id];
    UnsignedType **my_buf = out_buf_[my_id];
    size_t *my_buf_n = out_buf_n_[my_id];

    memset(my_buf_n, 0, sizeof(size_t) * (1 << Base));
    for (size_t i = my_bgn; i < my_end; ++i) {
      __builtin_prefetch(src + i + 1, 0, 1);

      size_t t = Encoder::extract(src[i], b, Base);
      my_buf[t][my_buf_n[t]] = src[i];
      value_manager_->Push(my_id, t, my_buf_n[t], i);
      ++my_buf_n[t];

      if (my_buf_n[t] == kOutBufferSize) {
        size_t p = my_histo[t];
        for (size_t j = 0; j < kOutBufferSize; ++j) {
          dst[p++] = my_buf[t][j];
        }
        value_manager_->Flush(my_id, t, kOutBufferSize, my_histo[t]);

        my_histo[t] += kOutBufferSize;
        my_buf_n[t] = 0;
      }
    }

    // Flush everything
    for (size_t i = 0; i < 1 << Base; ++i) {
      size_t p = my_histo[i];
      for (size_t j = 0; j < my_buf_n[i]; ++j) {
        dst[p++] = my_buf[i][j];
      }
      value_manager_->Flush(my_id, i, my_buf_n[i], my_histo[i]);
    }
  }
}
}  // namespace internal

// Encoders encode signed/unsigned integers and floating point numbers
// to correctly ordered unsigned integers
namespace encoder {
class EncoderUnsigned {
public:
  template<typename UnsignedType>
  inline static size_t extract(const UnsignedType &x, unsigned shift, unsigned Base) {
    return (x >> shift) & ((1 << Base) - 1);
  }
};

class EncoderSigned {
public:
  template<typename UnsignedType>
  inline static size_t extract(const UnsignedType &x, unsigned shift, unsigned Base) {
    x = x ^ (UnsignedType(1) << (CHAR_BIT * sizeof(UnsignedType) - 1));
    return (x >> shift) & ((1 << Base) - 1);
  }
};

class EncoderDecimal {
public:
  template<typename UnsignedType>
  inline static size_t extract(const UnsignedType &x, unsigned shift, unsigned Base) {
    static const int bits = CHAR_BIT * sizeof(UnsignedType);
    const UnsignedType a = x >> (bits - 1);
    const UnsignedType b = (-a) | (UnsignedType(1) << (bits - 1));
    x = x ^ b;
    return (x >> shift) & ((1 << Base) - 1);
  }
};
}  // namespace encoder

// Value managers are used to generalize the sorting algorithm
// to sorting of keys and sorting of pairs
namespace value_manager {
class DummyValueManager {
public:
  inline void Push(int thread __attribute__((unused)),
                   size_t bucket __attribute__((unused)),
                   size_t num __attribute__((unused)),
                   size_t from_pos __attribute__((unused))) {}

  inline void Flush(int thread __attribute__((unused)),
                    size_t bucket __attribute__((unused)),
                    size_t num __attribute__((unused)),
                    size_t to_pos __attribute__((unused))) {}

  void Next() {}
};

template<typename ValueType, int Base> class PairValueManager {
public:
  PairValueManager()
    : max_elems_(0), max_threads_(0), original_(NULL), tmp_(NULL),
      src_(NULL), dst_(NULL), out_buf_(NULL) {}

  ~PairValueManager() {
    DeleteAll();
  }

  void Init(size_t max_elems, int max_threads);

  void Start(ValueType *original, size_t num_elems, int num_threads) {
    assert(num_elems <= max_elems_);
    assert(num_threads <= max_threads_);
    src_ = original_ = original;
    dst_ = tmp_;
  }

  inline void Push(int thread, size_t bucket, size_t num, size_t from_pos) {
    out_buf_[thread][bucket][num] = src_[from_pos];
  }

  inline void Flush(int thread, size_t bucket, size_t num, size_t to_pos) {
    for (size_t i = 0; i < num; ++i) {
      dst_[to_pos++] = out_buf_[thread][bucket][i];
    }
  }

  void Next() {
    std::swap(src_, dst_);
  }

  ValueType *GetResult() {
    return src_;
  }
private:
  size_t max_elems_;
  int max_threads_;

  static const size_t kOutBufferSize = internal::kOutBufferSize;
  ValueType *original_, *tmp_;
  ValueType *src_, *dst_;
  ValueType ***out_buf_;

  void DeleteAll();
};

template<typename ValueType, int Base>
void PairValueManager<ValueType, Base>
::Init(size_t max_elems, int max_threads) {
  if (max_threads == -1) {
    max_threads = omp_get_max_threads();
  }
  assert(max_threads >= 1);

  DeleteAll();

  max_elems_ = max_elems;
  max_threads_ = max_threads;

  tmp_ = new ValueType[max_elems];

  out_buf_ = new ValueType**[max_threads];
  for (int i = 0; i < max_threads; ++i) {
    out_buf_[i] = new ValueType*[1 << Base];
    for (size_t j = 0; j < 1 << Base; ++j) {
      out_buf_[i][j] = new ValueType[kOutBufferSize];
    }
  }
}

template<typename ValueType, int Base>
void PairValueManager<ValueType, Base>
::DeleteAll() {
  delete [] tmp_;
  tmp_ = NULL;

  for (int i = 0; i < max_threads_; ++i) {
    for (size_t j = 0; j < 1 << Base; ++j) {
      delete [] out_buf_[i][j];
    }
    delete [] out_buf_[i];
  }
  delete [] out_buf_;
  out_buf_ = NULL;

  max_elems_ = 0;
  max_threads_ = 0;
}
}  // namespace value_manager

// Frontend class for sorting keys
template<typename PlainType, typename UnsignedType = PlainType,
         typename Encoder = encoder::EncoderUnsigned, int Base = 8>
class KeySort {
  typedef value_manager::DummyValueManager DummyValueManager;
  typedef internal::ParallelRadixSortInternal
  <PlainType, UnsignedType, Encoder, DummyValueManager, Base> Internal;

public:
  // In the following functions, when |max_threads| or |num_threads| is -1,
  // the default value given by OpenMP would be used.
  void Init(size_t max_elems, int max_threads = -1) {
    internal_.Init(max_elems, max_threads);
  }

  // Notice that the pointer returned by this
  // does not necessarily equal to |data|.
  PlainType *Sort(PlainType *data, size_t num_elems, int num_threads = -1) {
    return internal_.Sort(data, num_elems, num_threads, &dummy_value_manager_);
  }

  static void InitAndSort(PlainType *data, size_t num_elems, int num_threads = -1) {
    DummyValueManager dvm;
    Internal::InitAndSort(data, num_elems, num_threads, &dvm);
  }
private:
  Internal internal_;
  DummyValueManager dummy_value_manager_;
};

// Frontend class for sorting pairs
template<typename PlainType, typename ValueType,
         typename UnsignedType = PlainType,
         typename Encoder = encoder::EncoderUnsigned,
         int Base = 8>
class PairSort {
  typedef value_manager::PairValueManager
  <ValueType, Base> ValueManager;
  typedef internal::ParallelRadixSortInternal
  <PlainType, UnsignedType, Encoder, ValueManager, Base> Internal;

public:
  // In the following functions, when |max_threads| or |num_threads| is -1,
  // the default value given by OpenMP would be used.
  void Init(size_t max_elems, int max_threads = -1) {
    internal_.Init(max_elems, max_threads);
    value_manager_.Init(max_elems, max_threads);
  }

  // Notice that the pointers returned by this
  // do not necessarily equal to |keys| and |vals|.
  std::pair<PlainType*, ValueType*> Sort(PlainType *keys, ValueType *vals,
                                         size_t num_elems, int num_threads = -1) {
    value_manager_.Start(vals, num_elems, num_threads);
    PlainType *res_keys = internal_.Sort(keys, num_elems, num_threads, &value_manager_);
    ValueType *res_vals = value_manager_.GetResult();
    return std::make_pair(res_keys, res_vals);
  }

  static void InitAndSort(PlainType *keys, ValueType *vals,
                          size_t num_elems, int num_threads = -1) {
    ValueManager vm;
    vm.Init(num_elems, num_threads);
    vm.Start(vals, num_elems, num_threads);
    Internal::InitAndSort(keys, num_elems, num_threads, &vm);
    ValueType *res_vals = vm.GetResult();
    if (res_vals != vals) {
      for (size_t i = 0; i < num_elems; ++i) {
        vals[i] = res_vals[i];
      }
    }
  }
private:
  Internal internal_;
  ValueManager value_manager_;
};

#define TYPE_CASE(plain_type, unsigned_type, encoder_type)           \
  template<> class KeySort<plain_type>                               \
    : public KeySort<plain_type, unsigned_type,                      \
                     encoder::Encoder ## encoder_type> {};           \
  template<typename V> class PairSort<plain_type, V>                 \
    : public PairSort<plain_type, V, unsigned_type,                  \
                      encoder::Encoder ## encoder_type> {};          \

// Signed integers
TYPE_CASE(char,        unsigned char,      Signed);
TYPE_CASE(short,       unsigned short,     Signed);
TYPE_CASE(int,         unsigned int,       Signed);
TYPE_CASE(long,        unsigned long,      Signed);
TYPE_CASE(long long,   unsigned long long, Signed);

// |signed char| and |char| are treated as different types
TYPE_CASE(signed char, unsigned char,      Signed);

// Floating point numbers
TYPE_CASE(float,  uint32_t, Decimal);
TYPE_CASE(double, uint64_t, Decimal);

#undef TYPE_CASE

template<typename KeyType>
void SortKeys(KeyType *data, size_t num_elems, int num_threads = -1) {
  KeySort<KeyType>::InitAndSort(data, num_elems, num_threads);
}

template<typename KeyType, typename ValueType>
void SortPairs(KeyType *keys, ValueType *vals, size_t num_elems, int num_threads = -1) {
  PairSort<KeyType, ValueType>::InitAndSort(keys, vals, num_elems, num_threads);
}
};  // namespace parallel radix sort

#endif  // PARALLEL_RADIX_SORT_H_
