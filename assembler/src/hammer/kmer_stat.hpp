//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef HAMMER_KMERSTAT_HPP_
#define HAMMER_KMERSTAT_HPP_

#include "verify.hpp"

#include "sequence/seq.hpp"

#include <functional>
#include <vector>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <cstdint>
#include <cmath>

#include <sched.h>
#include <string.h>


namespace hammer {
const uint32_t K = 21;
typedef Seq<K> KMer;
};

class Read;
struct KMerStat;

static inline unsigned hamdistKMer(const hammer::KMer &x, const hammer::KMer &y,
                                   unsigned tau = hammer::K) {
  unsigned dist = 0;
  for (unsigned i = 0; i < hammer::K; ++i) {
    if (x[i] != y[i]) {
      ++dist; if (dist > tau) return dist;
    }
  }
  return dist;
}

struct QualBitSet {
  unsigned char q_[hammer::K];

  QualBitSet(const unsigned char *q = NULL) {
    if (q)
      memcpy(q_, q, hammer::K);
    else
      memset(q_, 0, hammer::K);
  }

  // QualBitSet is POD-like object, use default stuff.
  QualBitSet(const QualBitSet &qbs) = default;
  QualBitSet& operator=(const QualBitSet &qbs) = default;
  // Workaround gcc bug, it cannot generate default move ctor...
  QualBitSet& operator=(QualBitSet &&qbs) {
    memcpy(q_, qbs.q_, hammer::K);

    return *this;
  }

  QualBitSet& operator+=(const unsigned char *data) {
    for (size_t i = 0; i < hammer::K; ++i)
      q_[i] = (unsigned char)std::min(255, data[i] + q_[i]);

    return *this;
  }

  QualBitSet& operator+=(const QualBitSet &qbs) {
    const unsigned char* data = qbs.q_;
    return this->operator+=(data);
  }

  unsigned short operator[](size_t n) const {
    return (unsigned short)q_[n];
  }

  unsigned short at(size_t n) const {
    // FIXME: Bound checking
    return (unsigned short)q_[n];
  }

  void set(size_t n, unsigned short value) {
    q_[n] = (unsigned char)value;
  }

  void set(const char *value) {
    memcpy(q_, value, hammer::K);
  }
};

struct KMerStat {
  enum {
    Change             = 0,
    Bad                = 1,
    Good               = 2,
    GoodIter           = 3,
    GoodIterBad        = 4,
    MarkedForGoodIter  = 5
  } KMerStatus;

  KMerStat(uint32_t cnt, float kquality, const unsigned char *quality) : totalQual(kquality), count(cnt), qual(quality), lock_(0) {
    __sync_lock_release(&lock_);
  }
  KMerStat() : totalQual(1.0), count(0), qual(), lock_(0) {
    __sync_lock_release(&lock_);
  }

  union {
    struct {
      uint64_t changeto : 48;
      unsigned status : 3;
      unsigned res    : 13;
    };
    uint64_t raw_data;
  };

  float totalQual;
  uint32_t count;
  QualBitSet qual;
  uint8_t lock_; // FIXME: Turn into single bit of status field.

  void lock() {
    while (__sync_val_compare_and_swap(&lock_, 0, 1) == 1)
      sched_yield();
  }
  void unlock() {
    lock_ = 0;
    __sync_synchronize();
  }

  bool isGood() const { return status >= Good; }
  bool isGoodForIterative() const { return (status == GoodIter); }
  void makeGoodForIterative() { status = GoodIter; }
  void markGoodForIterative() { status = MarkedForGoodIter; }
  bool isMarkedGoodForIterative() { return (status == MarkedForGoodIter); }
  bool change() const { return status == Change; }
  void set_change(uint64_t kmer) {
    changeto = kmer & (((uint64_t) 1 << 48) - 1);
    status = Change;
  }
};

inline
std::ostream& operator<<(std::ostream &os, const KMerStat &kms) {
  os << /* kms.kmer().str() << */ " (" << std::setw(3) << kms.count << ", " << std::setprecision(6) << std::setw(8) << (1-kms.totalQual) << ')';

  return os;
}

template<class Writer>
inline Writer& binary_write(Writer &os, const QualBitSet &qbs) {
  os.write((char*)&qbs.q_[0], sizeof(qbs.q_));

  return os;
}

template<class Reader>
inline void binary_read(Reader &is, QualBitSet &qbs) {
  is.read((char*)&qbs.q_[0], sizeof(qbs.q_));
}

template<class Writer>
inline Writer& binary_write(Writer &os, const KMerStat &k) {
  os.write((char*)&k.count, sizeof(k.count));
  os.write((char*)&k.raw_data, sizeof(k.raw_data));
  os.write((char*)&k.totalQual, sizeof(k.totalQual));
  return binary_write(os, k.qual);
}

template<class Reader>
inline void binary_read(Reader &is, KMerStat &k) {
  is.read((char*)&k.count, sizeof(k.count));
  is.read((char*)&k.raw_data, sizeof(k.raw_data));
  is.read((char*)&k.totalQual, sizeof(k.totalQual));
  binary_read(is, k.qual);
}

inline unsigned char getQual(const KMerStat & kmc, size_t i) {
  return (unsigned char)kmc.qual[i];
}

inline double getProb(const KMerStat &kmc, size_t i, bool log);
inline double getRevProb(const KMerStat &kmc, size_t i, bool log);

namespace hammer {
typedef std::array<char, hammer::K> ExpandedSeq;

static inline unsigned hamdist(const ExpandedSeq &x, const ExpandedSeq &y,
                               unsigned tau = hammer::K) {
  unsigned dist = 0;
  for (unsigned i = 0; i < hammer::K; ++i) {
    if (x[i] != y[i]) {
      ++dist; if (dist > tau) return dist;
    }
  }
  return dist;
}

class ExpandedKMer {
 public:
  ExpandedKMer(const KMer k, const KMerStat &kmc) {
    for (unsigned i = 0; i < hammer::K; ++i) {
      s_[i] = k[i];
      for (unsigned j = 0; j < 4; ++j)
        lprobs_[4*i + j] = ((char)j != s_[i] ?
                            getRevProb(kmc, i, /* log */ true) - log(3) :
                            getProb(kmc, i, /* log */ true));
    }
    count_ = kmc.count;
  }

  double logL(const ExpandedSeq &center) const {
    double res = 0;
    for (unsigned i = 0; i < hammer::K; ++i)
      res += lprobs_[4*i + center[i]];

    return res;
  }

  double logL(const ExpandedKMer &center) const {
    return logL(center.s_);
  }

  unsigned hamdist(const ExpandedSeq &k,
                   unsigned tau = hammer::K) const {
    unsigned dist = 0;
    for (unsigned i = 0; i < hammer::K; ++i) {
      if (s_[i] != k[i]) {
        ++dist; if (dist > tau) return dist;
      }
    }

    return dist;
  }

  unsigned hamdist(const ExpandedKMer &k,
                   unsigned tau = hammer::K) const {
    return hamdist(k.s_, tau);
  }

  double logL(const KMer center) const {
    double res = 0;
    for (unsigned i = 0; i < hammer::K; ++i)
      res += lprobs_[4*i + center[i]];

    return res;
  }

  unsigned hamdist(const KMer &k,
                   unsigned tau = hammer::K) const {
    unsigned dist = 0;
    for (unsigned i = 0; i < hammer::K; ++i) {
      if (s_[i] != k[i]) {
        ++dist; if (dist > tau) return dist;
      }
    }

    return dist;
  }

  uint32_t count() const {
    return count_;
  }

  ExpandedSeq seq() const {
    return s_;
  }

 private:
  double lprobs_[4*hammer::K];
  uint32_t count_;
  ExpandedSeq s_;
};

inline
std::ostream& operator<<(std::ostream &os, const ExpandedSeq &seq) {
  for (auto s : seq)
    os << nucl(s);

  return os;
}

};

#endif //  HAMMER_KMERSTAT_HPP_
