//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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

#include <sched.h>
#include <string.h>
#include <stdint.h>


namespace hammer {
const uint32_t K = 21;
const uint32_t M = 55;
typedef Seq<K> KMer;
};

typedef uint64_t hint_t;


class Read;
class PositionRead;
class PositionKMer;
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

typedef std::map<PositionKMer, KMerStat> KMerStatMap;
typedef std::pair<hammer::KMer, std::pair<uint32_t, double> > StringCount;

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
      q_[i] = std::min(255, data[i] + q_[i]);

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

  KMerStat(uint32_t cnt, hammer::KMer k, float kquality, const unsigned char *quality) : kmer_(k), totalQual(kquality), count(cnt), qual(quality), lock_(0) {
    __sync_lock_release(&lock_);
  }
  KMerStat() : totalQual(1.0), count(0), qual(), lock_(0) {
    __sync_lock_release(&lock_);
  }

  union {
    struct {
      hint_t changeto : 48;
      unsigned status : 3;
      unsigned res    : 13;
    };
    uint64_t raw_data;
  };

  hammer::KMer kmer_;
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
  void set_change(hint_t kmer) {
    changeto = kmer;
    status = Change;
  }
  const hammer::KMer& kmer() const { return kmer_; }
};

inline
std::ostream& operator<<(std::ostream &os, const KMerStat &kms) {
  os << kms.kmer().str() << " (" << std::setw(3) << kms.count << ", " << std::setprecision(6) << std::setw(8) << (1-kms.totalQual) << ')';

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
  return kmc.qual[i];
}

#endif //  HAMMER_KMERSTAT_HPP_
