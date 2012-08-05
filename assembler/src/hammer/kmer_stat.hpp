//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef HAMMER_KMERSTAT_HPP_
#define HAMMER_KMERSTAT_HPP_

#include "verify.hpp"

#include <stdint.h>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <string.h>

const uint32_t K = 21;
const uint32_t M = 55;
typedef uint64_t hint_t;

#define BLOBKMER_UNDEFINED            2e12
#define KMERSTAT_CHANGE               2e12
#define KMERSTAT_BAD                  2e12 + 1
#define KMERSTAT_GOOD                 2e12 + 2
#define KMERSTAT_GOODITER             2e12 + 3
#define KMERSTAT_GOODITER_BAD         2e12 + 4
#define KMERSTAT_MARKED_FOR_GOODITER  2e12 + 5

#define MAX_SHORT 254

class Read;
class PositionRead;
class PositionKMer;
struct KMerStat;

typedef std::map<PositionKMer, KMerStat> KMerStatMap;
typedef std::pair<PositionKMer, KMerStat> KMerCount;
typedef std::pair<std::string, std::pair<uint32_t, double> > StringCount;

struct QualBitSet {
  unsigned char q_[K];

  QualBitSet(const unsigned char *q = NULL) {
    if (q)
      memcpy(q_, q, K);
  }

  // QualBitSet is POD-like object, use default stuff.
  QualBitSet(const QualBitSet &qbs) = default;
  QualBitSet& operator=(const QualBitSet &qbs) = default;
  // Workaround gcc bug, it cannot generate default move ctor...
  QualBitSet& operator=(QualBitSet &&qbs) {
    memcpy(q_, qbs.q_, K);

    return *this;
  }

  QualBitSet& operator+=(const unsigned char *data) {
    for (size_t i = 0; i < K; ++i)
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
    memcpy(q_, value, K);
  }
};

struct KMerStat {
  KMerStat(uint32_t cnt, hint_t cng, float kquality, const unsigned char *quality) : changeto(cng), totalQual(kquality), count(cnt), qual(quality) { }
  KMerStat() : changeto(KMERSTAT_BAD), totalQual(1.0), count(0), qual() { }

  hint_t changeto;
  float totalQual;
  uint32_t count;
  QualBitSet qual;

  bool isGood() const { return changeto >= KMERSTAT_GOOD; }
  bool isGoodForIterative() const { return (changeto == KMERSTAT_GOODITER); }
  void makeGoodForIterative() { changeto = KMERSTAT_GOODITER; }
  void markGoodForIterative() { changeto = KMERSTAT_MARKED_FOR_GOODITER; }
  bool isMarkedGoodForIterative() { return (changeto == KMERSTAT_MARKED_FOR_GOODITER); }
  bool change() const { return changeto < KMERSTAT_CHANGE; }
};

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
  os.write((char*)&k.changeto, sizeof(k.changeto));
  os.write((char*)&k.totalQual, sizeof(k.totalQual));
  return binary_write(os, k.qual);
}

template<class Reader>
inline void binary_read(Reader &is, KMerStat &k) {
  is.read((char*)&k.count, sizeof(k.count));
  is.read((char*)&k.changeto, sizeof(k.changeto));
  is.read((char*)&k.totalQual, sizeof(k.totalQual));
  binary_read(is, k.qual);
}

#endif //  HAMMER_KMERSTAT_HPP_
