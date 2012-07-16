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
  unsigned char* q;

  QualBitSet(bool empty = false):q(NULL) {
    if (!empty) {
      q = new unsigned char[K];
      memset(q, 0, K);
    }
  }
  ~QualBitSet() {
    delete[] q;
  }

  QualBitSet(const QualBitSet &qbs) {
    q = NULL;

    // Fill, if necessary.
    if (qbs.q) {
      q = new unsigned char[K];
      memcpy(q, qbs.q, K * sizeof(q[0]));
    }
  }
  QualBitSet& operator=(const QualBitSet &qbs) {
    if (&qbs != this) {
      delete[] q;
      q = NULL;

      // Fill, if necessary.
      if (qbs.q) {
        q = new unsigned char[K];
        memcpy(q, qbs.q, K * sizeof(q[0]));
      }
    }

    return *this;
  }

  unsigned short operator[](size_t n) const {
    return (unsigned short)q[n];
  }

  unsigned short at(size_t n) const {
    // FIXME: Bound checking
    return (unsigned short)q[n];
  }

  void set(size_t n, unsigned short value) {
    VERIFY(q != NULL);
    q[n] = (unsigned char)value;
  }
};

struct KMerStat {
  KMerStat (bool first, uint32_t cnt, hint_t cng, float quality) : changeto(cng), qual(first /* empty */), totalQual(quality), count(cnt) { }
  KMerStat () : changeto(KMERSTAT_BAD), qual(), totalQual(1.0), count(0) { }

  hint_t changeto;
  QualBitSet qual;
  float totalQual;
  uint32_t count;

  bool isGood() const { return changeto >= KMERSTAT_GOOD; }
  bool isGoodForIterative() const { return (changeto == KMERSTAT_GOODITER); }
  void makeGoodForIterative() { changeto = KMERSTAT_GOODITER; }
  void markGoodForIterative() { changeto = KMERSTAT_MARKED_FOR_GOODITER; }
  bool isMarkedGoodForIterative() { return (changeto == KMERSTAT_MARKED_FOR_GOODITER); }
  bool change() const { return changeto < KMERSTAT_CHANGE; }
};

inline std::ostream& binary_write(std::ostream &os, const QualBitSet &qbs) {
  size_t sz = (qbs.q ? K : 0);

  os.write((char*)&sz, sizeof(sz));
  if (sz)
    os.write((char*)&qbs.q[0], sz*sizeof(qbs.q[0]));

  return os;
}

inline void binary_read(std::istream &is, QualBitSet &qbs) {
  size_t sz;

  is.read((char*)&sz, sizeof(sz));
  delete[] qbs.q;
  qbs.q = NULL;
  
  if (sz) {
    VERIFY(sz == K);
    qbs.q = new unsigned char[sz];
    is.read((char*)&qbs.q[0], sz*sizeof(qbs.q[0]));
  }
}

inline std::ostream& binary_write(std::ostream &os, const KMerStat &k) {
  os.write((char*)&k.count, sizeof(k.count));
  os.write((char*)&k.changeto, sizeof(k.changeto));
  os.write((char*)&k.totalQual, sizeof(k.totalQual));
  return binary_write(os, k.qual);
}

inline void binary_read(std::istream &is, KMerStat &k) {
  is.read((char*)&k.count, sizeof(k.count));
  is.read((char*)&k.changeto, sizeof(k.changeto));
  is.read((char*)&k.totalQual, sizeof(k.totalQual));
  binary_read(is, k.qual);
}

#endif //  HAMMER_KMERSTAT_HPP_
