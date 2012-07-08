//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef HAMMER_KMERSTAT_HPP_
#define HAMMER_KMERSTAT_HPP_

#include <stdint.h>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <bitset>

const uint32_t K = 21;
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
  std::vector<unsigned char> q;
  QualBitSet() : q(K) {
    // q.reset();
  }
  QualBitSet(size_t n) : q(n) {
    // q.reset();
  }
  unsigned short operator[](size_t n) const {
    return (unsigned short)q[n];
  }

  unsigned short at(size_t n) const {
    return (unsigned short)q.at(n);
  }

  void set(size_t n, unsigned short value) {
    q[n] = (unsigned char)value;
  }
};

struct KMerStat {
  KMerStat (bool first, uint32_t cnt, hint_t cng, double quality) : count(cnt), changeto(cng), totalQual(quality), qual(first ? 0 : K) { }
  // KMerStat (uint32_t cnt, hint_t cng, double quality) : count(cnt), changeto(cng), totalQual(quality), qual() { }
  KMerStat () : count(0), changeto(KMERSTAT_BAD), totalQual(1), qual(0) { }

  uint32_t count;
  hint_t changeto;
  double totalQual;
  QualBitSet qual;

  bool isGood() const { return changeto >= KMERSTAT_GOOD; }
  bool isGoodForIterative() const { return (changeto == KMERSTAT_GOODITER); }
  void makeGoodForIterative() { changeto = KMERSTAT_GOODITER; }
  void markGoodForIterative() { changeto = KMERSTAT_MARKED_FOR_GOODITER; }
  bool isMarkedGoodForIterative() { return (changeto == KMERSTAT_MARKED_FOR_GOODITER); }
  bool change() const { return changeto < KMERSTAT_CHANGE; }
};

inline std::ostream& binary_write(std::ostream &os, const QualBitSet &qbs) {
  size_t sz = qbs.q.size();

  os.write((char*)&sz, sizeof(sz));
  os.write((char*)&qbs.q[0], sz*sizeof(qbs.q[0]));

  return os;
}

inline void binary_read(std::istream &is, QualBitSet &qbs) {
  size_t sz;

  is.read((char*)&sz, sizeof(sz));
  qbs.q.resize(sz);
  is.read((char*)&qbs.q[0], sz*sizeof(qbs.q[0]));
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

char getQual(const KMerCount & kmc, int i);

#endif //  HAMMER_KMERSTAT_HPP_
