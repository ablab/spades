//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef KMERNO_HPP_
#define KMERNO_HPP_

#include <vector>
#include <string>
#include <map>
#include <fstream>

#include <boost/unordered_map.hpp>
#include "google/sparse_hash_map"
#include "google/dense_hash_map"

//#define GOOGLE_SPARSE_MAP
#define BOOST_UNORDERED_MAP
//#define GOOGLE_DENSE_MAP

#include "kmer_stat.hpp"
#include "half.hpp"

const uint64_t KMERNO_HASH_MODULUS = 2305843009213693951;
const uint64_t KMERNO_HASH_Q = 3712758430079221;
const uint64_t KMERNO_HASH_Q_INV = 2250585152990002931;
const uint64_t KMERNO_HASH_Q_POW_K_MINUS_ONE = 412252044596125152;

class KMerNo {
public:
  explicit KMerNo(hint_t no = -1, float qual = 1.0) {
    index = no;
    errprob = prob_half::convert(qual);
  }

  bool equal(const KMerNo & kmerno) const;
  bool equal(const KMerCount & kmc) const;
  std::string str() const;
  bool less(const KMerNo &r) const;
  bool greater(const KMerNo &r) const;

  static uint64_t new_hash(hint_t index);
  static uint64_t next_hash(uint64_t old_hash, hint_t new_index);

  struct hash {
    uint64_t operator() (const KMerNo &kn) const;
  };

  struct string_hash {
    uint64_t operator() (const std::string &kn) const;
  };

  struct are_equal {
    bool operator() (const KMerNo &l, const KMerNo &r) const;
  };

  struct is_less {
    bool operator() (const KMerNo &l, const KMerNo &r) const;
  };

  struct is_less_kmercount {
    bool operator() (const KMerCount &l, const KMerCount &r) const;
  };

  hint_t getIndex() const { return index; }
  void setIndex(hint_t no) { index = no; }
  prob_half getQual() const { prob_half q; q.setBits(errprob); return q; }
  void setQual(float q) { errprob = prob_half::convert(q); }

private:
  hint_t index     : 48;
  uint16_t errprob : 16;
};

// FIXME: Eventually KMerNo should become POD-like class and thus can be possed by value
inline std::ostream& operator<<(std::ostream &os, const KMerNo &k) {
  os << k.getIndex() << '\t' << k.getQual() << '\n';
  return os;
}

inline std::istream& operator>>(std::istream &is, KMerNo &k) {
  std::string buf;
  hint_t pos; double prob;

  std::getline(is, buf);
  sscanf(buf.c_str(), "%llu\t%lf", &pos, &prob);
  k.setIndex(pos); k.setQual(prob);

  return is;
}

inline void binary_read(std::istream &is, KMerNo &k) {
  is.read((char*)&k, sizeof(k));
}

inline void binary_write(std::ostream &os, const KMerNo &k) {
  os.write((char*)&k, sizeof(k));
}

#ifdef GOOGLE_SPARSE_MAP
  typedef google::sparse_hash_map<KMerNo, KMerCount *, KMerNo::hash, KMerNo::are_equal> KMerNoHashMap;
#endif
#ifdef BOOST_UNORDERED_MAP
  typedef boost::unordered_map<KMerNo, KMerCount *, KMerNo::hash, KMerNo::are_equal> KMerNoHashMap;
#endif
#ifdef GOOGLE_DENSE_MAP
  typedef google::dense_hash_map<KMerNo, KMerCount *, KMerNo::hash, KMerNo::are_equal> KMerNoHashMap;
#endif


#endif
