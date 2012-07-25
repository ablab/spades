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

#include <unordered_map>

#include "kmer_stat.hpp"
#include "half.hpp"

class KMerNo {
public:
  explicit KMerNo(hint_t no = -1, float qual = 1.0, const Seq<K> s = Seq<K>()) {
    info.index = no;
    info.errprob = prob_half::convert(qual);
    seq = s;
  }

  bool operator==(const KMerNo &k) const {
    return seq == k.seq;
  }

  bool operator!=(const KMerNo &k) const {
    return seq != k.seq;
  }

  bool operator==(const KMerCount &kmc) const;

  bool less(const KMerNo &r) const {
    return Seq<K>::less2()(seq, r.seq);
  }

  std::string str() const;

  struct hash {
    uint64_t operator() (const KMerNo &kn) const;
  };

  struct are_equal {
    bool operator()(const KMerNo &l, const KMerNo &r) const {
      return l.seq == r.seq;
    }
  };

  struct is_less {
    bool operator()(const KMerNo &l, const KMerNo &r) const {
      return Seq<K>::less2()(l.seq, r.seq);
    }
  };

  hint_t getIndex() const { return info.index; }
  void setIndex(hint_t no) { info.index = no; }
  prob_half getQual() const { prob_half q; q.setBits(info.errprob); return q; }
  void setQual(float q) { info.errprob = prob_half::convert(q); }
  Seq<K> getSeq() const { return seq; }
  void setSeq(const Seq<K> s) { seq = s; }
private:
  struct {
    hint_t   index   : 48;
    uint16_t errprob : 16;
  } info;
  Seq<K>   seq;

  template<class Reader>
  friend void binary_read(Reader &is, KMerNo &s);
  template<class Writer>
  friend void binary_write(Writer &os, const KMerNo &s);
};

static_assert(sizeof(KMerNo) == 16, "Invalid size of KMerNo");

// FIXME: Unify with SubKMer

template<class Reader>
inline void binary_read(Reader &is, KMerNo &k) {
  Seq<K>::DataType seq_data[Seq<K>::DataSize];

  is.read((char*)&k.info, sizeof(k.info));
  is.read((char*)seq_data, sizeof(seq_data));

  k.seq = Seq<K>(seq_data);
}

template<class Writer>
inline Writer &binary_write(Writer &os, const KMerNo &k) {
  Seq<K>::DataType seq_data[Seq<K>::DataSize];
  k.seq.copy_data(seq_data);

  os.write((char*)&k.info, sizeof(k.info));
  os.write((char*)seq_data, sizeof(seq_data));

  return os;
}

typedef std::unordered_map<KMerNo, KMerCount *, KMerNo::hash, KMerNo::are_equal> KMerNoHashMap;

#endif
