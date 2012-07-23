//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef HAMMER_SUBKMER_SORTER_HPP
#define HAMMER_SUBKMER_SORTER_HPP

#include "kmer_stat.hpp"
#include "hammer_stats.hpp"
#include "mmapped_reader.hpp"
#include "position_kmer.hpp"

#include "sequence/seq.hpp"

#include <iostream>
#include <vector>

class unionFindClass;

struct SubKMer {
  size_t idx;
  Seq<K> data;
};

template<class Reader>
inline void binary_read(Reader &is, SubKMer &s) {
  Seq<K>::DataType seq_data[Seq<K>::DataSize];

  is.read((char*)&s.idx, sizeof(s.idx));
  is.read((char*)seq_data, sizeof(seq_data));

  s.data = Seq<K>(seq_data);
}

template<class Writer>
inline Writer &binary_write(Writer &os, const SubKMer &s) {
  Seq<K>::DataType seq_data[Seq<K>::DataSize];
  s.data.copy_data(seq_data);

  os.write((char*)&s.idx, sizeof(s.idx));
  os.write((char*)seq_data, sizeof(seq_data));

  return os;
}

static_assert(sizeof(SubKMer) == 16, "Too big SubKMer");


struct SubKMerDummySerializer {
  SubKMer serialize(hint_t idx, size_t fidx = -1ULL) const {
    SubKMer s;

    s.idx = (fidx == -1ULL ? idx : fidx);
    const char *seq = Globals::blob + idx;
    s.data = Seq<K>(seq, 0, K, /* raw */ true);

    // Yay for NRVO!
    return s;
  }
};

class SubKMerPartSerializer{
  size_t from_;
  size_t to_;

public:
  SubKMerPartSerializer(size_t from, size_t to)
      :from_(from), to_(to) { VERIFY(to_ - from_ <= K); }

  SubKMer serialize(hint_t idx, size_t fidx = -1ULL) const {
    SubKMer s;

    s.idx = (fidx == -1ULL ? idx : fidx);
    const char *seq = Globals::blob + idx + from_;
    s.data = Seq<K>(seq,
                    0, to_ - from_,
                    /* raw */ true);

    // Yay for NRVO!
    return s;
  }
};

class SubKMerStridedSerializer{
  size_t from_;
  size_t to_;
  size_t stride_;

public:
  SubKMerStridedSerializer(size_t from, size_t stride)
      :from_(from), stride_(stride) { VERIFY(from_ + stride_ <= K); }

  SubKMer serialize(hint_t idx, size_t fidx = -1ULL) const {
    SubKMer s;

    s.idx = (fidx == -1ULL ? idx : fidx);

    size_t sz = (K - from_ + stride_ - 1) / stride_;

    std::string str(sz, 'A');
    for (size_t i = from_, j = 0; i < K; i+= stride_, ++j)
      str[j] = Globals::blob[idx + i];

    s.data = Seq<K>(str, 0, sz);

    // Yay for NRVO!
    return s;
  }
};

class SubKMerBlockFile {
  MMappedReader ifs_;

 public:
  SubKMerBlockFile(const std::string &fname, bool unlink = false)
      : ifs_(fname, unlink) { }

  bool get_block(std::vector<size_t> &block) {
    block.clear();
#if 0
    block.shrink_to_fit();
#else
    std::vector<size_t>().swap(block);
#endif

    if (!ifs_.good())
      return false;

    size_t sz;
    ifs_.read((char*)&sz, sizeof(sz));
    block.resize(sz);
    for (size_t i = 0; i < sz; ++i) {
      SubKMer s;
      binary_read(ifs_, s);
      block[i] = s.idx;
    }

    return true;
  }
};

template<class Writer,
         class SubKMerSerializer = SubKMerDummySerializer>
void serialize(Writer &os,
               const std::vector<hint_t> &block, const std::vector<size_t> *fidx = NULL,
               const SubKMerSerializer &serializer = SubKMerSerializer()) {
  size_t sz = (fidx == NULL ? block.size() : fidx->size());
  os.write((char*)&sz, sizeof(sz));
  for (size_t i = 0, e = sz; i != e; ++i) {
    size_t idx = (fidx == NULL ? i : (*fidx)[i]);
    SubKMer s = serializer.serialize(block[idx], idx);
    binary_write(os, s);
  }
}

class SubKMerSplitter {
  const std::string ifname_;
  const std::string ofname_;

 public:
  SubKMerSplitter(const std::string &ifname, const std::string &ofname)
      : ifname_(ifname), ofname_(ofname) {}

  template<class Writer>
  void serialize(Writer &os,
                 const std::vector<SubKMer>::iterator &start,
                 const std::vector<SubKMer>::iterator &end) {
    size_t sz = end - start;

    os.write((char*)&sz, sizeof(sz));
    for (auto I = start, E = end; I != E; ++I)
      binary_write(os, *I);
  }

  template<class Reader>
  void deserialize(std::vector<SubKMer> &res,
                   Reader &is) {
    res.clear();
#if 0
    res.shrink_to_fit();
#else
    std::vector<SubKMer>().swap(res);
#endif

    size_t sz;
    is.read((char*)&sz, sizeof(sz));
    res.resize(sz);

    for (size_t i = 0, e = sz; i != e; ++i)
      binary_read(is, res[i]);
  }

  std::pair<size_t, size_t> split();
};

class KMerHamClusterer {
  unsigned tau_;

 public:
  KMerHamClusterer(unsigned tau)
      : tau_(tau) {}

  void cluster(const std::string &prefix,
               const std::vector<hint_t> &kmers,
               unionFindClass &uf);
};

#endif // HAMMER_SUBKMER_SORTER_HPP
