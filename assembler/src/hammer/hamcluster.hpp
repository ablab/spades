//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef HAMMER_SUBKMER_SORTER_HPP
#define HAMMER_SUBKMER_SORTER_HPP

#include "kmer_stat.hpp"
#include "kmer_data.hpp"
#include "io/mmapped_reader.hpp"

#include "logger/logger.hpp"
#include "sequence/seq.hpp"

#include <iostream>
#include <vector>

class ConcurrentDSU;

struct SubKMer {
  size_t idx;
  Seq<hammer::K> data;
};

template<class Reader>
inline void binary_read(Reader &is, SubKMer &s) {
  Seq<hammer::K>::DataType seq_data[Seq<hammer::K>::DataSize];

  is.read((char*)&s.idx, sizeof(s.idx));
  is.read((char*)seq_data, sizeof(seq_data));

  s.data = Seq<hammer::K>(seq_data);
}

template<class Writer>
inline Writer &binary_write(Writer &os, const SubKMer &s) {
  Seq<hammer::K>::DataType seq_data[Seq<hammer::K>::DataSize];
  s.data.copy_data(seq_data);

  os.write((char*)&s.idx, sizeof(s.idx));
  os.write((char*)seq_data, sizeof(seq_data));

  return os;
}

static_assert(sizeof(SubKMer) == 16, "Too big SubKMer");


struct SubKMerDummySerializer {
  SubKMer serialize(hammer::KMer k, size_t fidx) const {
    SubKMer s;

    s.idx = fidx;
    s.data = k;

    // Yay for NRVO!
    return s;
  }
};

class SubKMerPartSerializer{
  size_t from_;
  size_t to_;

public:
  SubKMerPartSerializer(size_t from, size_t to)
      :from_(from), to_(to) { VERIFY(to_ - from_ <= hammer::K); }

  SubKMer serialize(hammer::KMer k, size_t fidx) const {
    SubKMer s;

    s.idx = fidx;
    // FIXME: Get rid of string here!
    std::string seq = k.str();
    s.data = Seq<hammer::K>(seq.data(),
                            from_, to_ - from_,
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
      :from_(from), stride_(stride) { VERIFY(from_ + stride_ <= hammer::K); }

  SubKMer serialize(hammer::KMer k, size_t fidx) const {
    SubKMer s;

    s.idx = fidx;

    size_t sz = (hammer::K - from_ + stride_ - 1) / stride_;

    // FIXME: Get rid of strings here!
    std::string str(sz, 'A');
    for (size_t i = from_, j = 0; i < hammer::K; i+= stride_, ++j)
      str[j] = nucl(k[i]);

    s.data = Seq<hammer::K>(str, 0, sz);

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
               const KMerData &data, const std::vector<size_t> *block = NULL,
               const SubKMerSerializer &serializer = SubKMerSerializer()) {
  size_t sz = (block == NULL ? data.size() : block->size());
  os.write((char*)&sz, sizeof(sz));
  for (size_t i = 0, e = sz; i != e; ++i) {
    size_t idx = (block == NULL ? i : (*block)[i]);
    SubKMer s = serializer.serialize(data[idx].kmer(), idx);
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

  void cluster(const std::string &prefix, const KMerData &data, ConcurrentDSU &uf);
 private:
  DECL_LOGGER("Hamming Clustering");
};

#endif // HAMMER_SUBKMER_SORTER_HPP
