//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
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

typedef Seq<(hammer::K + 1) / 2, uint32_t> SubKMer;

struct SubKMerData {
  uint32_t idx;
  SubKMer data;
};

template<class Reader>
inline void binary_read(Reader &is, SubKMerData &s) {
  SubKMer::DataType seq_data[SubKMer::DataSize];

  is.read((char*)&s.idx, sizeof(s.idx));
  is.read((char*)seq_data, sizeof(seq_data));

  s.data = SubKMer(seq_data);
}

template<class Writer>
inline Writer &binary_write(Writer &os, const SubKMerData &s) {
  SubKMer::DataType seq_data[SubKMer::DataSize];
  s.data.copy_data(seq_data);

  os.write((char*)&s.idx, sizeof(s.idx));
  os.write((char*)seq_data, sizeof(seq_data));

  return os;
}

static_assert(sizeof(SubKMerData) == 8, "Too big SubKMer");

class SubKMerPartSerializer{
  size_t from_;
  size_t to_;

public:
  SubKMerPartSerializer(size_t from, size_t to)
      :from_(from), to_(to) { VERIFY(to_ - from_ <= hammer::K); }

  SubKMerData serialize(hammer::KMer k, size_t fidx) const {
    SubKMerData s;

    s.idx = (uint32_t)fidx;
    // FIXME: Get rid of string here!
    std::string seq = k.str();
    s.data = SubKMer(seq.data(),
                     from_, to_ - from_,
                     /* raw */ true);

    // Yay for NRVO!
    return s;
  }
};

class SubKMerStridedSerializer{
  size_t from_;
  size_t stride_;

public:
  SubKMerStridedSerializer(size_t from, size_t stride)
      :from_(from), stride_(stride) { VERIFY(from_ + stride_ <= hammer::K); }

  SubKMerData serialize(hammer::KMer k, size_t fidx) const {
    SubKMerData s;

    s.idx = (uint32_t)fidx;

    size_t sz = (hammer::K - from_ + stride_ - 1) / stride_;

    // FIXME: Get rid of strings here!
    std::string str(sz, 'A');
    for (size_t i = from_, j = 0; i < hammer::K; i+= stride_, ++j)
      str[j] = nucl(k[i]);

    s.data = SubKMer(str, 0, sz);

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
      SubKMerData s;
      binary_read(ifs_, s);
      block[i] = s.idx;
    }

    return true;
  }
};

template<class Writer,
         class SubKMerSerializer>
void serialize(Writer &os,
               const KMerData &data, const std::vector<size_t> *block = NULL,
               const SubKMerSerializer &serializer = SubKMerSerializer()) {
  size_t sz = (block == NULL ? data.size() : block->size());
  os.write((char*)&sz, sizeof(sz));
  for (size_t i = 0, e = sz; i != e; ++i) {
    size_t idx = (block == NULL ? i : (*block)[i]);
    SubKMerData s = serializer.serialize(data.kmer(idx), idx);
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
                 const std::vector<SubKMerData>::iterator &start,
                 const std::vector<SubKMerData>::iterator &end) {
    size_t sz = end - start;

    os.write((char*)&sz, sizeof(sz));
    for (auto I = start, E = end; I != E; ++I)
      binary_write(os, *I);
  }

  template<class Reader>
  void deserialize(std::vector<SubKMerData> &res,
                   Reader &is) {
    res.clear();
#if 0
    res.shrink_to_fit();
#else
    std::vector<SubKMerData>().swap(res);
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
