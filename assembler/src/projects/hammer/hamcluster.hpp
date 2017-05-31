//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef HAMMER_SUBKMER_SORTER_HPP
#define HAMMER_SUBKMER_SORTER_HPP

#include "kmer_stat.hpp"
#include "kmer_data.hpp"
#include "io/kmers/mmapped_reader.hpp"

#include "utils/logger/logger.hpp"
#include "sequence/seq.hpp"

#include <iostream>
#include <vector>
#include <common/adt/concurrent_dsu.hpp>


typedef Seq<(hammer::K + 1) / 2, uint32_t> SubKMer;

template<class Reader>
inline void binary_read(Reader &is, SubKMer &s) {
  SubKMer::DataType seq_data[SubKMer::DataSize];

  is.read((char*)seq_data, sizeof(seq_data));

  s = SubKMer(seq_data);
}

template<class Writer>
inline Writer &binary_write(Writer &os, const SubKMer &s) {
  SubKMer::DataType seq_data[SubKMer::DataSize];
  s.copy_data(seq_data);

  os.write((char*)seq_data, sizeof(seq_data));

  return os;
}

static_assert(sizeof(SubKMer) == 4, "Too big SubKMer");

class SubKMerPartSerializer{
  size_t from_;
  size_t to_;

public:
  SubKMerPartSerializer(size_t from, size_t to)
      :from_(from), to_(to) { VERIFY(to_ - from_ <= hammer::K); }

  SubKMer serialize(hammer::KMer k) const {
    SubKMer res;
    for (size_t i = 0; i < to_ - from_; ++i)
      res.set(i, k[from_ + i]);

    return res;
  }
};

class SubKMerStridedSerializer{
  size_t from_;
  size_t stride_;

public:
  SubKMerStridedSerializer(size_t from, size_t stride)
      :from_(from), stride_(stride) { VERIFY(from_ + stride_ <= hammer::K); }

  SubKMer serialize(hammer::KMer k) const {
    SubKMer res;

    for (size_t i = from_, j = 0; i < hammer::K; i+= stride_, ++j)
      res.set(j, k[i]);

    return res;
  }
};

template<class Writer,
         class SubKMerSerializer>
void serialize(Writer &blocks, Writer &kmers,
               const KMerData &data,
               const std::vector<size_t>::iterator *block = NULL, size_t sz = 0,
               const SubKMerSerializer &serializer = SubKMerSerializer()) {
  if (sz == 0)
    sz = data.size();

  blocks.write((char*)&sz, sizeof(sz));
  if (block) {
    blocks.write((char*)&**block, sz * sizeof((*block)[0]));
  } else {
    for (size_t i = 0, e = sz; i != e; ++i)
      blocks.write((char*)&i, sizeof(i));
  }

  for (size_t i = 0, e = sz; i != e; ++i) {
    size_t idx = (block == NULL ? i : (*block)[i]);
    SubKMer s = serializer.serialize(data.kmer(idx));
    binary_write(kmers, s);
  }
}

class SubKMerSplitter {
  const std::string bifname_, kifname_;

 public:
  SubKMerSplitter(const std::string &bifname, const std::string &kifname)
      : bifname_(bifname), kifname_(kifname) {}

  template<class Writer>
  void serialize(Writer &os,
                 const std::vector<size_t>::iterator &start,
                 size_t sz) {
    os.write((char*)&sz, sizeof(sz));
    os.write((char*)&*start, sz * sizeof(*start));
  }

  template<class Reader>
  void deserialize(std::vector<size_t> &blocks,
                   std::vector<SubKMer> &kmers,
                   Reader &bis, Reader &kis) {
    kmers.clear(); blocks.clear();

    size_t sz;
    bis.read((char*)&sz, sizeof(sz));
    blocks.resize(sz);
    bis.read((char*)blocks.data(), sz * sizeof(blocks[0]));

    kmers.resize(sz);
    for (size_t i = 0, e = sz; i != e; ++i)
      binary_read(kis, kmers[i]);
  }

  template<class Op>
  std::pair<size_t, size_t> split(Op &&op);
};

class KMerHamClusterer {
  unsigned tau_;

 public:
  KMerHamClusterer(unsigned tau)
      : tau_(tau) {}

  void cluster(const std::string &prefix, const KMerData &data, dsu::ConcurrentDSU &uf);
 private:
  DECL_LOGGER("Hamming Clustering");
};

class TauOneKMerHamClusterer {
 public:
  TauOneKMerHamClusterer() {} 
  void cluster(const std::string &prefix, const KMerData &data, dsu::ConcurrentDSU &uf);
 private:
  DECL_LOGGER("tau = 1 Hamming Clustering");
};


#endif // HAMMER_SUBKMER_SORTER_HPP
