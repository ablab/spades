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
    SubKMer s;

    // FIXME: Get rid of string here!
    std::string seq = k.str();
    return SubKMer(seq.data(),
                   from_, to_ - from_,
                   /* raw */ true);
  }
};

class SubKMerStridedSerializer{
  size_t from_;
  size_t stride_;

public:
  SubKMerStridedSerializer(size_t from, size_t stride)
      :from_(from), stride_(stride) { VERIFY(from_ + stride_ <= hammer::K); }

  SubKMer serialize(hammer::KMer k) const {
    SubKMer s;

    size_t sz = (hammer::K - from_ + stride_ - 1) / stride_;

    // FIXME: Get rid of strings here!
    std::string str(sz, 'A');
    for (size_t i = from_, j = 0; i < hammer::K; i+= stride_, ++j)
      str[j] = nucl(k[i]);

    return SubKMer(str, 0, sz);
  }
};

class BlockFile {
  MMappedReader ifs_;

 public:
  BlockFile(const std::string &fname, bool unlink = false)
      : ifs_(fname, unlink) { }

  bool read_block(std::vector<size_t> &block) {
    block.clear();
    block.shrink_to_fit();

    if (!ifs_.good())
      return false;

    size_t sz;
    ifs_.read((char*)&sz, sizeof(sz));
    block.resize(sz);
    ifs_.read((char*)block.data(), sz * sizeof(block[0]));

    return true;
  }
};

template<class Writer,
         class SubKMerSerializer>
void serialize(Writer &blocks, Writer &kmers,
               const KMerData &data, const std::vector<size_t> *block = NULL,
               const SubKMerSerializer &serializer = SubKMerSerializer()) {
  size_t sz = (block == NULL ? data.size() : block->size());
  blocks.write((char*)&sz, sizeof(sz));
  if (block) {
    blocks.write((char*)block->data(), sz * sizeof((*block)[0]));
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
  const std::string ofname_;

 public:
  SubKMerSplitter(const std::string &bifname, const std::string &kifname,
                  const std::string &ofname)
      : bifname_(bifname), kifname_(kifname), ofname_(ofname) {}

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
