#pragma once
//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io/kmers/mmapped_reader.hpp"
#include "utils/filesystem/temporary.hpp"

namespace utils {

template<class Seq>
struct kmer_index_traits {
  typedef Seq SeqType;
  typedef fs::TmpFile ResultFile;
  typedef MMappedRecordArrayReader<typename Seq::DataType> RawKMerStorage;
  typedef typename RawKMerStorage::iterator             raw_data_iterator;
  typedef typename RawKMerStorage::const_iterator       raw_data_const_iterator;
  typedef typename RawKMerStorage::iterator::value_type KMerRawData;
  typedef typename RawKMerStorage::iterator::reference  KMerRawReference;
  typedef typename RawKMerStorage::const_iterator::reference  KMerRawConstReference;

  struct raw_equal_to {
    bool operator()(const Seq &lhs, const KMerRawReference rhs) {
      return (adt::array_equal_to<typename Seq::DataType>()(lhs.data(), lhs.data_size(), rhs));
    }
  };

  struct raw_create {
    Seq operator()(unsigned K, const KMerRawReference kmer) {
      return Seq(K, kmer.data());
    }
    Seq operator()(unsigned K, const KMerRawConstReference kmer) {
      return Seq(K, kmer.data());
    }
  };

  struct hash_function {
    uint64_t operator()(const Seq &k, uint64_t seed = 0) const{
        return typename Seq::hash()(k, (uint32_t)seed);
    }
    uint64_t operator()(const KMerRawReference k, uint64_t seed = 0) const {
        return typename Seq::hash()(k.data(), k.size(), (uint32_t)seed);
    }
  };

  template<class Writer>
  static void raw_serialize(Writer &writer, RawKMerStorage *data) {
    size_t sz = data->data_size(), elcnt = data->elcnt();
    unsigned PageSize = getpagesize();
    writer.write((char*)&sz, sizeof(sz));
    writer.write((char*)&elcnt, sizeof(elcnt));
    // Make sure data is aligned to the page boundary
    size_t cpos = writer.tellp();
    size_t pos = (cpos + PageSize - 1 + sizeof(size_t)) / PageSize * PageSize;
    size_t off = pos - writer.tellp();
    writer.write((char*)&off, sizeof(off));
    writer.seekp(pos);
    writer.write((char*)data->data(), data->data_size());
  }

  template<class Writer>
  static void raw_serialize(Writer &writer, const std::unique_ptr<RawKMerStorage> &data) {
    raw_serialize(writer, data.get());
  }

  template<class Reader>
  static std::unique_ptr<RawKMerStorage> raw_deserialize(Reader &reader, const std::string &FileName) {
    size_t sz, off, elcnt;
    reader.read((char*)&sz, sizeof(sz));
    reader.read((char*)&elcnt, sizeof(elcnt));
    reader.read((char*)&off, sizeof(off));
    off -= sizeof(off);
    off += reader.tellg();

    return std::unique_ptr<RawKMerStorage>(new RawKMerStorage(FileName, elcnt, false, off, sz));
  }

};
}
