//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _HAMMER_KMERINDEX_HPP_
#define _HAMMER_KMERINDEX_HPP_

#include "mph_index.h"
#include "io/mmapped_reader.hpp"
#include "io/mmapped_writer.hpp"
#include "io/pointer_iterator.hpp"

#include "openmp_wrapper.h"

#include "logger/logger.hpp"
#include "path_helper.hpp"

#include <boost/lexical_cast.hpp>

#include <algorithm>
#ifdef USE_GLIBCXX_PARALLEL
#include <parallel/algorithm>
#endif
#include <vector>

#include <cmath>

template<class Seq>
class KMerIndexBuilder;

template<class Seq>
class KMerIndex {
 public:
  typedef typename Seq::hash hash_function;
  typedef typename MMappedRecordArrayReader<typename Seq::DataType>::iterator::value_type KMerRawData;
  typedef typename MMappedRecordArrayReader<typename Seq::DataType>::iterator::reference  KMerRawReference;

 private:
  // This is really the most fragile part of the whole story.  Basically, we're
  // building the PHM with "raw" data, but query the PHM with real Key-s!  We're
  // relying on fact that hashes are exactly the same in both cases (thus - two
  // almost equal implementations!).
  struct seeded_hash_function {
    cxxmph::h128 hash128(const KMerRawData &k, uint32_t seed) const {
      cxxmph::h128 h;
      MurmurHash3_x64_128(k.data(), k.data_size(), seed, &h);
      return h;
    }

    cxxmph::h128 hash128(const KMerRawReference k, uint32_t seed) const {
      cxxmph::h128 h;
      MurmurHash3_x64_128(k.data(), k.data_size(), seed, &h);
      return h;
    }

    cxxmph::h128 hash128(const Seq &k, uint32_t seed) const {
      cxxmph::h128 h;
      MurmurHash3_x64_128(k.data(), k.data_size() * sizeof(typename Seq::DataType), seed, &h);
      return h;
    }
  };

  typedef cxxmph::SimpleMPHIndex<Seq, seeded_hash_function> KMerDataIndex;

 public:
  KMerIndex(unsigned k):k_(k), index_(NULL),  num_buckets_(0) {}

  KMerIndex(const KMerIndex&) = delete;
  KMerIndex& operator=(const KMerIndex&) = delete;

  ~KMerIndex() { clear(); }

  void clear() {
    num_buckets_ = 0;
    bucket_starts_.clear();

    for (size_t i = 0; i < num_buckets_; ++i) {
      index_[i].clear();
    }

    delete[] index_;
    index_ = NULL;
  }

  size_t mem_size() {
    size_t sz = 0;
    for (size_t i = 0; i < num_buckets_; ++i)
      sz += index_[i].mem_size();

    return sz;
  }

  size_t seq_idx(const Seq &s) const {
    size_t bucket = seq_bucket(s);

    return bucket_starts_[bucket] + index_[bucket].index(s);
  }

  size_t raw_seq_idx(const KMerRawData &data) const {
    size_t bucket = raw_seq_bucket(data);

    return bucket_starts_[bucket] + index_[bucket].index(data);
  }

  template<class Writer>
  void serialize(Writer &os) const {
    os.write((char*)&num_buckets_, sizeof(num_buckets_));
    for (size_t i = 0; i < num_buckets_; ++i)
      index_[i].serialize(os);
    os.write((char*)&bucket_starts_[0], (num_buckets_ + 1) * sizeof(bucket_starts_[0]));
  }

  template<class Reader>
  void deserialize(Reader &is) {
    clear();

    is.read((char*)&num_buckets_, sizeof(num_buckets_));

    index_ = new KMerDataIndex[num_buckets_];
    for (size_t i = 0; i < num_buckets_; ++i)
      index_[i].deserialize(is);

    bucket_starts_.resize(num_buckets_ + 1);
    is.read((char*)&bucket_starts_[0], (num_buckets_ + 1) * sizeof(bucket_starts_[0]));
  }

 private:
  unsigned k_;
  KMerDataIndex *index_;

  size_t num_buckets_;
  std::vector<size_t> bucket_starts_;

  size_t seq_bucket(const Seq &s) const {
    return hash_function()(s) % num_buckets_;
  }
  size_t raw_seq_bucket(const KMerRawData &data) const {
    return hash_function()(data.data(), data.size()) % num_buckets_;
  }

  friend class KMerIndexBuilder<Seq>;
};

template<class Seq>
class KMerSplitter {
 public:
  typedef typename Seq::hash hash_function;

  KMerSplitter(const std::string &work_dir)
      : work_dir_(work_dir) {}

  virtual path::files_t Split(size_t num_files) = 0;

 protected:
  const std::string &work_dir_;
  hash_function hash_;

  std::string GetRawKMersFname(unsigned suffix) const {
    return path::append_path(work_dir_, "kmers.raw." + boost::lexical_cast<std::string>(suffix));
  }
  unsigned GetFileNumForSeq(const Seq &s, size_t total) const {
    return hash_(s) % total;
  }

  DECL_LOGGER("K-mer Splitting");
};

template<class Seq>
class KMerIndexBuilder {
  std::string work_dir_;
  unsigned num_buckets_;
  unsigned num_threads_;

 public:
  KMerIndexBuilder(const std::string &workdir,
                   unsigned num_buckets, unsigned num_threads)
      : work_dir_(workdir), num_buckets_(num_buckets), num_threads_(num_threads) {}
  size_t BuildIndex(KMerIndex<Seq> &out, KMerSplitter<Seq> &splitter,
                    bool save_final = false);

  unsigned num_buckets() const { return num_buckets_; }

  std::string GetFinalKMersFname() const {
    return path::append_path(work_dir_, "kmers.final");
  }

 private:
  size_t MergeKMers(const std::string &ifname, const std::string &ofname, unsigned K);
  std::string GetUniqueKMersFname(unsigned suffix) const {
    return path::append_path(work_dir_, "kmers.unique." + boost::lexical_cast<std::string>(suffix));
  }
  std::string GetMergedKMersFname(unsigned suffix) const {
    return path::append_path(work_dir_, "kmers.merged." + boost::lexical_cast<std::string>(suffix));
  }

  DECL_LOGGER("K-mer Index Building");
};

template<class Seq>
size_t KMerIndexBuilder<Seq>::MergeKMers(const std::string &ifname, const std::string &ofname,
                                         unsigned K) {
  MMappedRecordArrayReader<typename Seq::DataType> ins(ifname, Seq::GetDataSize(K), /* unlink */ true);

  // Sort the stuff
  std::sort(ins.begin(), ins.end(), array_less<typename Seq::DataType>());

  // FIXME: Use something like parallel version of unique_copy but with explicit
  // resizing.
  auto it = std::unique(ins.begin(), ins.end(), array_equal_to<typename Seq::DataType>());

  MMappedRecordArrayWriter<typename Seq::DataType> os(ofname, Seq::GetDataSize(K));
  os.resize(it - ins.begin());
  std::copy(ins.begin(), it, os.begin());

  return it - ins.begin();
}

template<class Seq>
size_t KMerIndexBuilder<Seq>::BuildIndex(KMerIndex<Seq> &index, KMerSplitter<Seq> &splitter,
                                         bool save_final) {
  index.clear();
  unsigned K = index.k_;

  INFO("Building kmer index ");

  // Split k-mers into buckets.
  path::files_t raw_kmers = splitter.Split(num_buckets_ * num_threads_);

  INFO("Starting k-mer counting.");
  size_t kmers = 0;
# pragma omp parallel for shared(raw_kmers) num_threads(num_threads_) schedule(dynamic)
  for (unsigned iFile = 0; iFile < raw_kmers.size(); ++iFile) {
    kmers += MergeKMers(raw_kmers[iFile], GetUniqueKMersFname(iFile), K);
  }
  INFO("K-mer counting done. There are " << kmers << " kmers in total. ");

  INFO("Merging temporary buckets.");
  for (unsigned i = 0; i < num_buckets_; ++i) {
    MMappedRecordArrayWriter<typename Seq::DataType> os(GetMergedKMersFname(i), Seq::GetDataSize(K));
    for (unsigned j = 0; j < num_threads_; ++j) {
      MMappedRecordArrayReader<typename Seq::DataType> ins(GetUniqueKMersFname(i + j * num_buckets_), Seq::GetDataSize(K), /* unlink */ true);
      size_t sz = ins.end() - ins.begin();
      os.reserve(sz);
      os.write(ins.data(), sz);
    }
  }

  index.num_buckets_ = num_buckets_;
  index.bucket_starts_.resize(num_buckets_ + 1);
  index.index_ = new typename KMerIndex<Seq>::KMerDataIndex[num_buckets_];

  INFO("Building perfect hash indices");
# pragma omp parallel for shared(index)
  for (unsigned iFile = 0; iFile < num_buckets_; ++iFile) {
    typename KMerIndex<Seq>::KMerDataIndex &data_index = index.index_[iFile];
    MMappedRecordArrayReader<typename Seq::DataType> ins(GetMergedKMersFname(iFile), Seq::GetDataSize(K), /* unlink */ !save_final);
    size_t sz = ins.end() - ins.begin();
    index.bucket_starts_[iFile + 1] = sz;
    if (!data_index.Reset(ins.begin(), ins.end(), sz)) {
      INFO("Something went really wrong (read = this should not happen). Try to restart and see if the problem will be fixed.");
      exit(-1);
    }
  }
  // Finally, record the sizes of buckets.
  for (unsigned iFile = 1; iFile < num_buckets_; ++iFile)
    index.bucket_starts_[iFile] += index.bucket_starts_[iFile - 1];

  if (save_final) {
    INFO("Merging final buckets.");
    MMappedRecordArrayWriter<typename Seq::DataType> os(GetFinalKMersFname(), Seq::GetDataSize(K));
    for (unsigned j = 0; j < num_buckets_; ++j) {
      MMappedRecordArrayReader<typename Seq::DataType> ins(GetMergedKMersFname(j), Seq::GetDataSize(K), /* unlink */ true);
      size_t sz = ins.end() - ins.begin();
      os.reserve(sz);
      os.write(ins.data(), sz);
    }
  }

  double bits_per_kmer = 8.0 * index.mem_size() / kmers;
  INFO("Index built. Total " << index.mem_size() << " bytes occupied (" << bits_per_kmer << " bits per kmer).");

  return kmers;
}

#endif // _HAMMER_KMERINDEX_HPP_
