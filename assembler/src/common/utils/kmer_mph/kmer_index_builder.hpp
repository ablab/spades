#pragma once
//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "kmer_index.hpp"

#include "io/kmers/mmapped_reader.hpp"
#include "io/kmers/mmapped_writer.hpp"
#include "common/adt/kmer_vector.hpp"

#include "utils/parallel/openmp_wrapper.h"

#include "utils/logger/logger.hpp"
#include "utils/filesystem/path_helper.hpp"

#include "utils/memory_limit.hpp"
#include "utils/filesystem/file_limit.hpp"

#include "adt/iterator_range.hpp"
#include "adt/loser_tree.hpp"

#include "boomphf/BooPHF.h"

#include <libcxx/sort.hpp>

#include <algorithm>
#ifdef USE_GLIBCXX_PARALLEL
#include <parallel/algorithm>
#endif

#include "config.hpp"

#ifdef SPADES_USE_JEMALLOC
# include <jemalloc/jemalloc.h>
#endif

#include <fstream>
#include <vector>
#include <cmath>

#include "kmer_splitters.hpp"

namespace utils {


template<class Seq, class traits = kmer_index_traits<Seq> >
class KMerCounter {
public:
  typedef typename traits::raw_data_iterator       iterator;
  typedef typename traits::raw_data_const_iterator const_iterator;
  typedef typename traits::ResultFile              ResultFile;
  typedef typename traits::RawKMerStorage          RawKMerStorage;

  KMerCounter()
      : kmers_(0), counted_(false) {}

  virtual size_t kmer_size() const = 0;

  virtual size_t Count(unsigned num_buckets, unsigned num_threads) = 0;
  virtual size_t CountAll(unsigned num_buckets, unsigned num_threads, bool merge = true) = 0;
  virtual void MergeBuckets() = 0;

  virtual std::unique_ptr<RawKMerStorage> GetBucket(size_t idx, bool unlink = true) = 0;

  virtual ~KMerCounter() {}

  unsigned num_buckets() const { return num_buckets_; }
  size_t kmers() const { return kmers_; }
  bool counted() const { return counted_; }

 protected:
  unsigned num_buckets_;
  size_t kmers_;
  bool counted_;

  DECL_LOGGER("K-mer Counting");
};

template<class Seq, class traits = kmer_index_traits<Seq> >
class KMerDiskCounter : public KMerCounter<Seq> {
  typedef KMerCounter<Seq, traits> __super;
  typedef typename traits::RawKMerStorage BucketStorage;
  typedef typename traits::ResultFile ResultFile;
public:
  KMerDiskCounter(fs::TmpDir work_dir,
                  KMerSplitter<Seq> &splitter)
      : work_dir_(work_dir), splitter_(splitter), k_(splitter.K()) {
    kmer_prefix_ = work_dir_->tmp_file("kmers");
  }

  KMerDiskCounter(const std::string &work_dir,
                  KMerSplitter<Seq> &splitter)
      : KMerDiskCounter(fs::tmp::make_temp_dir(work_dir, "kmer_counter"), splitter) {}

  ~KMerDiskCounter() {}

  unsigned k() const { return k_; }

  size_t kmer_size() const override {
    return Seq::GetDataSize(k_) * sizeof(typename Seq::DataType);
  }

  std::unique_ptr<BucketStorage> GetBucket(size_t idx, bool unlink = true) override {
    VERIFY_MSG(this->counted_, "k-mers were not counted yet");
    return std::unique_ptr<BucketStorage>(new BucketStorage(GetMergedKMersFname((unsigned)idx), Seq::GetDataSize(k_), unlink));
  }

  size_t Count(unsigned num_buckets, unsigned num_threads) override {
    this->num_buckets_ = num_buckets;
    unsigned num_files = num_buckets * num_threads;

    // Split k-mers into buckets.
    INFO("Splitting kmer instances into " << num_files << " files using " << num_threads << " threads. This might take a while.");
    auto raw_kmers = splitter_.Split(num_files, num_threads);

    INFO("Starting k-mer counting.");
    size_t kmers = 0;
#   pragma omp parallel for shared(raw_kmers) num_threads(num_threads) schedule(dynamic) reduction(+:kmers)
    for (unsigned i = 0; i < raw_kmers.size(); ++i) {
      kmers += MergeKMers(*raw_kmers[i], GetUniqueKMersFname(i));
      raw_kmers[i].reset();
    }
    INFO("K-mer counting done. There are " << kmers << " kmers in total. ");
    if (!kmers) {
      FATAL_ERROR("No kmers were extracted from reads. Check the read lengths and k-mer length settings");
      exit(-1);
    }

    INFO("Merging temporary buckets.");
    for (unsigned i = 0; i < num_buckets; ++i) {
      std::string ofname = GetMergedKMersFname(i);
      std::ofstream ofs(ofname.c_str(), std::ios::out | std::ios::binary);
      for (unsigned j = 0; j < num_threads; ++j) {
        BucketStorage ins(GetUniqueKMersFname(i + j * num_buckets), Seq::GetDataSize(k_), /* unlink */ true);
        ofs.write((const char*)ins.data(), ins.data_size());
      }
    }

    this->kmers_ = kmers;
    this->counted_ = true;

    return kmers;
  }

  void MergeBuckets() override {
    INFO("Merging final buckets.");

    final_kmers_ = work_dir_->tmp_file("final_kmers");
    std::ofstream ofs(*final_kmers_, std::ios::out | std::ios::binary);
    for (unsigned j = 0; j < this->num_buckets_; ++j) {
      auto bucket = GetBucket(j, /* unlink */ true);
      ofs.write((const char*)bucket->data(), bucket->data_size());
    }
    ofs.close();
  }

  size_t CountAll(unsigned num_buckets, unsigned num_threads, bool merge = true) override {
    size_t kmers = Count(num_buckets, num_threads);
    if (merge)
      MergeBuckets();

    return kmers;
  }

  std::string GetMergedKMersFname(unsigned suffix) const {
    return kmer_prefix_->file() + ".merged." + std::to_string(suffix);
  }

  ResultFile final_kmers_file() {
    VERIFY_MSG(this->final_kmers_, "k-mers were not counted yet");
    return final_kmers_;
  }

private:
  fs::TmpDir work_dir_;
  fs::TmpFile kmer_prefix_;
  fs::TmpFile final_kmers_;
  KMerSplitter<Seq> &splitter_;
  unsigned k_;

  std::string GetUniqueKMersFname(unsigned suffix) const {
    return kmer_prefix_->file() + ".unique." + std::to_string(suffix);
  }

  size_t MergeKMers(const std::string &ifname, const std::string &ofname) {
    MMappedRecordArrayReader<typename Seq::DataType> ins(ifname, Seq::GetDataSize(k_), /* unlink */ true);

    std::string IdxFileName = ifname + ".idx";
    if (FILE *f = fopen(IdxFileName.c_str(), "rb")) {
      fclose(f);
      MMappedRecordReader<size_t> index(ifname + ".idx", true, -1ULL);

      // INFO("Total runs: " << index.size());

      // Prepare runs
      std::vector<adt::iterator_range<decltype(ins.begin())>> ranges;
      auto beg = ins.begin();
      for (size_t sz : index) {
        auto end = std::next(beg, sz);
        ranges.push_back(adt::make_range(beg, end));
        VERIFY(std::is_sorted(beg, end, adt::array_less<typename Seq::DataType>()));
        beg = end;
      }

      // Construct tree on top entries of runs
      adt::loser_tree<decltype(beg),
              adt::array_less<typename Seq::DataType>> tree(ranges);

      if (tree.empty()) {
        FILE *g = fopen(ofname.c_str(), "ab");
        if (!g)
          FATAL_ERROR("Cannot open temporary file " << ofname << " for writing");
        fclose(g);
        return 0;
      }

      // Write it down!
      adt::KMerVector<Seq> buf(k_, 1024*1024);
      auto pval = tree.pop();
      size_t total = 0;
      while (!tree.empty()) {
          buf.clear();
          for (size_t cnt = 0; cnt < buf.capacity() && !tree.empty(); ) {
              auto cval = tree.pop();
              if (!adt::array_equal_to<typename Seq::DataType>()(pval, cval)) {
                  buf.push_back(pval);
                  pval = cval;
                  cnt += 1;
              }
          }
          total += buf.size();

          FILE *g = fopen(ofname.c_str(), "ab");
          if (!g)
            FATAL_ERROR("Cannot open temporary file " << ofname << " for writing");
          size_t res = fwrite(buf.data(), buf.el_data_size(), buf.size(), g);
          if (res != buf.size())
            FATAL_ERROR("I/O error! Incomplete write! Reason: " << strerror(errno) << ". Error code: " << errno);
          fclose(g);
      }

      // Handle very last value
      {
        FILE *g = fopen(ofname.c_str(), "ab");
        if (!g)
          FATAL_ERROR("Cannot open temporary file " << ofname << " for writing");
        size_t res = fwrite(pval.data(), pval.data_size(), 1, g);
        if (res != 1)
          FATAL_ERROR("I/O error! Incomplete write! Reason: " << strerror(errno) << ". Error code: " << errno);
        fclose(g);
        total += 1;
      }

      return total;
    } else {
      // Sort the stuff
      libcxx::sort(ins.begin(), ins.end(), adt::array_less<typename Seq::DataType>());

      // FIXME: Use something like parallel version of unique_copy but with explicit
      // resizing.
      auto it = std::unique(ins.begin(), ins.end(), adt::array_equal_to<typename Seq::DataType>());

      MMappedRecordArrayWriter<typename Seq::DataType> os(ofname, Seq::GetDataSize(k_));
      os.resize(it - ins.begin());
      std::copy(ins.begin(), it, os.begin());

      return it - ins.begin();
    }
  }
};

template<class Index>
class KMerIndexBuilder {
  typedef typename Index::KMerSeq Seq;
  typedef typename Index::kmer_index_traits kmer_index_traits;

  unsigned num_buckets_;
  unsigned num_threads_;

 public:
  KMerIndexBuilder(unsigned num_buckets, unsigned num_threads)
      : num_buckets_(num_buckets), num_threads_(num_threads) {}
  size_t BuildIndex(Index &out, KMerCounter<Seq> &counter,
                    bool save_final = false);

 private:

  DECL_LOGGER("K-mer Index Building");
};

// FIXME: Should be able to build index just from the set of files
template<class Index>
size_t KMerIndexBuilder<Index>::BuildIndex(Index &index, KMerCounter<Seq> &counter,
                                           bool save_final) {
  index.clear();

  INFO("Building kmer index ");

  // First, count the unique k-mers
  if (!counter.counted())
    counter.Count(num_buckets_, num_threads_);

  size_t kmers = counter.kmers();
  unsigned buckets = counter.num_buckets();

  index.num_buckets_ = buckets;
  index.bucket_starts_.resize(buckets + 1);
  index.index_ = new typename KMerIndex<kmer_index_traits>::KMerDataIndex[buckets];

  INFO("Building perfect hash indices");

# pragma omp parallel for shared(index) num_threads(num_threads_)
  for (unsigned i = 0; i < buckets; ++i) {
    typename KMerIndex<kmer_index_traits>::KMerDataIndex &data_index = index.index_[i];
    auto bucket = counter.GetBucket(i, !save_final);
    size_t sz = bucket->end() - bucket->begin();
    index.bucket_starts_[i + 1] = sz;

    data_index = typename Index::KMerDataIndex(sz,
                                               boomphf::range(bucket->begin(), bucket->end()),
                                               1, 2.0, false, false);
  }

  // Finally, record the sizes of buckets.
  for (unsigned i = 1; i < buckets; ++i)
    index.bucket_starts_[i] += index.bucket_starts_[i - 1];

  if (save_final)
    counter.MergeBuckets();

  double bits_per_kmer = 8.0 * (double)index.mem_size() / (double)kmers;
  INFO("Index built. Total " << index.mem_size() << " bytes occupied (" << bits_per_kmer << " bits per kmer).");
  index.count_size();
  return kmers;
}
}
