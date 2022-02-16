#pragma once
//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "kmer_splitter.hpp"
#include "kmer_index.hpp"

#include "io/kmers/mmapped_reader.hpp"
#include "io/kmers/mmapped_writer.hpp"

#include "utils/parallel/openmp_wrapper.h"
#include "utils/memory_limit.hpp"
#include "utils/logger/logger.hpp"
#include "utils/filesystem/path_helper.hpp"
#include "utils/filesystem/file_limit.hpp"
#include "utils/perf/timetracer.hpp"

#include "adt/kmer_vector.hpp"
#include "adt/iterator_range.hpp"
#include "adt/loser_tree.hpp"

#include <boomphf/BooPHF.h>

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

namespace kmers {

// TODO: Make interface
template<class S, class traits = kmer_index_traits<S>>
class KMerDiskStorage {
  typedef typename traits::RawKMerStorage BucketStorage;
 public:
  typedef S Seq;
  typedef typename std::vector<fs::DependentTmpFile> Buckets;
  typedef typename kmer::KMerSegmentPolicy<Seq>       KMerSegmentPolicy;
  typedef typename std::pair<const typename Seq::DataType*, size_t> KMerRawData;

  class kmer_iterator :
      public boost::iterator_facade<kmer_iterator,
                                    KMerRawData,
                                    std::input_iterator_tag,
                                    KMerRawData> {
   public:
    // Default ctor, used to implement "end" iterator
    kmer_iterator()
        : inner_iterator_(),
          k_(0), kmer_bytes_(0) { }

    kmer_iterator(const std::string &FileName, unsigned k)
        : inner_iterator_(FileName, Seq::GetDataSize(k)),
          k_(k), kmer_bytes_(Seq::GetDataSize(k_) * sizeof(typename Seq::DataType)) {}

    void operator+=(size_t n) {
      inner_iterator_ += n;
    }

   private:
    friend class boost::iterator_core_access;

    void increment() {
      ++inner_iterator_;
    }

    bool equal(const kmer_iterator &other) const {
        return inner_iterator_ == other.inner_iterator_;
    }

    KMerRawData dereference() const {
      return { *inner_iterator_, kmer_bytes_ };
    }

    MMappedFileRecordArrayIterator<typename Seq::DataType> inner_iterator_;
    unsigned k_;
    size_t kmer_bytes_;
  };

  static_assert(std::is_nothrow_move_constructible<kmer_iterator>::value, "kmer_iterator must be nonthrow move constructible");

  KMerDiskStorage() {}
  
  KMerDiskStorage(fs::TmpDir work_dir, unsigned k,
                  KMerSegmentPolicy policy)
      : work_dir_(work_dir), k_(k), segment_policy_(std::move(policy)) {
    kmer_prefix_ = work_dir_->tmp_file("kmers");
    resize(policy.num_segments());
  }

  KMerDiskStorage(KMerDiskStorage &&) = default;
  KMerDiskStorage &operator=(KMerDiskStorage &&) = default;
  
  fs::DependentTmpFile create() {
    fs::DependentTmpFile res;
#pragma omp critical
    {
      buckets_.emplace_back(kmer_prefix_->CreateDep(std::to_string(buckets_.size())));
      res = buckets_.back();
    }
    return res;
  }

  fs::DependentTmpFile create(size_t idx) {
    fs::DependentTmpFile res = kmer_prefix_->CreateDep(std::to_string(idx));
    buckets_.at(idx) = res;
    return res;
  }

  void resize(size_t n) {
    buckets_.resize(n);
  }

  unsigned k() const { return k_; }

  size_t total_kmers() const {
    size_t fsize = 0;
    if (all_kmers_) {
      fsize = fs::filesize(*all_kmers_);
    } else {
      for (const auto &file : buckets_)
        fsize += fs::filesize(*file);
    }

    return fsize / (Seq::GetDataSize(k_) * sizeof(typename Seq::DataType));
  }

  fs::TmpFile final_kmers() {
    VERIFY_MSG(all_kmers_, "k-mers were not merged yet");
    return all_kmers_;
  }

  size_t bucket_size(size_t i) const {
    return fs::filesize(*buckets_.at(i)) / (Seq::GetDataSize(k_) * sizeof(typename Seq::DataType));
  }

  kmer_iterator bucket_begin(size_t i) const {
    return kmer_iterator(*buckets_.at(i), k_);
  }

  kmer_iterator bucket_end(size_t) const {
    return kmer_iterator();
  }

  auto bucket(size_t i) const {
    return adt::make_range(bucket_begin(i), bucket_end(i));
  }

  size_t num_buckets() const { return buckets_.size(); }
  KMerSegmentPolicy segment_policy() const { return segment_policy_; }

  void merge() {
    INFO("Merging final buckets.");
    TIME_TRACE_SCOPE("KMerDiskStorage::MergeFinal");

    all_kmers_ = work_dir_->tmp_file("final_kmers");
    std::ofstream ofs(*all_kmers_, std::ios::out | std::ios::binary);
    for (auto &entry : buckets_) {
      BucketStorage bucket(*entry, Seq::GetDataSize(k_), false);
      ofs.write((const char*)bucket.data(), bucket.data_size());
      entry.reset();
    }
    buckets_.clear();
    ofs.close();
  }


 private:
  fs::TmpDir work_dir_;
  fs::TmpFile kmer_prefix_;
  fs::TmpFile all_kmers_;
  unsigned k_;
  Buckets buckets_;
  KMerSegmentPolicy segment_policy_;
};


template<class S, class traits = kmer_index_traits<S> >
class KMerCounter {
public:
  typedef typename traits::raw_data_iterator       iterator;
  typedef typename traits::raw_data_const_iterator const_iterator;
  typedef typename traits::RawKMerStorage          RawKMerStorage;
  typedef S Seq;

  KMerCounter(unsigned k)
      : k_(k) { }

  virtual size_t kmer_size() const = 0;
  unsigned k() const { return k_; }

  virtual KMerDiskStorage<Seq> Count(unsigned num_buckets, unsigned num_threads) = 0;
  virtual KMerDiskStorage<Seq> CountAll(unsigned num_buckets, unsigned num_threads, bool merge = true) = 0;

  virtual ~KMerCounter() {}

 protected:
  unsigned k_;

  DECL_LOGGER("K-mer Counting");
};

template<class Seq, class traits = kmer_index_traits<Seq> >
class KMerDiskCounter : public KMerCounter<Seq> {
  typedef KMerCounter<Seq, traits> __super;
  typedef typename traits::RawKMerStorage BucketStorage;
  typedef typename traits::ResultFile ResultFile;
public:
  template<class Splitter>
  KMerDiskCounter(fs::TmpDir work_dir,
                  Splitter splitter)
      : __super(splitter.K()), splitter_(new Splitter{std::move(splitter)}), work_dir_(work_dir) {}

  template<class Splitter>
  KMerDiskCounter(const std::string &work_dir,
                  Splitter splitter)
      : KMerDiskCounter(fs::tmp::make_temp_dir(work_dir, "kmer_counter"), std::move(splitter)) {}

  ~KMerDiskCounter() {}

  size_t kmer_size() const override {
    return Seq::GetDataSize(this->k()) * sizeof(typename Seq::DataType);
  }

  KMerDiskStorage<Seq> Count(unsigned num_buckets, unsigned num_threads) override {
    // Split k-mers into buckets.
    INFO("Splitting kmer instances into " << num_buckets << " files using " << num_threads << " threads. This might take a while.");
    TIME_TRACE_BEGIN("KMerDiskCounter::Split");
    auto raw_kmers = splitter_->Split(num_buckets, num_threads);
    VERIFY(raw_kmers.size() == num_buckets);
    TIME_TRACE_END;

    INFO("Starting k-mer counting.");
    KMerDiskStorage<Seq> res(work_dir_, this->k(), splitter_->bucket_policy());
    size_t kmers = 0;
    {
        TIME_TRACE_SCOPE("KMerDiskCounter::Count");
#       pragma omp parallel for shared(raw_kmers) num_threads(num_threads) schedule(dynamic) reduction(+:kmers)
        for (size_t i = 0; i < raw_kmers.size(); ++i) {
          kmers += MergeKMers(*raw_kmers[i], *res.create(i));
          raw_kmers[i].reset();
        }
    }
    INFO("K-mer counting done. There are " << kmers << " kmers in total. ");
    if (!kmers) {
      FATAL_ERROR("No kmers were extracted from reads. Check the read lengths and k-mer length settings");
      exit(-1);
    }

    return res;
  }

  KMerDiskStorage<Seq> CountAll(unsigned num_buckets, unsigned num_threads, bool merge = true) override {
    auto storage = Count(num_buckets, num_threads);
    if (merge)
      storage.merge();

    return storage;
  }

private:
  std::unique_ptr<kmers::KMerSplitter<Seq>> splitter_;
  fs::TmpDir work_dir_;

  size_t MergeKMers(const std::string &ifname, const std::string &ofname) {
    MMappedRecordArrayReader<typename Seq::DataType> ins(ifname, Seq::GetDataSize(this->k()), /* unlink */ true);

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
      adt::KMerVector<Seq> buf(this->k(), 1024*1024);
      size_t total = 0;
      while (!tree.empty()) {
          buf.clear();
          buf.push_back(tree.pop());
          size_t cnt = 1;

          while (cnt < buf.capacity()) {
            while (!tree.empty() &&
                   adt::array_equal_to<typename Seq::DataType>()(buf.back(), tree.top()))
              tree.replay();

            if (tree.empty())
              break;

            buf.push_back(tree.top());
            tree.replay();
            cnt += 1;
          }

          // Handle the last value
          while (!tree.empty() &&
                 adt::array_equal_to<typename Seq::DataType>()(buf.back(), tree.top()))
            tree.replay();

          total += buf.size();

          FILE *g = fopen(ofname.c_str(), "ab");
          if (!g)
            FATAL_ERROR("Cannot open temporary file " << ofname << " for writing");
          size_t res = fwrite(buf.data(), buf.el_data_size(), buf.size(), g);
          if (res != buf.size())
            FATAL_ERROR("I/O error! Incomplete write! Reason: " << strerror(errno) << ". Error code: " << errno);
          fclose(g);
      }

      return total;
    } else {
      // Sort the stuff
      libcxx::sort(ins.begin(), ins.end(), adt::array_less<typename Seq::DataType>());

      // FIXME: Use something like parallel version of unique_copy but with explicit
      // resizing.
      auto it = std::unique(ins.begin(), ins.end(), adt::array_equal_to<typename Seq::DataType>());

      MMappedRecordArrayWriter<typename Seq::DataType> os(ofname, Seq::GetDataSize(this->k()));
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
  KMerIndexBuilder(unsigned num_threads)
      : num_buckets_(0), num_threads_(num_threads) {}

  KMerIndexBuilder(unsigned num_buckets, unsigned num_threads)
      : num_buckets_(num_buckets), num_threads_(num_threads) {}

  template<class KMerStorage>
  void BuildIndex(Index &index, const KMerStorage &kmer_storage) {
    TIME_TRACE_SCOPE("KMerIndexBuilder::BuildIndex(storage)");

    index.clear();

    auto segment_policy = kmer_storage.segment_policy();
    size_t segments = segment_policy.num_segments();
    index.segment_starts_.resize(segments + 1, 0);
    index.segment_policy_ = segment_policy;
    index.num_segments_ = segments;

    INFO("Building perfect hash indices");

    if (segments == kmer_storage.num_buckets()) {
      // Build segmented index joining all buckets
      TIME_TRACE_SCOPE("KMerIndexBuilder::BuildIndex(storage, segmented)");
      for (size_t i = 0; i < kmer_storage.num_buckets(); ++i)
        // Use of gamma = 4 implies ~5.7 bits per k-mer, however, faster construction and lookup
        index.index_.emplace_back(kmer_storage.bucket_size(i),
                                  Index::KMerDataIndex::ConflictPolicy::Ignore,
                                  /* gamma */ 4.0);
         
#     pragma omp parallel for shared(index) num_threads(num_threads_)
      for (size_t i = 0; i < kmer_storage.num_buckets(); ++i) {
          index.segment_starts_[i + 1] = kmer_storage.bucket_size(i);
          index.index_[i].build(boomphf::range(kmer_storage.bucket_begin(i), kmer_storage.bucket_end(i)));
      }
    } else {
      // Build single index parallel over buckets
      TIME_TRACE_SCOPE("KMerIndexBuilder::BuildIndex(storage, parallel)");
      VERIFY(segments == 1);
      // Use of gamma = 4 implies ~5.7 bits per k-mer, however, faster construction and lookup
      index.index_.emplace_back(kmer_storage.total_kmers(),
                                Index::KMerDataIndex::ConflictPolicy::Ignore,
                                /* gamma */ 4.0);
      using Range = decltype(boomphf::range(kmer_storage.bucket_begin(0), kmer_storage.bucket_end(0)));
      std::vector<Range> ranges;
      for (size_t i = 0; i < kmer_storage.num_buckets(); ++i)
        ranges.emplace_back(boomphf::range(kmer_storage.bucket_begin(i), kmer_storage.bucket_end(i)));
      index.index_[0].build(ranges, num_threads_);
    }

    // Finally, record the sizes of buckets.
    for (unsigned i = 1; i < segments; ++i)
      index.segment_starts_[i] += index.segment_starts_[i - 1];

    double bits_per_kmer = 8.0 * (double)index.mem_size() / (double)kmer_storage.total_kmers();
    INFO("Index built. Total " << kmer_storage.total_kmers() << " kmers, " << index.mem_size() << " bytes occupied (" << bits_per_kmer << " bits per kmer).");
    index.count_size();
  }

  auto BuildIndex(Index &index, KMerCounter<Seq> &counter, bool save_final = false) {
    TIME_TRACE_SCOPE("KMerIndexBuilder::BuildIndex(counter)");

    INFO("Building kmer index");
    // First, count the unique k-mers
    auto kmer_storage = counter.Count(num_buckets_, num_threads_);

    // Now, build the index itself
    BuildIndex(index, kmer_storage);

    if (save_final)
      kmer_storage.merge();

    return kmer_storage;
  }

 private:
  DECL_LOGGER("K-mer Index Building");
};

}
