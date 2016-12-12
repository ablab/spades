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
#include "common/adt/pointer_iterator.hpp"
#include "common/adt/kmer_vector.hpp"

#include "utils/openmp_wrapper.h"

#include "utils/logger/logger.hpp"
#include "utils/path_helper.hpp"

#include "utils/memory_limit.hpp"
#include "utils/file_limit.hpp"

#include "adt/iterator_range.hpp"
#include "adt/loser_tree.hpp"

#include "mphf.hpp"
#include "base_hash.hpp"
#include "hypergraph.hpp"
#include "hypergraph_sorter_seq.hpp"

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

template<class Seq>
class KMerSplitter {
 public:
  typedef typename Seq::hash hash_function;

  KMerSplitter(const std::string &work_dir, unsigned K, uint32_t seed = 0)
      : work_dir_(work_dir), K_(K), seed_(seed) {}

  virtual ~KMerSplitter() {}

  virtual path::files_t Split(size_t num_files) = 0;

  size_t kmer_size() const {
    return Seq::GetDataSize(K_) * sizeof(typename Seq::DataType);
  }

  unsigned K() const { return K_; }

 protected:
  const std::string &work_dir_;
  hash_function hash_;
  unsigned K_;
  uint32_t seed_;

  DECL_LOGGER("K-mer Splitting");
};

template<class Seq>
class KMerSortingSplitter : public KMerSplitter<Seq> {
 public:
  KMerSortingSplitter(const std::string &work_dir, unsigned K, uint32_t seed = 0)
      : KMerSplitter<Seq>(work_dir, K, seed), cell_size_(0), num_files_(0) {}

 protected:
  using SeqKMerVector = KMerVector<Seq>;
  using KMerBuffer = std::vector<SeqKMerVector>;

  std::vector<KMerBuffer> kmer_buffers_;
  size_t cell_size_;
  size_t num_files_;

  path::files_t PrepareBuffers(size_t num_files, unsigned nthreads, size_t reads_buffer_size) {
    num_files_ = num_files;
    
    // Determine the set of output files
    path::files_t out;
    for (unsigned i = 0; i < num_files_; ++i)
      out.push_back(this->GetRawKMersFname(i));

    size_t file_limit = num_files_ + 2*nthreads;
    size_t res = limit_file(file_limit);
    if (res < file_limit) {
      WARN("Failed to setup necessary limit for number of open files. The process might crash later on.");
      WARN("Do 'ulimit -n " << file_limit << "' in the console to overcome the limit");
    }

    if (reads_buffer_size == 0) {
      reads_buffer_size = 536870912ull;
      size_t mem_limit =  (size_t)((double)(get_free_memory()) / (nthreads * 3));
      INFO("Memory available for splitting buffers: " << (double)mem_limit / 1024.0 / 1024.0 / 1024.0 << " Gb");
      reads_buffer_size = std::min(reads_buffer_size, mem_limit);
    }
    cell_size_ = reads_buffer_size / (num_files_ * this->kmer_size());
    // Set sane minimum cell size
    if (cell_size_ < 16384)
      cell_size_ = 16384;

    INFO("Using cell size of " << cell_size_);
    kmer_buffers_.resize(nthreads);
    for (unsigned i = 0; i < nthreads; ++i) {
      KMerBuffer &entry = kmer_buffers_[i];
      entry.resize(num_files_, KMerVector<Seq>(this->K_, (size_t) (1.1 * (double) cell_size_)));
    }

    return out;
  }
  
  bool push_back_internal(const Seq &seq, unsigned thread_id) {
    KMerBuffer &entry = kmer_buffers_[thread_id];

    size_t idx = this->GetFileNumForSeq(seq, (unsigned)num_files_);
    entry[idx].push_back(seq);
    return entry[idx].size() > cell_size_;
  }
  
  void DumpBuffers(const path::files_t &ostreams) {
    VERIFY(ostreams.size() == num_files_ && kmer_buffers_[0].size() == num_files_);

#   pragma omp parallel for
    for (unsigned k = 0; k < num_files_; ++k) {
      // Below k is thread id!
      
      size_t sz = 0;
      for (size_t i = 0; i < kmer_buffers_.size(); ++i)
        sz += kmer_buffers_[i][k].size();

      KMerVector<Seq> SortBuffer(this->K_, sz);
      for (auto & entry : kmer_buffers_) {
        const auto &buffer = entry[k];
        for (size_t j = 0; j < buffer.size(); ++j)
          SortBuffer.push_back(buffer[j]);
      }
      libcxx::sort(SortBuffer.begin(), SortBuffer.end(), typename KMerVector<Seq>::less2_fast());
      auto it = std::unique(SortBuffer.begin(), SortBuffer.end(), typename KMerVector<Seq>::equal_to());

#     pragma omp critical
      {
        size_t cnt =  it - SortBuffer.begin();

        // Write k-mers
        FILE *f = fopen(ostreams[k].c_str(), "ab");
        VERIFY_MSG(f, "Cannot open temporary file to write");
        fwrite(SortBuffer.data(), SortBuffer.el_data_size(), cnt, f);
        fclose(f);

        // Write index
        f = fopen((ostreams[k] + ".idx").c_str(), "ab");
        VERIFY_MSG(f, "Cannot open temporary file to write");
        fwrite(&cnt, sizeof(cnt), 1, f);
        fclose(f);
      }
    }

    for (auto & entry : kmer_buffers_)
      for (auto & eentry : entry)
        eentry.clear();
  }

  void ClearBuffers() {
    for (auto & entry : kmer_buffers_)
      for (auto & eentry : entry) {
        eentry.clear();
        eentry.shrink_to_fit();
      }
  }
  
  std::string GetRawKMersFname(unsigned suffix) const {
    return path::append_path(this->work_dir_, "kmers.raw." + std::to_string(suffix));
  }

  unsigned GetFileNumForSeq(const Seq &s, unsigned total) const {
    return (unsigned)(this->hash_(s, this->seed_) % total);
  }

};

template<class Seq, class traits = kmer_index_traits<Seq> >
class KMerCounter {
 public:
  typedef typename traits::raw_data_iterator       iterator;
  typedef typename traits::raw_data_const_iterator const_iterator;
  typedef typename traits::RawKMerStorage          RawKMerStorage;
  typedef typename traits::FinalKMerStorage        FinalKMerStorage;

  virtual size_t kmer_size() const = 0;

  virtual size_t Count(unsigned num_buckets, unsigned num_threads) = 0;
  virtual size_t CountAll(unsigned num_buckets, unsigned num_threads, bool merge = true) = 0;
  virtual void MergeBuckets(unsigned num_buckets) = 0;

  virtual std::unique_ptr<RawKMerStorage> GetBucket(size_t idx, bool unlink = true) = 0;
  virtual std::unique_ptr<FinalKMerStorage> GetFinalKMers() = 0;

  virtual ~KMerCounter() {}

protected:
  DECL_LOGGER("K-mer Counting");
};

template<class Seq, class traits = kmer_index_traits<Seq> >
class KMerDiskCounter : public KMerCounter<Seq> {
  typedef KMerCounter<Seq, traits> __super;
  typedef typename traits::RawKMerStorage BucketStorage;
public:
  KMerDiskCounter(const std::string &work_dir, KMerSplitter<Seq> &splitter)
      : work_dir_(work_dir), splitter_(splitter) {
    std::string prefix = path::append_path(work_dir, "kmers_XXXXXX");
    char *tempprefix = strcpy(new char[prefix.length() + 1], prefix.c_str());
    VERIFY_MSG(-1 != (fd_ = ::mkstemp(tempprefix)), "Cannot create temporary file");
    kmer_prefix_ = tempprefix;
    delete[] tempprefix;
  }

  ~KMerDiskCounter() {
    ::close(fd_);
    ::unlink(kmer_prefix_.c_str());
  }

  size_t kmer_size() const override {
    return Seq::GetDataSize(splitter_.K()) * sizeof(typename Seq::DataType);
  }

  std::unique_ptr<BucketStorage> GetBucket(size_t idx, bool unlink = true) override {
    unsigned K = splitter_.K();
    return std::unique_ptr<BucketStorage>(new BucketStorage(GetMergedKMersFname((unsigned)idx), Seq::GetDataSize(K), unlink));
  }

  size_t Count(unsigned num_buckets, unsigned num_threads) override {
    unsigned K = splitter_.K();

    // Split k-mers into buckets.
    path::files_t raw_kmers = splitter_.Split(num_buckets * num_threads);

    INFO("Starting k-mer counting.");
    size_t kmers = 0;
#   pragma omp parallel for shared(raw_kmers) num_threads(num_threads) schedule(dynamic) reduction(+:kmers)
    for (unsigned iFile = 0; iFile < raw_kmers.size(); ++iFile) {
      kmers += MergeKMers(raw_kmers[iFile], GetUniqueKMersFname(iFile), K);
    }
    INFO("K-mer counting done. There are " << kmers << " kmers in total. ");

    INFO("Merging temporary buckets.");
    for (unsigned i = 0; i < num_buckets; ++i) {
      std::string ofname = GetMergedKMersFname(i);
      std::ofstream ofs(ofname.c_str(), std::ios::out | std::ios::binary);
      for (unsigned j = 0; j < num_threads; ++j) {
        BucketStorage ins(GetUniqueKMersFname(i + j * num_buckets), Seq::GetDataSize(K), /* unlink */ true);
        ofs.write((const char*)ins.data(), ins.data_size());
      }
    }

    return kmers;
  }

  void MergeBuckets(unsigned num_buckets) override {
    unsigned K = splitter_.K();

    INFO("Merging final buckets.");

    MMappedRecordArrayWriter<typename Seq::DataType> os(GetFinalKMersFname(), Seq::GetDataSize(K));
    std::string ofname = GetFinalKMersFname();
    std::ofstream ofs(ofname.c_str(), std::ios::out | std::ios::binary);
    for (unsigned j = 0; j < num_buckets; ++j) {
      auto bucket = GetBucket(j, /* unlink */ true);
      ofs.write((const char*)bucket->data(), bucket->data_size());
    }
    ofs.close();
  }

  size_t CountAll(unsigned num_buckets, unsigned num_threads, bool merge = true) override {
    size_t kmers = Count(num_buckets, num_threads);
    if (merge)
      MergeBuckets(num_buckets);

    return kmers;
  }

  std::unique_ptr<typename __super::FinalKMerStorage> GetFinalKMers() override {
    unsigned K = splitter_.K();
    return std::unique_ptr<typename __super::FinalKMerStorage>(new typename __super::FinalKMerStorage(GetFinalKMersFname(), Seq::GetDataSize(K), /* unlink */ true));
  }

  std::string GetMergedKMersFname(unsigned suffix) const {
    return kmer_prefix_ + ".merged." + std::to_string(suffix);
  }

  std::string GetFinalKMersFname() const {
    return kmer_prefix_ + ".final";
  }

private:
  std::string work_dir_;
  KMerSplitter<Seq> &splitter_;
  int fd_;
  std::string kmer_prefix_;

  std::string GetUniqueKMersFname(unsigned suffix) const {
    return kmer_prefix_ + ".unique." + std::to_string(suffix);
  }

  size_t MergeKMers(const std::string &ifname, const std::string &ofname,
                    unsigned K) {
    MMappedRecordArrayReader<typename Seq::DataType> ins(ifname, Seq::GetDataSize(K), /* unlink */ true);

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
        VERIFY(std::is_sorted(beg, end, array_less<typename Seq::DataType>()));
        beg = end;
      }

      // Construct tree on top entries of runs
      adt::loser_tree<decltype(beg),
                      array_less<typename Seq::DataType>> tree(ranges);

      if (tree.empty()) {
        FILE *g = fopen(ofname.c_str(), "ab");
        VERIFY_MSG(g, "Cannot open temporary file to write");
        fclose(g);
        return 0;
      }

      // Write it down!
      KMerVector<Seq> buf(K, 1024*1024);
      auto pval = tree.pop();
      size_t total = 0;
      while (!tree.empty()) {
          buf.clear();
          for (size_t cnt = 0; cnt < buf.capacity() && !tree.empty(); ) {
              auto cval = tree.pop();
              if (!array_equal_to<typename Seq::DataType>()(pval, cval)) {
                  buf.push_back(pval);
                  pval = cval;
                  cnt += 1;
              }
          }
          total += buf.size();

          FILE *g = fopen(ofname.c_str(), "ab");
          VERIFY_MSG(g, "Cannot open temporary file to write");
          fwrite(buf.data(), buf.el_data_size(), buf.size(), g);
          fclose(g);
      }

      // Handle very last value
      {
        FILE *g = fopen(ofname.c_str(), "ab");
        VERIFY_MSG(g, "Cannot open temporary file to write");
        fwrite(pval.data(), pval.data_size(), 1, g);
        fclose(g);
        total += 1;
      }
      
      return total;
    } else {
      // Sort the stuff
      libcxx::sort(ins.begin(), ins.end(), array_less<typename Seq::DataType>());

      // FIXME: Use something like parallel version of unique_copy but with explicit
      // resizing.
      auto it = std::unique(ins.begin(), ins.end(), array_equal_to<typename Seq::DataType>());

      MMappedRecordArrayWriter<typename Seq::DataType> os(ofname, Seq::GetDataSize(K));
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

  std::string work_dir_;
  unsigned num_buckets_;
  unsigned num_threads_;

 public:
  KMerIndexBuilder(const std::string &workdir,
                   unsigned num_buckets, unsigned num_threads)
      : work_dir_(workdir), num_buckets_(num_buckets), num_threads_(num_threads) {}
  size_t BuildIndex(Index &out, KMerCounter<Seq> &counter,
                    bool save_final = false);

  unsigned num_buckets() const { return num_buckets_; }

 private:

  DECL_LOGGER("K-mer Index Building");
};

template<class Index>
size_t KMerIndexBuilder<Index>::BuildIndex(Index &index, KMerCounter<Seq> &counter,
                                           bool save_final) {
  index.clear();

  INFO("Building kmer index ");

  // First, count the unique k-mers
  size_t kmers = counter.Count(num_buckets_, num_threads_);

  index.num_buckets_ = num_buckets_;
  index.bucket_starts_.resize(num_buckets_ + 1);
  index.index_ = new typename KMerIndex<kmer_index_traits>::KMerDataIndex[num_buckets_];

  INFO("Building perfect hash indices");

  // Index building requires up to 40 bytes per k-mer. Limit number of threads depending on the memory limit.
  unsigned num_threads = num_threads_;
# ifdef SPADES_USE_JEMALLOC
  const size_t *cmem = 0;
  size_t clen = sizeof(cmem);

  je_mallctl("stats.cactive", &cmem, &clen, NULL, 0);
  size_t bucket_size = (36 * kmers + kmers * counter.kmer_size()) / num_buckets_;
  num_threads = std::min<unsigned>((unsigned) ((get_memory_limit() - *cmem) / bucket_size), num_threads);
  if (num_threads < 1)
    num_threads = 1;
  if (num_threads < num_threads_)
    WARN("Number of threads was limited down to " << num_threads << " in order to fit the memory limits during the index construction");
# endif

# pragma omp parallel for shared(index) num_threads(num_threads)
  for (unsigned iFile = 0; iFile < num_buckets_; ++iFile) {
    typename KMerIndex<kmer_index_traits>::KMerDataIndex &data_index = index.index_[iFile];
    auto bucket = counter.GetBucket(iFile, !save_final);
    size_t sz = bucket->end() - bucket->begin();
    index.bucket_starts_[iFile + 1] = sz;
    typename kmer_index_traits::KMerRawReferenceAdaptor adaptor;
    size_t max_nodes = (size_t(std::ceil(double(sz) * 1.23)) + 2) / 3 * 3;
    if (max_nodes >= uint64_t(1) << 32) {
        emphf::hypergraph_sorter_seq<emphf::hypergraph<uint64_t> > sorter;
        typename KMerIndex<kmer_index_traits>::KMerDataIndex(sorter,
                                                             sz, emphf::range(bucket->begin(), bucket->end()),
                                                             adaptor).swap(data_index);
    } else {
        emphf::hypergraph_sorter_seq<emphf::hypergraph<uint32_t> > sorter;
        typename KMerIndex<kmer_index_traits>::KMerDataIndex(sorter,
                                                             sz, emphf::range(bucket->begin(), bucket->end()),
                                                             adaptor).swap(data_index);
    }
  }

  // Finally, record the sizes of buckets.
  for (unsigned iFile = 1; iFile < num_buckets_; ++iFile)
    index.bucket_starts_[iFile] += index.bucket_starts_[iFile - 1];

  if (save_final)
    counter.MergeBuckets(num_buckets_);

  double bits_per_kmer = 8.0 * (double)index.mem_size() / (double)kmers;
  INFO("Index built. Total " << index.mem_size() << " bytes occupied (" << bits_per_kmer << " bits per kmer).");
  index.count_size();
  return kmers;
}
