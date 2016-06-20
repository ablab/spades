#pragma once
//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io/kmers_io/mmapped_reader.hpp"
#include "io/kmers_io/mmapped_writer.hpp"
#include "utils/adt/pointer_iterator.hpp"
#include "utils/adt/kmer_vector.hpp"

#include "mphf.hpp"
#include "base_hash.hpp"
#include "hypergraph.hpp"
#include "hypergraph_sorter_seq.hpp"

#include "dev_support/openmp_wrapper.h"

#include "dev_support/logger/logger.hpp"
#include "dev_support/path_helper.hpp"

#include "dev_support/memory_limit.hpp"

#include <libcxx/sort.hpp>

#include <algorithm>
#ifdef USE_GLIBCXX_PARALLEL
#include <parallel/algorithm>
#endif
#include <fstream>
#include <vector>
#include <cmath>

#include "config.hpp"

#ifdef SPADES_USE_JEMALLOC
# include <jemalloc/jemalloc.h>
#endif

template<class Index>
class KMerIndexBuilder;

template<class Seq>
struct kmer_index_traits {
  typedef Seq SeqType;
  typedef MMappedRecordArrayReader<typename Seq::DataType> RawKMerStorage;
  typedef MMappedRecordArrayReader<typename Seq::DataType> FinalKMerStorage;
  typedef typename RawKMerStorage::iterator             raw_data_iterator;
  typedef typename RawKMerStorage::const_iterator       raw_data_const_iterator;
  typedef typename RawKMerStorage::iterator::value_type KMerRawData;
  typedef typename RawKMerStorage::iterator::reference  KMerRawReference;
  typedef typename RawKMerStorage::const_iterator::reference  KMerRawConstReference;

  struct raw_equal_to {
    bool operator()(const Seq &lhs, const KMerRawReference rhs) {
      return (array_equal_to<typename Seq::DataType>()(lhs.data(), lhs.data_size(), rhs));
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
    uint64_t operator()(const Seq &k) const{
      return typename Seq::hash()(k);
    }
    uint64_t operator()(const KMerRawReference k) const {
      return typename Seq::hash()(k.data(), k.size());
    }
  };

  struct KMerRawReferenceAdaptor {
      emphf::byte_range_t operator()(const KMerRawReference k) const {
          const uint8_t * data = (const uint8_t*)k.data();
          return std::make_pair(data, data + k.data_size());
      }
  };

  struct KMerSeqAdaptor {
      emphf::byte_range_t operator()(const Seq &k) const {
          const uint8_t * data = (const uint8_t*)k.data();
          return std::make_pair(data, data + k.data_size() * sizeof(typename Seq::DataType));
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

template<class traits>
class KMerIndex {
 public:
  typedef traits kmer_index_traits;
  typedef typename traits::SeqType          KMerSeq;
  typedef typename traits::hash_function    hash_function;
  typedef typename traits::KMerRawData      KMerRawData;
  typedef typename traits::KMerRawReference KMerRawReference;
  typedef size_t IdxType;

 private:
  using KMerDataIndex = emphf::mphf<emphf::city_hasher>;
  typedef KMerIndex __self;

 public:
  KMerIndex(): index_(NULL), num_buckets_(0), size_(0) {}

  KMerIndex(const KMerIndex&) = delete;
  KMerIndex& operator=(const KMerIndex&) = delete;

  ~KMerIndex() { clear(); }

  void clear() {
    num_buckets_ = 0;
    bucket_starts_.clear();

    delete[] index_;
    index_ = NULL;
  }

  size_t mem_size() {
    size_t sz = 0;
    for (size_t i = 0; i < num_buckets_; ++i)
      sz += index_[i].mem_size();

    return sz;
  }

  void count_size() {
      if (index_ == NULL)
          return;
      size_ = 0;
      for (size_t i = 0; i < num_buckets_; i++)
          size_ += index_[i].size();
  }

  size_t size() const {
      return size_;
  }

  size_t seq_idx(const KMerSeq &s) const {
    size_t bucket = seq_bucket(s);

    return bucket_starts_[bucket] +
            index_[bucket].lookup(s, typename traits::KMerSeqAdaptor());
  }

  size_t raw_seq_idx(const KMerRawReference data) const {
    size_t bucket = raw_seq_bucket(data);

    return bucket_starts_[bucket] +
            index_[bucket].lookup(data, typename traits::KMerRawReferenceAdaptor());
  }

  template<class Writer>
  void serialize(Writer &os) const {
    os.write((char*)&num_buckets_, sizeof(num_buckets_));
    for (size_t i = 0; i < num_buckets_; ++i)
      index_[i].save(os);
    os.write((char*)&bucket_starts_[0], (num_buckets_ + 1) * sizeof(bucket_starts_[0]));
  }

  template<class Reader>
  void deserialize(Reader &is) {
    clear();

    is.read((char*)&num_buckets_, sizeof(num_buckets_));

    index_ = new KMerDataIndex[num_buckets_];
    for (size_t i = 0; i < num_buckets_; ++i)
      index_[i].load(is);

    bucket_starts_.resize(num_buckets_ + 1);
    is.read((char*)&bucket_starts_[0], (num_buckets_ + 1) * sizeof(bucket_starts_[0]));
    count_size();
  }

  void swap(KMerIndex<traits> &other) {
    std::swap(index_, other.index_);
    std::swap(num_buckets_, other.num_buckets_);
    std::swap(size_, other.size_);
    std::swap(bucket_starts_, other.bucket_starts_);
  }

 private:
  KMerDataIndex *index_;

  size_t num_buckets_;
  std::vector<size_t> bucket_starts_;
  size_t size_;

  size_t seq_bucket(const KMerSeq &s) const {
    return hash_function()(s) % num_buckets_;
  }
  size_t raw_seq_bucket(const KMerRawReference data) const {
    return hash_function()(data) % num_buckets_;
  }

  friend class KMerIndexBuilder<__self>;
};

template<class Seq>
class KMerSplitter {
 public:
  typedef typename Seq::hash hash_function;

  KMerSplitter(const std::string &work_dir, unsigned K, uint32_t seed = 0)
      : work_dir_(work_dir), K_(K), seed_(seed) {}

  virtual ~KMerSplitter() {}

  virtual path::files_t Split(size_t num_files) = 0;

  unsigned K() const { return K_; }

 protected:
  const std::string &work_dir_;
  hash_function hash_;
  unsigned K_;
  uint32_t seed_;

  std::string GetRawKMersFname(unsigned suffix) const {
    return path::append_path(work_dir_, "kmers.raw." + std::to_string(suffix));
  }

  unsigned GetFileNumForSeq(const Seq &s, unsigned total) const {
    return (unsigned)(hash_(s, seed_) % total);
  }

  DECL_LOGGER("K-mer Splitting");
};

template<class Seq>
class KMerSortingSplitter : public KMerSplitter<Seq> {
 public:
  KMerSortingSplitter(const std::string &work_dir, unsigned K, uint32_t seed = 0)
      : KMerSplitter<Seq>(work_dir, K, seed) {}

 protected:
  using SeqKMerVector = KMerVector<Seq>;
  using KMerBuffer = std::vector<SeqKMerVector>;

  void DumpBuffers(size_t num_files, size_t nthreads,
                   std::vector<KMerBuffer> &buffers,
                   const path::files_t &ostreams) const {
#   pragma omp parallel for
    for (unsigned k = 0; k < num_files; ++k) {
      size_t sz = 0;
      for (size_t i = 0; i < nthreads; ++i)
        sz += buffers[i][k].size();

      KMerVector<Seq> SortBuffer(this->K_, sz);
      for (size_t i = 0; i < nthreads; ++i) {
        KMerBuffer &entry = buffers[i];
        for (size_t j = 0; j < entry[k].size(); ++j)
          SortBuffer.push_back(entry[k][j]);
      }
      libcxx::sort(SortBuffer.begin(), SortBuffer.end(), typename KMerVector<Seq>::less2_fast());
      auto it = std::unique(SortBuffer.begin(), SortBuffer.end(), typename KMerVector<Seq>::equal_to());

#     pragma omp critical
      {
        FILE *f = fopen(ostreams[k].c_str(), "ab");
        VERIFY_MSG(f, "Cannot open temporary file to write");
        fwrite(SortBuffer.data(), SortBuffer.el_data_size(), it - SortBuffer.begin(), f);
        fclose(f);
      }
    }

    for (unsigned i = 0; i < nthreads; ++i) {
      for (unsigned j = 0; j < num_files; ++j) {
        buffers[i][j].clear();
      }
    }
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
