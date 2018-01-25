//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "kmer_data.hpp"
#include "io/reads/read_processor.hpp"
#include "valid_kmer_generator.hpp"

#include "io/reads/ireadstream.hpp"
#include "config_struct_hammer.hpp"

#include "utils/kmer_mph/kmer_index_builder.hpp"
#include "utils/logger/logger.hpp"

#include "io/kmers/kmer_iterator.hpp"
#include "adt/cqf.hpp"
#include "adt/hll.hpp"

using namespace hammer;

class BufferFiller;

struct KMerComparator {
    bool operator()(const KMer &l, const KMer &r) const {
      for (size_t i = 0; i < KMer::DataSize ; ++i) {
        if (l.data()[i] != r.data()[i]) {
          return (l.data()[i] < r.data()[i]);
        }
      }

      return false;
    }
};


class HammerFilteringKMerSplitter : public utils::KMerSortingSplitter<hammer::KMer> {
 public:
  using typename utils::KMerSortingSplitter<hammer::KMer>::RawKMers;
  typedef std::function<bool(const KMer&)> KMerFilter;

  HammerFilteringKMerSplitter(std::string &work_dir,
                              KMerFilter filter = [](const KMer&) { return true; })
      : KMerSortingSplitter<hammer::KMer>(work_dir, hammer::K),
      filter_(std::move(filter)) {}

  RawKMers Split(size_t num_files, unsigned nthreads) override;

 private:
  KMerFilter filter_;

  friend class BufferFiller;
};

class BufferFiller {
  HammerFilteringKMerSplitter &splitter_;

 public:
  BufferFiller(HammerFilteringKMerSplitter &splitter)
      : splitter_(splitter) {}

  bool operator()(std::unique_ptr<Read> r) {
    int trim_quality = cfg::get().input_trim_quality;

    Read cr = *r;
    size_t sz = cr.trimNsAndBadQuality(trim_quality);
  
    if (sz < hammer::K)
      return false;
    
    unsigned thread_id = omp_get_thread_num();
    ValidKMerGenerator<hammer::K> gen(cr);
    bool stop = false;
    for (; gen.HasMore(); gen.Next()) {
      KMer seq = gen.kmer();
      if (!splitter_.filter_(seq))
        continue;

      stop |= splitter_.push_back_internal( seq, thread_id);
      stop |= splitter_.push_back_internal(!seq, thread_id);
    }

    return stop;
  }
};

HammerFilteringKMerSplitter::RawKMers HammerFilteringKMerSplitter::Split(size_t num_files, unsigned nthreads) {
  size_t reads_buffer_size = cfg::get().count_split_buffer;

  auto out = PrepareBuffers(num_files, nthreads, reads_buffer_size);

  size_t n = 15, processed = 0;
  BufferFiller filler(*this);
  for (const auto &reads : cfg::get().dataset.reads()) {
    INFO("Processing " << reads);
    ireadstream irs(reads, cfg::get().input_qvoffset);
    while (!irs.eof()) {
      hammer::ReadProcessor rp(nthreads);
      rp.Run(irs, filler);
      DumpBuffers(out);
      VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
      processed += rp.processed();

      if (processed >> n) {
        INFO("Processed " << processed << " reads");
        n += 1;
      }
    }
  }
  INFO("Total " << processed << " reads processed");

  this->ClearBuffers();

  return out;
}

static inline void Merge(KMerStat &lhs, const KMerStat &rhs) {
  lhs.set_count(lhs.count() + rhs.count());
  lhs.total_qual *= rhs.total_qual;
  lhs.qual += rhs.qual;
}

static void PushKMer(KMerData &data,
                     KMer kmer, const unsigned char *q, double prob) {
  size_t idx = data.checking_seq_idx(kmer);
  if (idx == -1ULL)
      return;
  KMerStat &kmc = data[idx];
  kmc.lock();
  Merge(kmc,
        KMerStat(1, (float)prob, q));
  kmc.unlock();
}

static void PushKMerRC(KMerData &data,
                       KMer kmer, const unsigned char *q, double prob) {
  unsigned char rcq[K];

  // Prepare RC kmer with quality.
  kmer = !kmer;
  for (unsigned i = 0; i < K; ++i)
    rcq[K - i - 1] = q[i];

  size_t idx = data.checking_seq_idx(kmer);
  if (idx == -1ULL)
      return;
  KMerStat &kmc = data[idx];
  kmc.lock();
  Merge(kmc,
        KMerStat(1, (float)prob, rcq));
  kmc.unlock();
}

class KMerDataFiller {
  KMerData &data_;

 public:
  KMerDataFiller(KMerData &data)
      : data_(data) {}

  bool operator()(std::unique_ptr<Read> r) {
    uint8_t trim_quality = (uint8_t)cfg::get().input_trim_quality;

    // FIXME: Get rid of this
    Read cr = *r;
    size_t sz = cr.trimNsAndBadQuality(trim_quality);

    if (sz < hammer::K)
      return false;

    ValidKMerGenerator<hammer::K> gen(cr);
    const char *q = cr.getQualityString().data();
    while (gen.HasMore()) {
      KMer kmer = gen.kmer();
      const unsigned char *kq = (const unsigned char*)(q + gen.pos() - 1);

      PushKMer(data_, kmer, kq, 1 - gen.correct_probability());
      PushKMerRC(data_, kmer, kq, 1 - gen.correct_probability());

      gen.Next();
    }

    return false;
  }
};

class KMerMultiplicityCounter {
  qf::cqf_with_hasher<KMer> cqf_;

  public:
  KMerMultiplicityCounter(size_t size)
      : cqf_(size, [](const KMer &k) { return k.GetHash(); }) {}

  ~KMerMultiplicityCounter() {}

    bool operator()(std::unique_ptr<Read> r) {
      uint8_t trim_quality = (uint8_t)cfg::get().input_trim_quality;

      // FIXME: Get rid of this
      Read cr = *r;
      size_t sz = cr.trimNsAndBadQuality(trim_quality);

      if (sz < hammer::K)
        return false;

      ValidKMerGenerator<hammer::K> gen(cr);
      for (; gen.HasMore(); gen.Next()) {
          KMer kmer = gen.kmer();

          cqf_.add(kmer);
          cqf_.add(!kmer);
      }

      return false;
  }

  size_t count(const KMer &k) const {
      return cqf_.lookup(k);
  }
};

class KMerCountEstimator {
  std::vector<hll::hll_with_hasher<KMer>> hll_;

  public:
  KMerCountEstimator(unsigned thread_num) {
      hll_.reserve(thread_num);
      for (unsigned i = 0; i < thread_num; ++i)
          hll_.emplace_back([](const KMer &k) { return k.GetHash(); });
  }

  ~KMerCountEstimator() {}

    bool operator()(std::unique_ptr<Read> r) {
      uint8_t trim_quality = (uint8_t)cfg::get().input_trim_quality;

      // FIXME: Get rid of this
      Read cr = *r;
      size_t sz = cr.trimNsAndBadQuality(trim_quality);

      if (sz < hammer::K)
        return false;

      ValidKMerGenerator<hammer::K> gen(cr);
      for (; gen.HasMore(); gen.Next()) {
          KMer kmer = gen.kmer();
          auto &hll = hll_[omp_get_thread_num()];

          hll.add(kmer);
          hll.add(!kmer);
      }

      return false;
  }

  std::pair<double, bool> cardinality() const {
      return hll_[0].cardinality();
  }

  void merge() {
      for (size_t i = 1; i < hll_.size(); ++i) {
          hll_[0].merge(hll_[i]);
          hll_[i].clear();
      }
  }
};

void KMerDataCounter::BuildKMerIndex(KMerData &data) {
  // Build the index
  std::string workdir = cfg::get().input_working_dir;

  // Optionally perform a filtering step
  size_t kmers = 0;
  typename KMerData::traits::ResultFile final_kmers;
  if (cfg::get().count_filter_singletons) {
      size_t buffer_size;
      {
          INFO("Estimating k-mer count");

          size_t n = 15, processed = 0;
          KMerCountEstimator mcounter(omp_get_max_threads());
          for (const auto &reads : cfg::get().dataset.reads()) {
              INFO("Processing " << reads);
              ireadstream irs(reads, cfg::get().input_qvoffset);
              while (!irs.eof()) {
                  hammer::ReadProcessor rp(omp_get_max_threads());
                  rp.Run(irs, mcounter);
                  VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
                  processed += rp.processed();

                  if (processed >> n) {
                      INFO("Processed " << processed << " reads");
                      n += 1;
                  }
              }
          }
          INFO("Total " << processed << " reads processed");
          mcounter.merge();
          std::pair<double, bool> res = mcounter.cardinality();
          if (res.second == false) {
              buffer_size = cfg::get().count_split_buffer;
              if (buffer_size == 0) buffer_size = 512ull * 1024 * 1024;
          } else {
              INFO("Estimated " << size_t(res.first) << " distinct kmers");
              buffer_size = size_t(res.first);
          }
      }

      INFO("Filtering singleton k-mers");

      KMerMultiplicityCounter mcounter(buffer_size);

      size_t n = 15, processed = 0;
      for (const auto &reads : cfg::get().dataset.reads()) {
          INFO("Processing " << reads);
          ireadstream irs(reads, cfg::get().input_qvoffset);
          while (!irs.eof()) {
              hammer::ReadProcessor rp(omp_get_max_threads());
              rp.Run(irs, mcounter);
              VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
              processed += rp.processed();

              if (processed >> n) {
                  INFO("Processed " << processed << " reads");
                  n += 1;
              }
          }
      }
      INFO("Total " << processed << " reads processed");

      // FIXME: Reduce code duplication
      HammerFilteringKMerSplitter splitter(workdir,
                                           [&] (const KMer &k) { return mcounter.count(k) > 1; });
      utils::KMerDiskCounter<hammer::KMer> counter(workdir, splitter);

      kmers = utils::KMerIndexBuilder<HammerKMerIndex>(num_files_, omp_get_max_threads()).BuildIndex(data.index_, counter, /* save final */ true);
      final_kmers = counter.final_kmers_file();
  } else {
      HammerFilteringKMerSplitter splitter(workdir);
      utils::KMerDiskCounter<hammer::KMer> counter(workdir, splitter);

      kmers = utils::KMerIndexBuilder<HammerKMerIndex>(num_files_, omp_get_max_threads()).BuildIndex(data.index_, counter, /* save final */ true);
      final_kmers = counter.final_kmers_file();
  }

  // Check, whether we'll ever have enough memory for running BH and bail out earlier
  double needed = 1.25 * (double)kmers * (sizeof(KMerStat) + sizeof(hammer::KMer));
  if (needed > (double) utils::get_memory_limit())
      FATAL_ERROR("The reads contain too many k-mers to fit into available memory. You need approx. "
                  << needed / 1024.0 / 1024.0 / 1024.0
                  << "GB of free RAM to assemble your dataset");

  {
    INFO("Arranging kmers in hash map order");
    data.kmers_.set_size(kmers);
    data.kmers_.set_data(new hammer::KMer::DataType[kmers * hammer::KMer::GetDataSize(hammer::K)]);

    unsigned nthreads = std::min(cfg::get().count_merge_nthreads, cfg::get().general_max_nthreads);
    auto kmers_its = io::make_kmer_iterator<hammer::KMer>(*final_kmers, hammer::K, 16*nthreads);

#   pragma omp parallel for num_threads(nthreads) schedule(guided)
    for (size_t i = 0; i < kmers_its.size(); ++i) {
        auto &kmer_it = kmers_its[i];
        for (; kmer_it.good(); ++kmer_it) {
            size_t kidx = data.index_.seq_idx(hammer::KMer(hammer::K, *kmer_it));
            memcpy(data.kmers_[kidx].data(), *kmer_it, hammer::KMer::TotalBytes);
        }
    }
  }
}

void KMerDataCounter::FillKMerData(KMerData &data) {
  // Now use the index to fill the kmer quality information.
  INFO("Collecting K-mer information, this takes a while.");
  data.data_.resize(data.kmers_.size());

  KMerDataFiller filler(data);
  const auto& dataset = cfg::get().dataset;
  for (auto I = dataset.reads_begin(), E = dataset.reads_end(); I != E; ++I) {
    INFO("Processing " << *I);
    ireadstream irs(*I, cfg::get().input_qvoffset);
    hammer::ReadProcessor rp(omp_get_max_threads());
    rp.Run(irs, filler);
    VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
  }

  INFO("Collection done, postprocessing.");

  size_t singletons = 0;
  for (size_t i = 0; i < data.size(); ++i) {
    VERIFY(data[i].count());

    // Make sure all the kmers are marked as 'Bad' in the beginning
    data[i].mark_bad();

    if (data[i].count() == 1)
      singletons += 1;
  }

  INFO("There are " << data.size() << " kmers in total. "
       "Among them " << singletons << " (" <<  100.0 * (double)singletons / (double)data.size() << "%) are singletons.");
}
