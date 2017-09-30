//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "kmer_data.hpp"
#include "config_struct.hpp"
#include "valid_hkmer_generator.hpp"

#include "utils/kmer_mph/kmer_index_builder.hpp"

#include <mutex>
#include <random>
#include "io/kmers/mmapped_writer.hpp"
#include "io/reads/file_reader.hpp"
#include "io/reads/read_processor.hpp"

using namespace hammer;

class BufferFiller;

class HammerKMerSplitter : public utils::KMerSortingSplitter<HKMer> {
 public:
  using typename utils::KMerSortingSplitter<HKMer>::RawKMers;

  HammerKMerSplitter(const std::string &work_dir)
      : KMerSortingSplitter<HKMer>(work_dir, hammer::K) {}

  RawKMers Split(size_t num_files, unsigned nthreads) override;

  friend class BufferFiller;
};

class BufferFiller {
  size_t processed_;
  HammerKMerSplitter &splitter_;

 public:
  BufferFiller(HammerKMerSplitter &splitter)
      : processed_(0), splitter_(splitter) {}

  size_t processed() const { return processed_; }

  bool operator()(std::unique_ptr<io::SingleRead> r) {
    ValidHKMerGenerator<hammer::K> gen(*r);
    unsigned thread_id = omp_get_thread_num();

#pragma omp atomic
    processed_ += 1;

    bool stop = false;
    while (gen.HasMore()) {
      HKMer seq = gen.kmer();

      stop |= splitter_.push_back_internal(seq, thread_id);
      stop |= splitter_.push_back_internal(!seq, thread_id);

      gen.Next();
    }

    return stop;
  }
};

HammerKMerSplitter::RawKMers HammerKMerSplitter::Split(size_t num_files, unsigned nthreads) {
  size_t reads_buffer_size = cfg::get().count_split_buffer;

  auto out = PrepareBuffers(num_files, nthreads, reads_buffer_size);

  size_t n = 15;
  BufferFiller filler(*this);
  for (const auto &reads : cfg::get().dataset.reads()) {
    INFO("Processing " << reads);
    io::FileReadStream irs(reads, io::PhredOffset);
    hammer::ReadProcessor rp(nthreads);
    while (!irs.eof()) {
      rp.Run(irs, filler);
      DumpBuffers(out);
      VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");

      if (filler.processed() >> n) {
        INFO("Processed " << filler.processed() << " reads");
        n += 1;
      }
    }
  }
  INFO("Processed " << filler.processed() << " reads");

  this->ClearBuffers();

  return out;
}

static inline void Merge(KMerStat &lhs, const KMerStat &rhs) {
  if (lhs.count == 0) lhs.kmer = rhs.kmer;

  lhs.count += rhs.count;
  lhs.qual += rhs.qual;
}

static void PushKMer(KMerData &data, HKMer kmer, double qual) {
  KMerStat &kmc = data[kmer];
  kmc.lock();
  Merge(kmc, KMerStat(1, kmer, (float)qual));
  kmc.unlock();
}

static void PushKMerRC(KMerData &data, HKMer kmer, double qual) {
  PushKMer(data, !kmer, qual);
}

class KMerDataFiller {
  KMerData &Data;
  mutable std::default_random_engine RandomEngine;
  mutable std::uniform_real_distribution<double> UniformRandGenerator;
  mutable std::mutex Lock;
  double SampleRate;

 public:
  KMerDataFiller(KMerData &data, double sampleRate = 1.0)
      : Data(data),
        RandomEngine(42),
        UniformRandGenerator(0, 1),
        SampleRate(sampleRate) {}

  double NextUniform() const {
    std::lock_guard<std::mutex> guard(Lock);
    return UniformRandGenerator(RandomEngine);
  }

  bool operator()(std::unique_ptr<io::SingleRead> &&r) const {
    ValidHKMerGenerator<hammer::K> gen(*r);

    // tiny quality regularization
    const double decay = 0.9999;
    double prior = 1.0;

    bool skipRead = SampleRate < 1.0 && (NextUniform() > SampleRate);

    if (skipRead) {
      return false;
    }

    while (gen.HasMore()) {
      const HKMer kmer = gen.kmer();
      const double p = gen.correct_probability();
      gen.Next();

      assert(p < 1.0);
      assert(p >= 0);
      const double correct = p * prior;

      prior *= decay;
      {
        PushKMer(Data, kmer, log(1 - correct));

        PushKMerRC(Data, kmer, log(1 - correct));
      }
    }
    // Do not stop
    return false;
  }
};

void KMerDataCounter::FillKMerData(KMerData &data) {
  HammerKMerSplitter splitter(cfg::get().working_dir);
  utils::KMerDiskCounter<hammer::HKMer> counter(cfg::get().working_dir, splitter);

  size_t sz = utils::KMerIndexBuilder<HammerKMerIndex>(num_files_, cfg::get().max_nthreads).BuildIndex(data.index_, counter);

  // Now use the index to fill the kmer quality information.
  INFO("Collecting K-mer information, this takes a while.");
  data.data_.resize(sz);

  const auto &dataset = cfg::get().dataset;
  for (auto it = dataset.reads_begin(), et = dataset.reads_end(); it != et;
       ++it) {
    INFO("Processing " << *it);
    io::FileReadStream irs(*it, io::PhredOffset);
    KMerDataFiller filler(data, cfg::get().sample_rate);
    hammer::ReadProcessor(cfg::get().max_nthreads).Run(irs, filler);
  }

  INFO("Collection done, postprocessing.");

  size_t singletons = 0;
  size_t skipped = 0;
  for (size_t i = 0; i < data.size(); ++i) {
    if (data[i].count == 1) {
      singletons += 1;
    }
    if (data[i].count == 0) {
      skipped += 1;
    }
  }

  INFO("Merge done. There are "
       << data.size()
       << " kmers in total. "
          "Among them "
       << singletons << " (" << 100.0 * double(singletons) / double(data.size())
       << "%) are singletons."
       << "Among them " << skipped << " ("
       << 100.0 * double(skipped) / double(data.size())
       << "%) are skipped during sampling.");
}
