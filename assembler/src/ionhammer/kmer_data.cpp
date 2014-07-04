//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "kmer_data.hpp"
#include "config_struct.hpp"
#include "valid_hkmer_generator.hpp"

#include "io/mmapped_writer.hpp"
#include "io/file_reader.hpp"
#include "io/read_processor.hpp"

#include <libcxx/sort.hpp>

using namespace hammer;

class BufferFiller;

class HammerKMerSplitter : public KMerSplitter<hammer::HKMer> {
  typedef std::vector<std::vector<HKMer> > KMerBuffer;

  void DumpBuffers(size_t num_files, size_t nthreads,
                   std::vector<KMerBuffer> &buffers,
                   const path::files_t &ostreams) const;

 public:
  HammerKMerSplitter(const std::string &work_dir)
      : KMerSplitter<hammer::HKMer>(work_dir, hammer::K) {}

  virtual path::files_t Split(size_t num_files);

  friend class BufferFiller;
};

void HammerKMerSplitter::DumpBuffers(size_t num_files, size_t nthreads,
                                     std::vector<KMerBuffer> &buffers,
                                     const path::files_t &ostreams) const {
# pragma omp parallel for num_threads(nthreads)
  for (unsigned k = 0; k < num_files; ++k) {
    size_t sz = 0;
    for (size_t i = 0; i < nthreads; ++i)
      sz += buffers[i][k].size();

    std::vector<HKMer> SortBuffer;
    SortBuffer.reserve(sz);
    for (size_t i = 0; i < nthreads; ++i) {
      KMerBuffer &entry = buffers[i];
      SortBuffer.insert(SortBuffer.end(), entry[k].begin(), entry[k].end());
    }
    libcxx::sort(SortBuffer.begin(), SortBuffer.end(), HKMer::less2_fast());
    auto it = std::unique(SortBuffer.begin(), SortBuffer.end());

#   pragma omp critical
    {
      FILE *f = fopen(ostreams[k].c_str(), "ab");
      VERIFY_MSG(f, "Cannot open temporary file to write");
      fwrite(SortBuffer.data(), sizeof(HKMer), it - SortBuffer.begin(), f);
      fclose(f);
    }
  }

  for (unsigned i = 0; i < nthreads; ++i) {
    for (unsigned j = 0; j < num_files; ++j) {
      buffers[i][j].clear();
    }
  }
}

class BufferFiller {
  std::vector<HammerKMerSplitter::KMerBuffer> &tmp_entries_;
  unsigned num_files_;
  size_t cell_size_;
  size_t processed_;
  const HammerKMerSplitter &splitter_;

 public:
  BufferFiller(std::vector<HammerKMerSplitter::KMerBuffer> &tmp_entries, size_t cell_size, const HammerKMerSplitter &splitter):
      tmp_entries_(tmp_entries), num_files_((unsigned)tmp_entries[0].size()), cell_size_(cell_size), processed_(0), splitter_(splitter) {}

  size_t processed() const { return processed_; }

  bool operator()(const io::SingleRead &r) {
    ValidHKMerGenerator<hammer::K> gen(r);
    HammerKMerSplitter::KMerBuffer &entry = tmp_entries_[omp_get_thread_num()];

#   pragma omp atomic
    processed_ += 1;

    bool stop = false;
    while (gen.HasMore()) {
      HKMer seq = gen.kmer(); size_t idx;

      idx = splitter_.GetFileNumForSeq(seq, num_files_);
      entry[idx].push_back(seq);
      stop |= entry[idx].size() > cell_size_;

      seq = !seq;

      idx = splitter_.GetFileNumForSeq(seq, num_files_);
      entry[idx].push_back(seq);
      stop |= entry[idx].size() > cell_size_;

      gen.Next();
    }

    return stop;
  }
};

path::files_t HammerKMerSplitter::Split(size_t num_files) {
  unsigned nthreads = cfg::get().max_nthreads;

  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  // Determine the set of output files
  path::files_t out;
  for (unsigned i = 0; i < num_files; ++i)
    out.push_back(GetRawKMersFname(i));

  size_t reads_buffer_size = cfg::get().count_split_buffer;
  if (reads_buffer_size == 0) {
    reads_buffer_size = 536870912ull;
    size_t mem_limit =  (size_t)((double)(get_free_memory()) / (nthreads * 3));
    INFO("Memory available for splitting buffers: " << (double)mem_limit / 1024.0 / 1024.0 / 1024.0 << " Gb");
    reads_buffer_size = std::min(reads_buffer_size, mem_limit);
  }
  size_t cell_size = reads_buffer_size / (num_files * sizeof(HKMer));
  // Set sane minimum cell size
  if (cell_size < 16384)
    cell_size = 16384;

  INFO("Using cell size of " << cell_size);
  std::vector<KMerBuffer> tmp_entries(nthreads);
  for (unsigned i = 0; i < nthreads; ++i) {
    KMerBuffer &entry = tmp_entries[i];
    entry.resize(num_files);
    for (unsigned j = 0; j < num_files; ++j) {
      entry[j].reserve((size_t)(1.1 * (double)cell_size));
    }
  }

  size_t n = 15;
  const auto& dataset = cfg::get().dataset;
  BufferFiller filler(tmp_entries, cell_size, *this);
  for (auto it = dataset.reads_begin(), et = dataset.reads_end(); it != et; ++it) {
    INFO("Processing " << *it);
    io::FileReadStream irs(*it, io::PhredOffset);
    hammer::ReadProcessor rp(nthreads);
    while (!irs.eof()) {
      rp.Run(irs, filler);
      DumpBuffers(num_files, nthreads, tmp_entries, out);
      VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");

      if (filler.processed() >> n) {
        INFO("Processed " << filler.processed() << " reads");
        n += 1;
      }
    }
  }
  INFO("Processed " << filler.processed() << " reads");

  return out;
}

static inline void Merge(KMerStat &lhs, const KMerStat &rhs) {
  if (lhs.count == 0)
    lhs.kmer = rhs.kmer;

  lhs.count += rhs.count;
  lhs.qual *= rhs.qual;
}

static void PushKMer(KMerData &data, HKMer kmer, double qual) {
  KMerStat &kmc = data[kmer];
  kmc.lock();
  Merge(kmc, KMerStat(1, kmer, qual));
  kmc.unlock();
}

static void PushKMerRC(KMerData &data, HKMer kmer, double qual) {
  kmer = !kmer;

  KMerStat &kmc = data[kmer];
  kmc.lock();
  Merge(kmc, KMerStat(1, kmer, qual));
  kmc.unlock();
}

class KMerDataFiller {
  KMerData &data_;

 public:
  KMerDataFiller(KMerData &data)
      : data_(data) {}

  bool operator()(const io::SingleRead &r) const {
    ValidHKMerGenerator<hammer::K> gen(r);
    while (gen.HasMore()) {
      HKMer kmer = gen.kmer();
      double correct = gen.correct_probability();

      PushKMer(data_, kmer, 1 - correct);
      PushKMerRC(data_, kmer, 1 - correct);

      gen.Next();
    }

    // Do not stop
    return false;
  }
};

void KMerDataCounter::FillKMerData(KMerData &data) {
  HammerKMerSplitter splitter(cfg::get().working_dir);
  KMerDiskCounter<hammer::HKMer> counter(cfg::get().working_dir, splitter);
  size_t sz = KMerIndexBuilder<HammerKMerIndex>(cfg::get().working_dir, num_files_, cfg::get().max_nthreads).BuildIndex(data.index_, counter);

  // Now use the index to fill the kmer quality information.
  INFO("Collecting K-mer information, this takes a while.");
  data.data_.resize(sz);

  const auto& dataset = cfg::get().dataset;
  for (auto it = dataset.reads_begin(), et = dataset.reads_end(); it != et; ++it) {
    INFO("Processing " << *it);
    io::FileReadStream irs(*it, io::PhredOffset);
    KMerDataFiller filler(data);
    hammer::ReadProcessor(cfg::get().max_nthreads).Run(irs, filler);
  }

  INFO("Collection done, postprocessing.");

  size_t singletons = 0;
  for (size_t i = 0; i < data.size(); ++i) {
    VERIFY(data[i].count);

    if (data[i].count == 1)
      singletons += 1;
  }

  INFO("Merge done. There are " << data.size() << " kmers in total. "
       "Among them " << singletons << " (" <<  100.0 * double(singletons) / double(data.size()) << "%) are singletons.");
}
