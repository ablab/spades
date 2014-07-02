//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "kmer_data.hpp"
#include "io/read_processor.hpp"
#include "valid_kmer_generator.hpp"

#include "io/mmapped_writer.hpp"
#include "io/ireadstream.hpp"
#include "config_struct_hammer.hpp"

#include "file_limit.hpp"

#include <libcxx/sort.hpp>

using namespace hammer;

class BufferFiller;

class HammerKMerSplitter : public KMerSplitter<hammer::KMer> {
  typedef std::vector<std::vector<KMer> > KMerBuffer;

  void DumpBuffers(size_t num_files, size_t nthreads,
                   std::vector<KMerBuffer> &buffers,
                   const path::files_t &ostreams) const;

 public:
  HammerKMerSplitter(std::string &work_dir)
      : KMerSplitter<hammer::KMer>(work_dir, hammer::K) {}

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

    if (!sz)
      continue;

    std::vector<KMer> SortBuffer;
    SortBuffer.reserve(sz);
    for (size_t i = 0; i < nthreads; ++i) {
      KMerBuffer &entry = buffers[i];
      SortBuffer.insert(SortBuffer.end(), entry[k].begin(), entry[k].end());
    }
    libcxx::sort(SortBuffer.begin(), SortBuffer.end(), KMer::less2_fast());
    auto it = std::unique(SortBuffer.begin(), SortBuffer.end());

#   pragma omp critical
    {
        FILE *f = fopen(ostreams[k].c_str(), "ab");
        VERIFY_MSG(f, "Cannot open temporary file to write");
        fwrite(SortBuffer.data(), sizeof(KMer), it - SortBuffer.begin(), f);
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

  bool operator()(const Read &r) {
    int trim_quality = cfg::get().input_trim_quality;

    // FIXME: Get rid of this
    Read cr = r;
    size_t sz = cr.trimNsAndBadQuality(trim_quality);

    #pragma omp atomic
    processed_ += 1;

    if (sz < hammer::K)
      return false;

    HammerKMerSplitter::KMerBuffer &entry = tmp_entries_[omp_get_thread_num()];
    ValidKMerGenerator<hammer::K> gen(cr);
    bool stop = false;
    while (gen.HasMore()) {
      KMer seq = gen.kmer();
      size_t idx = splitter_.GetFileNumForSeq(seq, num_files_);
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
  unsigned nthreads = std::min(cfg::get().count_merge_nthreads, cfg::get().general_max_nthreads);

  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  // Determine the set of output files
  path::files_t out;
  for (unsigned i = 0; i < num_files; ++i)
    out.push_back(GetRawKMersFname(i));

  size_t file_limit = num_files + 2*nthreads;
  size_t res = limit_file(file_limit);
  if (res < file_limit) {
    WARN("Failed to setup necessary limit for number of open files. The process might crash later on.");
    WARN("Do 'ulimit -n " << file_limit << "' in the console to overcome the limit");
  }

  size_t reads_buffer_size = cfg::get().count_split_buffer;
  if (reads_buffer_size == 0) {
    reads_buffer_size = 536870912ull;
    size_t mem_limit =  (size_t)((double)(get_free_memory()) / (nthreads * 3));
    INFO("Memory available for splitting buffers: " << (double)mem_limit / 1024.0 / 1024.0 / 1024.0 << " Gb");
    reads_buffer_size = std::min(reads_buffer_size, mem_limit);
  }
  size_t cell_size = reads_buffer_size / (num_files * sizeof(KMer));
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
  BufferFiller filler(tmp_entries, cell_size, *this);
  const auto& dataset = cfg::get().dataset;
  for (auto I = dataset.reads_begin(), E = dataset.reads_end(); I != E; ++I) {
    INFO("Processing " << *I);
    ireadstream irs(*I, cfg::get().input_qvoffset);
    while (!irs.eof()) {
      hammer::ReadProcessor rp(nthreads);
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
  lhs.count += rhs.count;
  lhs.totalQual *= rhs.totalQual;
  lhs.qual += rhs.qual;
}

static void PushKMer(KMerData &data,
                     KMer kmer, const unsigned char *q, double prob) {
  size_t idx = data.seq_idx(kmer);
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

  size_t idx = data.seq_idx(kmer);
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

  bool operator()(const Read &r) {
    int trim_quality = cfg::get().input_trim_quality;

    // FIXME: Get rid of this
    Read cr = r;
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

void KMerDataCounter::BuildKMerIndex(KMerData &data) {
  // Build the index
  std::string workdir = cfg::get().input_working_dir;
  HammerKMerSplitter splitter(workdir);
  KMerDiskCounter<hammer::KMer> counter(workdir, splitter);
  size_t kmers = KMerIndexBuilder<HammerKMerIndex>(workdir, num_files_, omp_get_max_threads()).BuildIndex(data.index_, counter, /* save finall */ true);

  // Check, whether we'll ever have enough memory for running BH and bail out earlier
  if (1.25 * kmers * (sizeof(KMerStat) + sizeof(hammer::KMer)) > (double) get_memory_limit())
      FATAL_ERROR("The reads contain too many k-mers to fit into available memory limit. Increase memory limit and restart")

  data.kmers_ = counter.GetFinalKMers();

  size_t swaps = 0;
  INFO("Arranging kmers in hash map order");
  for (auto I = data.kmers_->begin(), E = data.kmers_->end(); I != E; ++I) {
    size_t cidx = I - data.kmers_->begin();
    size_t kidx = data.index_.raw_seq_idx(*I);
    while (cidx != kidx) {
      auto J = data.kmers_->begin() + kidx;
      using std::swap;
      swap(*I, *J);
      swaps += 1;

      kidx = data.index_.raw_seq_idx(*I);
    }
  }
  INFO("Done. Total swaps: " << swaps);
}

void KMerDataCounter::FillKMerData(KMerData &data) {
  // Now use the index to fill the kmer quality information.
  INFO("Collecting K-mer information, this takes a while.");
  data.data_.resize(data.kmers_->size());

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
    VERIFY(data[i].count);
    // Make sure all the kmers are marked as 'Bad' in the beginning
    data[i].status = KMerStat::Bad;

    if (data[i].count == 1)
      singletons += 1;
  }

  INFO("There are " << data.size() << " kmers in total. "
       "Among them " << singletons << " (" <<  100.0 * (double)singletons / (double)data.size() << "%) are singletons.");
}
