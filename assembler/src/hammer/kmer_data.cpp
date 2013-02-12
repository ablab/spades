//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "kmer_data.hpp"
#include "io/read_processor.hpp"
#include "valid_kmer_generator.hpp"

#include "io/mmapped_writer.hpp"
#include "read/ireadstream.hpp"
#include "config_struct_hammer.hpp"

#include <libcxx/sort.hpp>

using namespace hammer;

class BufferFiller;

class HammerKMerSplitter : public KMerSplitter<hammer::KMer> {
  typedef std::vector<std::vector<KMer> > KMerBuffer;

  void DumpBuffers(size_t num_files, size_t nthreads,
                   std::vector<KMerBuffer> &buffers,
                   MMappedRecordWriter<KMer> *ostreams) const;

 public:
  HammerKMerSplitter(std::string &work_dir)
      : KMerSplitter<hammer::KMer>(work_dir, hammer::K) {}

  virtual path::files_t Split(size_t num_files);

  friend class BufferFiller;
};

void HammerKMerSplitter::DumpBuffers(size_t num_files, size_t nthreads,
                                     std::vector<KMerBuffer> &buffers,
                                     MMappedRecordWriter<KMer> *ostreams) const {
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
      size_t osz = it - SortBuffer.begin();
      ostreams[k].reserve(osz);
      ostreams[k].write(&SortBuffer[0], osz);
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
  size_t num_files_;
  size_t cell_size_;
  size_t processed_;
  const HammerKMerSplitter &splitter_;

 public:
  BufferFiller(std::vector<HammerKMerSplitter::KMerBuffer> &tmp_entries, size_t cell_size, const HammerKMerSplitter &splitter):
      tmp_entries_(tmp_entries), num_files_(tmp_entries[0].size()), cell_size_(cell_size), processed_(0), splitter_(splitter) {}

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

  MMappedRecordWriter<KMer>* ostreams = new MMappedRecordWriter<KMer>[num_files];
  for (unsigned i = 0; i < num_files; ++i)
    ostreams[i].open(out[i]);

  size_t read_buffer = cfg::get().count_split_buffer;
  size_t cell_size = (read_buffer / (nthreads * num_files * sizeof(KMer)));
  // Set sane minimum cell size
  if (cell_size < 16384)
    cell_size = 16384;

  INFO("Using cell size of " << cell_size);
  std::vector<KMerBuffer> tmp_entries(nthreads);
  for (unsigned i = 0; i < nthreads; ++i) {
    KMerBuffer &entry = tmp_entries[i];
    entry.resize(num_files);
    for (unsigned j = 0; j < num_files; ++j) {
      entry[j].reserve(1.1 * cell_size);
    }
  }

  size_t n = 15;
  BufferFiller filler(tmp_entries, cell_size, *this);
  for (auto I = Globals::input_filenames.begin(), E = Globals::input_filenames.end(); I != E; ++I) {
    ireadstream irs(*I, cfg::get().input_qvoffset);
    while (!irs.eof()) {
      hammer::ReadProcessor rp(nthreads);
      rp.Run(irs, filler);
      DumpBuffers(num_files, nthreads, tmp_entries, ostreams);
      VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");

      if (filler.processed() >> n) {
        INFO("Processed " << filler.processed() << " reads");
        n += 1;
      }
    }
  }
  INFO("Processed " << filler.processed() << " reads");

  delete[] ostreams;

  return out;
}

static inline void Merge(KMerStat &lhs, const KMerStat &rhs) {
  if (lhs.kmer() != rhs.kmer()) {
    lhs.kmer_ = rhs.kmer_;
  }

  lhs.count += rhs.count;
  lhs.totalQual *= rhs.totalQual;
  lhs.qual += rhs.qual;
}

static void PushKMer(KMerData &data,
                     KMer kmer, const unsigned char *q, double prob) {
  KMerStat &kmc = data[kmer];
  kmc.lock();
  Merge(kmc,
        KMerStat(1, kmer, prob, q));
  kmc.unlock();
}

static void PushKMerRC(KMerData &data,
                       KMer kmer, const unsigned char *q, double prob) {
  unsigned char rcq[K];

  // Prepare RC kmer with quality.
  kmer = !kmer;
  for (unsigned i = 0; i < K; ++i)
    rcq[K - i - 1] = q[i];

  KMerStat &kmc = data[kmer];
  kmc.lock();
  Merge(kmc,
        KMerStat(1, kmer, prob, rcq));
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

void KMerDataCounter::FillKMerData(KMerData &data) {
  // Build the index
  std::string workdir = cfg::get().input_working_dir;
  HammerKMerSplitter splitter(workdir);
  KMerDiskCounter<hammer::KMer> counter(workdir, splitter);
  size_t sz = KMerIndexBuilder<HammerKMerIndex>(workdir, num_files_, omp_get_max_threads()).BuildIndex(data.index_, counter);

  // Now use the index to fill the kmer quality information.
  INFO("Collecting K-mer information, this takes a while.");
  data.data_.resize(sz);

  KMerDataFiller filler(data);
  for (auto I = Globals::input_filenames.begin(), E = Globals::input_filenames.end(); I != E; ++I) {
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

  INFO("Merge done. There are " << data.size() << " kmers in total. "
       "Among them " << singletons << " (" <<  100.0 * singletons / data.size() << "%) are singletons.");
}
