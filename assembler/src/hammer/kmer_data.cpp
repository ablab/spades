//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "kmer_data.hpp"
#include "valid_kmer_generator.hpp"

#include "io/mmapped_writer.hpp"

#include "config_struct_hammer.hpp"

using namespace hammer;

class HammerKMerSplitter : public KMerSplitter<KMer> {
  typedef std::vector<std::vector<KMer> > KMerBuffer;

  void DumpBuffers(size_t num_files, size_t nthreads,
                   std::vector<KMerBuffer> &buffers,
                   MMappedRecordWriter<KMer> *ostreams) const;

 public:
  HammerKMerSplitter(std::string &work_dir)
      : KMerSplitter<KMer>(work_dir) {}

  virtual path::files_t Split(size_t num_files);
};

void HammerKMerSplitter::DumpBuffers(size_t num_files, size_t nthreads,
                                     std::vector<KMerBuffer> &buffers,
                                     MMappedRecordWriter<KMer> *ostreams) const {
# pragma omp parallel for num_threads(nthreads)
  for (unsigned k = 0; k < num_files; ++k) {
    size_t sz = 0;
    for (size_t i = 0; i < nthreads; ++i)
      sz += buffers[i][k].size();

    std::vector<KMer> SortBuffer;
    SortBuffer.reserve(sz);
    for (size_t i = 0; i < nthreads; ++i) {
      KMerBuffer &entry = buffers[i];
      SortBuffer.insert(SortBuffer.end(), entry[k].begin(), entry[k].end());
    }
    std::sort(SortBuffer.begin(), SortBuffer.end(), KMer::less2_fast());
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


path::files_t HammerKMerSplitter::Split(size_t num_files) {
  unsigned nthreads = min(cfg::get().count_merge_nthreads, cfg::get().general_max_nthreads);

  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  // Determine the set of output files
  path::files_t out;
  for (unsigned i = 0; i < num_files; ++i)
    out.push_back(GetRawKMersFname(i));

  MMappedRecordWriter<KMer>* ostreams = new MMappedRecordWriter<KMer>[num_files];
  for (unsigned i = 0; i < num_files; ++i)
    ostreams[i].open(out[i]);

  size_t read_buffer = cfg::get().count_split_buffer;
  size_t cell_size = (read_buffer * (Globals::read_length - K + 1)) /
                      (nthreads * num_files * sizeof(KMer));

  INFO("Using cell size of " << cell_size);
  std::vector<KMerBuffer> tmp_entries(nthreads);
  for (unsigned i = 0; i < nthreads; ++i) {
    KMerBuffer &entry = tmp_entries[i];
    entry.resize(num_files);
    for (unsigned j = 0; j < num_files; ++j) {
      entry[j].reserve(1.25 * cell_size);
    }
  }

  size_t cur_i = 0, cur_limit = 0, n = 15;
  size_t cur_fileindex = 0;
  while (cur_i < Globals::pr->size()) {
    cur_limit = std::min(cur_limit + read_buffer, Globals::pr->size());

    #pragma omp parallel for shared(tmp_entries) num_threads(nthreads)
    for (size_t i = cur_i; i < cur_limit; ++i) {
      const PositionRead &pr = Globals::pr->at(i);
      // Skip opaque reads
      if (!pr.valid())
        continue;

      size_t cpos = pr.start(), csize = pr.size();
      const char *s = Globals::blob + cpos;
      const char *q;
      std::string q_common;
      if (Globals::use_common_quality) {
        q_common.assign(csize, (char)Globals::common_quality);
        q = q_common.data();
      } else {
        q = Globals::blobquality + cpos;
      }

      ValidKMerGenerator<K> gen(s, q, csize);
      while (gen.HasMore()) {
        KMer seq = gen.kmer();
        tmp_entries[omp_get_thread_num()][GetFileNumForSeq(seq, num_files)].push_back(seq);
        seq = !seq;
        tmp_entries[omp_get_thread_num()][GetFileNumForSeq(seq, num_files)].push_back(seq);
        gen.Next();
      }
    }
    cur_i = cur_limit;

    ++cur_fileindex;

    if (cur_i >> n) {
      INFO("Processed " << cur_i << " reads");
      n += 1;
    }

    DumpBuffers(num_files, nthreads, tmp_entries, ostreams);
  }

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

void KMerCounter::FillKMerData(KMerData &data) {
  // Build the index
  std::string workdir = cfg::get().input_working_dir;
  HammerKMerSplitter splitter(workdir);
  size_t sz = KMerIndexBuilder<KMer>(workdir, num_files_, omp_get_max_threads()).BuildIndex(data.index_, splitter);

  // Now use the index to fill the kmer quality information.
  INFO("Collecting K-mer information, this takes a while.");
  data.data_.resize(sz);

# pragma omp parallel for shared(data)
  for (size_t readno = 0; readno < Globals::pr->size(); ++readno) {
    PositionRead &pr = Globals::pr->at(readno);

    // skip opaque reads w/o kmers
    if (!pr.valid()) continue;

    size_t cpos = pr.start(), csize = pr.size();
    const char *s = Globals::blob + cpos;
    const char *q;
    std::string q_common;
    if (Globals::use_common_quality) {
      q_common.assign(csize, (char)Globals::common_quality);
      q = q_common.data();
    } else {
      q = Globals::blobquality + cpos;
    }

    ValidKMerGenerator<K> gen(s, q, csize);
    while (gen.HasMore()) {
      const KMer &kmer = gen.kmer();
      const unsigned char *kq = (const unsigned char*)(q + gen.pos() - 1);

      PushKMer(data, kmer, kq, 1 - gen.correct_probability());
      PushKMerRC(data, kmer, kq, 1 - gen.correct_probability());

      gen.Next();
    }
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
