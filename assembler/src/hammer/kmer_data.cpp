//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "kmer_data.hpp"
#include "mmapped_writer.hpp"
#include "valid_kmer_generator.hpp"

#include "config_struct_hammer.hpp"

class HammerKMerSplitter : public KMerSplitter {
 public:
  HammerKMerSplitter(std::string &work_dir, unsigned num_files)
      : KMerSplitter(work_dir, num_files) {}

  virtual path::files_t Split();
};

path::files_t HammerKMerSplitter::Split() {
  unsigned count_num_threads = min(cfg::get().count_merge_nthreads, cfg::get().general_max_nthreads);

  INFO("Splitting kmer instances into " << num_files_ << " buckets. This might take a while.");

  // Determine the set of output files
  path::files_t out;
  for (unsigned i = 0; i < num_files_; ++i)
    out.push_back(GetRawKMersFname(i));

  MMappedRecordWriter<KMer>* ostreams = new MMappedRecordWriter<KMer>[num_files_];
  for (unsigned i = 0; i < num_files_; ++i)
    ostreams[i].open(out[i]);

  KMerIndex::hash_function hash_function;
  size_t readbuffer = cfg::get().count_split_buffer;
  std::vector<std::vector< std::vector<KMer> > > tmp_entries(count_num_threads);
  for (unsigned i = 0; i < count_num_threads; ++i) {
    tmp_entries[i].resize(num_files_);
    for (unsigned j = 0; j < num_files_; ++j) {
      tmp_entries[i][j].reserve(1.25 * readbuffer / count_num_threads);
    }
  }

  size_t cur_i = 0, cur_limit = 0;
  size_t cur_fileindex = 0;
  while (cur_i < Globals::pr->size()) {
    cur_limit = min(cur_limit + readbuffer, Globals::pr->size());

    #pragma omp parallel for shared(tmp_entries) num_threads(count_num_threads)
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
        tmp_entries[omp_get_thread_num()][hash_function(seq) % num_files_].push_back(seq);
        seq = !seq;
        tmp_entries[omp_get_thread_num()][hash_function(seq) % num_files_].push_back(seq);
        gen.Next();
      }
    }
    cur_i = cur_limit;

    ++cur_fileindex;

    for (unsigned k = 0; k < num_files_; ++k) {
      size_t sz = 0;
      for (size_t i = 0; i < count_num_threads; ++i)
        sz += tmp_entries[i][k].size();

      ostreams[k].reserve(sz);
      for (size_t i = 0; i < count_num_threads; ++i) {
        ostreams[k].write(&tmp_entries[i][k][0], tmp_entries[i][k].size());
      }
    }

    for (unsigned i = 0; i < count_num_threads; ++i) {
      for (unsigned j = 0; j < num_files_; ++j) {
        tmp_entries[i][j].clear();
      }
    }
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
  HammerKMerSplitter splitter(workdir, num_files_);
  size_t sz = KMerIndexBuilder(workdir).BuildIndex(data.index_, splitter);

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
