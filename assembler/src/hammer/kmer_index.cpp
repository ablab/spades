//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "kmer_index.hpp"

#include "config_struct_hammer.hpp"
#include "hammer_tools.hpp"
#include "mmapped_reader.hpp"
#include "mmapped_writer.hpp"
#include "pointer_iterator.hpp"
#include "valid_kmer_generator.hpp"

#include <algorithm>
#ifdef USE_GLIBCXX_PARALLEL
#include <parallel/algorithm>
#endif

#include "globals.hpp"

std::string KMerIndexBuilder::GetRawKMersFname(unsigned suffix) const {
  // FIXME: This is ugly!
  std::ostringstream tmp;
  tmp.str("");
  tmp << work_dir_ << "/" << "kmers.raw." << suffix;

  return tmp.str();
}

std::string KMerIndexBuilder::GetUniqueKMersFname(unsigned suffix) const {
  // FIXME: This is ugly!
  std::ostringstream tmp;
  tmp.str("");
  tmp << work_dir_ << "/" << "kmers.unique." << suffix;

  return tmp.str();
}

void KMerIndexBuilder::Split(size_t num_files) {
  unsigned count_num_threads = min(cfg::get().count_merge_nthreads, cfg::get().general_max_nthreads);
  unsigned num_minimizers = 0;
  if (cfg::get().general_num_minimizers) {
    num_minimizers = *cfg::get().general_num_minimizers;
  }
  bool use_minimizers = HammerTools::doingMinimizers();
  int which_first = Globals::iteration_no % 4;

  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while");

  MMappedRecordWriter<KMer>* ostreams = new MMappedRecordWriter<KMer>[num_files];
  for (unsigned i = 0; i < num_files; ++i) {
    ostreams[i].open(GetRawKMersFname(i));
  }

  KMerIndex::hash_function hash_function;
  size_t readbuffer = cfg::get().count_split_buffer;
  std::vector<std::vector< std::vector<KMer> > > tmp_entries(count_num_threads);
  for (unsigned i = 0; i < count_num_threads; ++i) {
    tmp_entries[i].resize(num_files);
    for (unsigned j = 0; j < num_files; ++j) {
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
      if (use_minimizers) {
        // M is the larger k-mer size -- every m-mer should have several minimizer k-mers inside
        ValidKMerGenerator<M> gen_m(s, q, csize);
        std::vector<hint_t> mmers;
        while (gen_m.HasMore()) {
          mmers.push_back(cpos + gen_m.pos() - 1);
          gen_m.Next();
        }

        std::vector< std::pair<hint_t, std::pair< double, size_t > > > kmers;
        while (gen.HasMore()) {
          KMer seq = gen.kmer();
          kmers.push_back(std::make_pair(cpos + gen.pos() - 1,
                                         std::make_pair(1 - gen.correct_probability(), hash_function(seq))));
          gen.Next();
        }
        HammerTools::findMinimizers(kmers, num_minimizers, mmers, which_first);

        for (vector< pair<hint_t, pair< double, size_t > > >::const_iterator it = kmers.begin(); it != kmers.end(); ++it ) {
          const char *seq = Globals::blob + it->first;
          KMer kmer = KMer(seq, 0, K, /* raw */ true);
          tmp_entries[omp_get_thread_num()][it->second.second % num_files].push_back(kmer);
        }
      } else {
        while (gen.HasMore()) {
          KMer seq = gen.kmer();
          tmp_entries[omp_get_thread_num()][hash_function(seq) % num_files].push_back(seq);
          seq = !seq;
          tmp_entries[omp_get_thread_num()][hash_function(seq) % num_files].push_back(seq);
          gen.Next();
        }
      }
    }
    cur_i = cur_limit;

    ++cur_fileindex;

    #pragma omp parallel for shared(tmp_entries) num_threads(count_num_threads)
    for (unsigned k = 0; k < num_files; ++k) {
      size_t sz = 0;
      for (size_t i = 0; i < count_num_threads; ++i)
        sz += tmp_entries[i][k].size();

      ostreams[k].reserve(sz);
      for (size_t i = 0; i < count_num_threads; ++i) {
        ostreams[k].write(&tmp_entries[i][k][0], tmp_entries[i][k].size());
      }
    }

    for (unsigned i = 0; i < count_num_threads; ++i) {
      for (unsigned j = 0; j < num_files; ++j) {
        tmp_entries[i][j].clear();
      }
    }
  }

  delete[] ostreams;
}

size_t KMerIndexBuilder::MergeKMers(const std::string &ifname, const std::string &ofname) {
  MMappedRecordReader<KMer> ins(ifname, /* unlink */ true, -1ULL);
#ifdef USE_GLIBCXX_PARALLEL
  // Explicitly force a call to parallel sort routine.
  __gnu_parallel::sort(ins.begin(), ins.end(), KMer::less2());
#else
  std::sort(ins.begin(), ins.end(), KMer::less2());
#endif
  INFO("Sorting done, starting unification.");

  // FIXME: Use something like parallel version of unique_copy but with explicit
  // resizing.
  auto it = std::unique(ins.begin(), ins.end());

  MMappedRecordWriter<KMer> os(ofname);
  os.resize(it - ins.begin());
  std::copy(ins.begin(), it, os.begin());

  return it - ins.begin();
}

static inline void Merge(KMerStat &lhs, const KMerStat &rhs) {
  if (lhs.kmer() != rhs.kmer()) {
    lhs.kmer_ = rhs.kmer_;
  }

  lhs.count += rhs.count;
  lhs.totalQual *= rhs.totalQual;
  lhs.qual += rhs.qual;
}

size_t KMerIndexBuilder::BuildIndex(KMerIndex &index, size_t num_buckets) {
  index.clear();

  INFO("Building kmer index");

  // Split k-mers into buckets.
  Split(num_buckets);

  unsigned count_num_threads = std::min(cfg::get().count_merge_nthreads, cfg::get().general_max_nthreads);
  INFO("Starting k-mer counting in " << count_num_threads << " threads.");

  size_t sz = 0;
  index.bucket_starts_.push_back(sz);
  for (unsigned iFile = 0; iFile < num_buckets; ++iFile) {
    INFO("Processing file " << iFile);
    sz += MergeKMers(GetRawKMersFname(iFile), GetUniqueKMersFname(iFile));
    index.bucket_starts_.push_back(sz);
  }
  INFO("K-mer counting done. There are " << sz << " kmers in total. ");

  index.num_buckets_ = num_buckets;
  index.index_ = new KMerIndex::KMerDataIndex[num_buckets];
  index.bucket_locks_ = new omp_lock_t[num_buckets];

  for (size_t i = 0; i < num_buckets; ++i)
    omp_init_lock(index.bucket_locks_ + i);

  INFO("Building perfect hash indices");
# pragma omp parallel for shared(index)
  for (unsigned iFile = 0; iFile < num_buckets; ++iFile) {
    KMerIndex::KMerDataIndex &data_index = index.index_[iFile];
    MMappedRecordReader<KMer> ins(GetUniqueKMersFname(iFile), /* unlink */ true, -1ULL);
    if (!data_index.Reset(ins.begin(), ins.end(), ins.end() - ins.begin())) {
      INFO("Something went really wrong (read = this should not happen). Try to restart and see if the problem will be fixed.");
      exit(-1);
    }
  }

  INFO("Index built. Total " << index.mem_size() << " bytes occupied.");

  return sz;
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
  size_t sz = KMerIndexBuilder(cfg::get().input_working_dir).BuildIndex(data.index_, num_files_);

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
