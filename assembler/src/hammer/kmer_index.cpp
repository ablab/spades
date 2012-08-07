//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "kmer_index.hpp"
#include "hammer_tools.hpp"
#include "valid_kmer_generator.hpp"
#include "mmapped_reader.hpp"
#include "mmapped_writer.hpp"

#include <algorithm>
#ifdef USE_GLIBCXX_PARALLEL
#include <parallel/algorithm>
#endif

#include "globals.hpp"

void KMerIndex::push_back(const KMerCount &k) {
  const char* s = Globals::blob + k.first.start();
  index_.insert(std::make_pair(Seq<K>(s, 0, K, /* raw */ true), data_.size()));
  data_.push_back(k);
}

KMerIndex &KMerIndex::operator+=(const KMerIndex &rhs) {
  for (auto I = rhs.seq_begin(), E = rhs.seq_end(); I != E; ++I) {
    const KMerCount &kmer = rhs[I->second];
    // Check, whether we already have this kmer.
    auto it = index_.find(I->first);
    if (it != index_.end()) {
      // If yes - merge it in.
      Merge(data_[it->second], kmer);
    } else {
      // Otherwise, just inser as-is.
      push_back(kmer);
    }
  }

  return *this;
}

static size_t my_hash(hint_t pos) {
  size_t hash = 877;
  for (size_t i = 0; i < K; i++) {
    hash = ((hash << 5) - hash) + ((int)Globals::blob[pos + i]) * 13;
  }
  return hash;
}

void KMerCounter::Split() {
  unsigned count_num_threads = min(cfg::get().count_merge_nthreads, cfg::get().general_max_nthreads);
  unsigned num_minimizers = 0;
  if (cfg::get().general_num_minimizers) {
    num_minimizers = *cfg::get().general_num_minimizers;
  }
  bool use_minimizers = HammerTools::doingMinimizers();
  int which_first = Globals::iteration_no % 4;

  INFO("Splitting kmer instances into " << num_files_ << " buckets. This might take a while");

  MMappedWriter* ostreams = new MMappedWriter[num_files_];
  for (unsigned i = 0; i < num_files_; ++i) {
    std::string filename = HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "tmp.kmers", i);
    ostreams[i].open(filename);
  }

  Seq<K>::hash hash_function;
  size_t readbuffer = cfg::get().count_split_buffer;
  std::vector< std::vector< std::vector< KMerNo > > > tmp_entries(count_num_threads);
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
          kmers.push_back(std::make_pair(cpos + gen.pos() - 1,
                                         std::make_pair(1 - gen.correct_probability(), my_hash(cpos + gen.pos() - 1))));
          gen.Next();
        }
        HammerTools::findMinimizers(kmers, num_minimizers, mmers, which_first);

        for (vector< pair<hint_t, pair< double, size_t > > >::const_iterator it = kmers.begin(); it != kmers.end(); ++it ) {
          tmp_entries[omp_get_thread_num()][it->second.second % num_files_].push_back(KMerNo(it->first, it->second.first));
        }
      } else {
        while (gen.HasMore()) {
          Seq<K> seq = gen.kmer();
          tmp_entries[omp_get_thread_num()][hash_function(seq) % num_files_].push_back(KMerNo(cpos + gen.pos() - 1,
                                                                                              1 - gen.correct_probability(),
                                                                                              seq));
          gen.Next();
        }
      }
    }
    cur_i = cur_limit;

    ++cur_fileindex;

    #pragma omp parallel for shared(tmp_entries) num_threads(count_num_threads)
    for (unsigned k = 0; k < num_files_; ++k) {
      size_t sz = 0;
      for (size_t i = 0; i < count_num_threads; ++i)
        sz += tmp_entries[i][k].size() * sizeof(tmp_entries[i][k][0]);

      ostreams[k].reserve(sz);
      for (size_t i = 0; i < count_num_threads; ++i) {
        ostreams[k].write(&tmp_entries[i][k][0], tmp_entries[i][k].size() * sizeof(tmp_entries[i][k][0]));
      }
    }

    for (unsigned i = 0; i < count_num_threads; ++i) {
      for (unsigned j = 0; j < num_files_; ++j) {
        tmp_entries[i][j].clear();
      }
    }
  }

  delete[] ostreams;
}

static void EquallySplit(size_t size, unsigned num_threads, size_t *borders) {
  size_t s = size / num_threads;
  size_t rem = size - s * num_threads;
  size_t current = 0;
  borders[0] = 0;
  for (unsigned i = 0; i < num_threads; i++) {
    if (rem > 0) {
      current += 1;
      rem -= 1;
    }
    current += s;
    borders[i + 1] = current;
  }
}

static void KmerHashUnique(const std::vector<KMerNo>::const_iterator first,
                           const std::vector<KMerNo>::const_iterator last,
                           std::vector<KMerCount> &result) {
  size_t size = last - first;
  if (size == 0)
    return;

  // counter    - number of unic kmers in block
  // borders    - borders of blocks
  // overlaps     - number of overlaping kmers
  // result_borders- borders of the result part for each thread
  // merged_overlaps-merge of overlap part for each thread

  size_t *counter, *borders, *overlaps, *result_borders;
  KMerCount *merged_overlaps;
  size_t total_unique = 0;
  unsigned num_threads = omp_get_max_threads();
# pragma omp parallel num_threads(num_threads)
  {
#   pragma omp single
    {
      num_threads = omp_get_num_threads();
      borders = new size_t[num_threads + 1];
      EquallySplit(size, num_threads, borders);
      counter = new size_t[num_threads];
      overlaps = new size_t[num_threads];
      merged_overlaps = new KMerCount[num_threads];
    }

    // Count unique and overlapping kmers for thread
    unsigned iam = omp_get_thread_num();

    size_t cnt = 0;
    size_t overlap_len = 0;
    if (iam == 0) {
      size_t begin = borders[0] + 1; // == 1
      size_t end = borders[iam + 1];

      overlap_len = 0;
      cnt = 1;
      auto I = first + begin, E = first + end;
      for (; I != E; ++I) {
        if (*I != *(I - 1))
          cnt += 1;
      }
    } else {
      size_t begin = borders[iam];
      size_t end = borders[iam + 1];
      // Count overlaps first (the ones which need to be merged into preceding
      // chunk).

      auto I = first + begin, E = first + end;
      for (; I != E; ++I) {
        if (*I == *(I - 1))
          overlap_len += 1;
        else
          break;
      }

      // Now count remaining entries
      I = first + begin + overlap_len, E = first + end;
      for (; I != E; ++I) {
        if (*I != *(I - 1))
          cnt += 1;
      }
    }
    overlaps[iam] = overlap_len;
    counter[iam] = cnt;

#   pragma omp barrier

#  pragma omp single
    {
      // Now we can deduce result_borders
      result_borders = new size_t[num_threads + 1];
      result_borders[0] = 0;
      for (unsigned i = 0; i < num_threads; i++) {
        total_unique += counter[i];
        result_borders[i + 1] = total_unique;
      }
      result.resize(total_unique);
    }

    // Merge overlaps to merged_overlaps if any
    size_t begin = borders[iam];
    size_t end = borders[iam] + overlaps[iam];

    if (overlaps[iam]) {
      hint_t bidx = (first + begin)->getIndex();
      const unsigned char *qdata = (const unsigned char*)(Globals::blobquality + bidx);
      merged_overlaps[iam] = KMerCount(PositionKMer(bidx),
                                       KMerStat(1, KMERSTAT_GOODITER,
                                                (first + begin)->getQual(), qdata));
      auto I = first + begin + 1, E = first + end;
      for (; I != E; ++I) {
        Merge(merged_overlaps[iam], *I);
      }
    }

    // Count unique kmers
    begin = borders[iam] + overlaps[iam];
    end = borders[iam + 1];
    size_t out_idx = result_borders[iam];

    // If block is one big overlap, do nothing
    if (begin != end) {
      hint_t bidx = (first + begin)->getIndex();
      const unsigned char *qdata = (const unsigned char*)(Globals::blobquality + bidx);
      result[out_idx] = KMerCount(PositionKMer(bidx),
                                  KMerStat(1, KMERSTAT_GOODITER, (first + begin)->getQual(), qdata));
      auto I = first + begin + 1, E = first + end;
      for (; I != E; ++I) {
        if (*I != *(I - 1)) {
          hint_t cidx = I->getIndex();
          const unsigned char *qdata = (const unsigned char*)(Globals::blobquality + cidx);
          result[++out_idx] = KMerCount(PositionKMer(cidx),
                                        KMerStat(1, KMERSTAT_GOODITER, I->getQual(), qdata));
        } else {
          Merge(result[out_idx], *I);
        }
      }
    }

#   pragma omp barrier
    // End of parallel part
  }

  // Merge merged overlaps to result
  for (size_t i = 0; i < num_threads; i++) {
    if (overlaps[i]) {
      Merge(result[result_borders[i] - 1], merged_overlaps[i]);
    }
  }

  delete[] counter;
  delete[] borders;
  delete[] overlaps;
  delete[] result_borders;
  delete[] merged_overlaps;
}

static size_t MergeKMers(const std::string &ifname, const std::string &ofname) {
  std::vector<KMerNo> vec;
  std::vector<KMerCount> vkmc;

  // Make sure memory mapping is released as soon as possible
  {
    MMappedRecordReader<KMerNo> ins(ifname, /* unlink */ true);

    vec.resize(ins.size());
    ins.read(&vec[0], vec.size());
    VERIFY(!ins.good());
  }

#ifdef USE_GLIBCXX_PARALLEL
  // Explicitly force a call to parallel sort routine.
  __gnu_parallel::sort(vec.begin(), vec.end(), KMerNo::is_less());
#else
  std::sort(vec.begin(), vec.end(), KMerNo::is_less());
#endif
  INFO("Sorting done, starting unification.");

  KmerHashUnique(vec.begin(), vec.end(), vkmc);

  MMappedWriter os(ofname);
  os.reserve(vkmc.size() * sizeof(vkmc[0]));
  os.write(&vkmc[0], vkmc.size() * sizeof(vkmc[0]));

  return vkmc.size();
}

void KMerCounter::BuildIndex(KMerIndex &index) {
  index.clear();

  // First, split k-mers into buckets.
  Split();

  unsigned count_num_threads = std::min(cfg::get().count_merge_nthreads, cfg::get().general_max_nthreads);
  INFO("Starting merge in " << count_num_threads << " threads.");

  size_t sz = 0;
  for (unsigned iFile = 0; iFile < num_files_; ++iFile) {
    INFO("Processing file " << iFile);
    std::string ifname = HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "tmp.kmers", iFile);
    std::string ofname = HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmer.idx", iFile);

    sz += MergeKMers(ifname, ofname);
  }

  INFO("Building kmer index");
  index.reserve(sz);
  for (unsigned iFile = 0; iFile < num_files_; ++iFile) {
    INFO("Processing file " << iFile);
    std::string ifname = HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmer.idx", iFile);
    MMappedRecordReader<KMerCount> ins(ifname, /* unlink */ true, -1ULL);

    index.push_back(ins.begin(), ins.end());
  }

  size_t singletons = 0;
  for (size_t i = 0; i < index.size(); ++i) {
    if (index[i].second.count == 1)
      singletons += 1;
  }

  INFO("Merge done. There are " << index.size() << " kmers in total. "
       "Among them " << singletons << " (" <<  100.0 * singletons / index.size() << "%) singletons.");
}
