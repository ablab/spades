//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "kmer_index.hpp"

#include "mmapped_reader.hpp"
#include "mmapped_writer.hpp"
#include "pointer_iterator.hpp"

#include "path_helper.hpp"

#include <boost/lexical_cast.hpp>

#include <algorithm>
#ifdef USE_GLIBCXX_PARALLEL
#include <parallel/algorithm>
#endif

#include "globals.hpp"

std::string KMerSplitter::GetRawKMersFname(unsigned suffix) const {
  return path::append_path(work_dir_, "kmers.raw." + boost::lexical_cast<std::string>(suffix));
}

std::string KMerIndexBuilder::GetUniqueKMersFname(unsigned suffix) const {
  return path::append_path(work_dir_, "kmers.unique." + boost::lexical_cast<std::string>(suffix));
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

size_t KMerIndexBuilder::BuildIndex(KMerIndex &index, KMerSplitter &splitter) {
  index.clear();

  INFO("Building kmer index");

  // Split k-mers into buckets.
  path::files_t buckets = splitter.Split();
  unsigned num_buckets = buckets.size();

  INFO("Starting k-mer counting.");
  size_t sz = 0;
  index.bucket_starts_.push_back(sz);
  for (unsigned iFile = 0; iFile < num_buckets; ++iFile) {
    INFO("Processing file " << iFile);
    sz += MergeKMers(buckets[iFile], GetUniqueKMersFname(iFile));
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
