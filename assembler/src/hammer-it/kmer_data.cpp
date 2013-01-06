//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "kmer_data.hpp"
#include "valid_hkmer_generator.hpp"

#include "io/mmapped_writer.hpp"

#include "read/read.hpp"
#include "read/ireadstream.hpp"

#include <libcxx/sort.hpp>

using namespace hammer;

class HammerKMerSplitter : public KMerSplitter<hammer::HKMer> {
  typedef std::vector<std::vector<HKMer> > KMerBuffer;

  void DumpBuffers(size_t num_files, size_t nthreads,
                   std::vector<KMerBuffer> &buffers,
                   MMappedRecordWriter<HKMer> *ostreams) const;

 public:
  HammerKMerSplitter(std::string &work_dir)
      : KMerSplitter<hammer::HKMer>(work_dir, hammer::K) {}

  virtual path::files_t Split(size_t num_files);
};

void HammerKMerSplitter::DumpBuffers(size_t num_files, size_t nthreads,
                                     std::vector<KMerBuffer> &buffers,
                                     MMappedRecordWriter<HKMer> *ostreams) const {
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
  unsigned nthreads = 1; // min(cfg::get().count_merge_nthreads, cfg::get().general_max_nthreads);

  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  // Determine the set of output files
  path::files_t out;
  for (unsigned i = 0; i < num_files; ++i)
    out.push_back(GetRawKMersFname(i));

  MMappedRecordWriter<HKMer>* ostreams = new MMappedRecordWriter<HKMer>[num_files];
  for (unsigned i = 0; i < num_files; ++i)
    ostreams[i].open(out[i]);

  size_t read_buffer = 0; // cfg::get().count_split_buffer;
  size_t cell_size = (read_buffer * (400 - hammer::K + 1)) /
                      (nthreads * num_files * sizeof(HKMer));
  // Set sane minimum cell size
  if (cell_size < 16384)
    cell_size = 1638400;

  INFO("Using cell size of " << cell_size);
  std::vector<KMerBuffer> tmp_entries(nthreads);
  for (unsigned i = 0; i < nthreads; ++i) {
    KMerBuffer &entry = tmp_entries[i];
    entry.resize(num_files);
    for (unsigned j = 0; j < num_files; ++j) {
      entry[j].reserve(1.25 * cell_size);
    }
  }

  ireadstream irs("test.fastq", 33);
  while (!irs.eof()) {
    Read r;
    irs >> r;

    ValidHKMerGenerator<hammer::K> gen(r);
    while (gen.HasMore()) {
      HKMer seq = gen.kmer();
      tmp_entries[omp_get_thread_num()][GetFileNumForSeq(seq, num_files)].push_back(seq);
      seq = !seq;
      tmp_entries[omp_get_thread_num()][GetFileNumForSeq(seq, num_files)].push_back(seq);
      gen.Next();
    }
  }

  DumpBuffers(num_files, nthreads, tmp_entries, ostreams);

  delete[] ostreams;

  return out;
}

void KMerDataCounter::FillKMerData(KMerData &data) {
  std::string workdir(".");

  HammerKMerSplitter splitter(workdir);
  KMerDiskCounter<hammer::HKMer> counter(workdir, splitter);
  size_t sz = KMerIndexBuilder<HammerKMerIndex>(workdir, num_files_, omp_get_max_threads()).BuildIndex(data.index_, counter);

  // Now use the index to fill the kmer quality information.
  INFO("Collecting K-mer information, this takes a while.");
  data.data_.resize(sz);

  ireadstream irs("test.fastq", 33);
  while (!irs.eof()) {
    Read r;
    irs >> r;

    ValidHKMerGenerator<hammer::K> gen(r);
    while (gen.HasMore()) {
      HKMer seq = gen.kmer();
      data[seq].count += 1;
      data[!seq].count += 1;
      
      gen.Next();
    }
  }

    INFO("Collection done, postprocessing.");

  size_t singletons = 0;
  for (size_t i = 0; i < data.size(); ++i) {
    VERIFY(data[i].count);

    if (data[i].count == 1)
      singletons += 1;
  }
  
  INFO("Merge done. There are " << data.size() << " kmers in total. "
       "Among them " << singletons << " (" <<  100.0 * singletons / data.size() << "%) are singletons.");
}

