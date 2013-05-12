
//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "adapter_index.hpp"
#include "io/read_processor.hpp"
#include "valid_kmer_generator.hpp"

#include "io/mmapped_writer.hpp"
#include "io/ireadstream.hpp"
#include "config_struct_cclean.hpp"

#include <libcxx/sort.hpp>

using namespace cclean;

class BufferFiller;

class AdapterKMerSplitter : public KMerSplitter<cclean::KMer> {
  typedef std::vector<std::vector<KMer> > KMerBuffer;

  void DumpBuffers(size_t num_files, size_t nthreads,
                   std::vector<KMerBuffer> &buffers,
                   MMappedRecordWriter<KMer> *ostreams) const;

 public:
  AdapterKMerSplitter(std::string &work_dir)
      : KMerSplitter<cclean::KMer>(work_dir, cclean::K) {}

  virtual path::files_t Split(size_t num_files);

  friend class BufferFiller;
};

void AdapterKMerSplitter::DumpBuffers(size_t num_files, size_t nthreads,
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
  std::vector<AdapterKMerSplitter::KMerBuffer> &tmp_entries_;
  size_t num_files_;
  size_t cell_size_;
  size_t processed_;
  const AdapterKMerSplitter &splitter_;

 public:
  BufferFiller(std::vector<AdapterKMerSplitter::KMerBuffer> &tmp_entries, size_t cell_size, const AdapterKMerSplitter &splitter):
      tmp_entries_(tmp_entries), num_files_(tmp_entries[0].size()), cell_size_(cell_size), processed_(0), splitter_(splitter) {}

  size_t processed() const { return processed_; }

  bool operator()(const Read &r) {
    AdapterKMerSplitter::KMerBuffer &entry = tmp_entries_[omp_get_thread_num()];
    const std::string &seq = r.getSequenceString();
    ValidKMerGenerator<cclean::K> gen(seq.c_str(), NULL, seq.size());
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

#   pragma omp atomic
    processed_ += 1;

    return stop;
  }
};

path::files_t AdapterKMerSplitter::Split(size_t num_files) {
  unsigned nthreads = cfg::get().nthreads;

  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  // Determine the set of output files
  path::files_t out;
  for (unsigned i = 0; i < num_files; ++i)
    out.push_back(GetRawKMersFname(i));

  MMappedRecordWriter<KMer>* ostreams = new MMappedRecordWriter<KMer>[num_files];
  for (unsigned i = 0; i < num_files; ++i)
    ostreams[i].open(out[i]);

  size_t read_buffer = cfg::get().count_split_buffer;
  size_t cell_size = (read_buffer / (num_files * sizeof(KMer)));
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
  // for (auto I = Globals::input_filenames.begin(), E = Globals::input_filenames.end(); I != E; ++I)
  {
    //ireadstream irs(*I);
    ireadstream irs("adapters.fasta.gz");
    while (!irs.eof()) {
      hammer::ReadProcessor rp(nthreads);
      rp.Run(irs, filler);
      DumpBuffers(num_files, nthreads, tmp_entries, ostreams);
      VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");

      if (filler.processed() >> n) {
        INFO("Processed " << filler.processed() << " sequences");
        n += 1;
      }
    }
  }
  INFO("Processed " << filler.processed() << " sequences");

  delete[] ostreams;

  return out;
}

class AdapterIndexFiller {
  AdapterIndex &data_;

 public:
  AdapterIndexFiller(AdapterIndex &data)
      : data_(data) {}

  bool operator()(const Read &r) {
    const std::string &seq = r.getSequenceString();
    ValidKMerGenerator<cclean::K> gen(seq.c_str(), NULL, seq.size());
    while (gen.HasMore()) {
      KMer kmer = gen.kmer();

      auto& entry = data_[kmer];
      entry.lock();
      entry.kmer_ = kmer;
      entry.seqs_.insert(seq);
      entry.unlock();

      kmer = !kmer;
      auto& rcentry = data_[kmer];
      rcentry.lock();
      rcentry.kmer_ = kmer;
      rcentry.seqs_.insert(ReverseComplement(seq));
      rcentry.unlock();

      gen.Next();
    }

    return false;
  }
};


void AdapterIndexBuilder::FillAdapterIndex(AdapterIndex &data) {
  // Build the index
  std::string workdir = cfg::get().input_working_dir;
  AdapterKMerSplitter splitter(workdir);
  KMerDiskCounter<cclean::KMer> counter(workdir, splitter);
  size_t sz = KMerIndexBuilder<Index>(workdir, num_files_, omp_get_max_threads()).BuildIndex(data.index_, counter);

  // Now use the index to fill the kmer quality information.
  INFO("Collecting K-mer information, this takes a while.");
  data.data_.resize(sz);

  AdapterIndexFiller filler(data);
  // for (auto I = Globals::input_filenames.begin(), E = Globals::input_filenames.end(); I != E; ++I)
  {
    //ireadstream irs(*I);
    ireadstream irs("adapters.fasta.gz");
    hammer::ReadProcessor rp(omp_get_max_threads());
    rp.Run(irs, filler);
    VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
  }
}
