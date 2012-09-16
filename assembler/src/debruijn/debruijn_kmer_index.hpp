//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef DEBRUIJN_KMER_INDEX_HPP
#define DEBRUIJN_KMER_INDEX_HPP

#include "openmp_wrapper.h"

#include "io/multifile_reader.hpp"
#include "read_converter.hpp"

#include "mph_index/kmer_index.hpp"

namespace debruijn_graph {

// used for temporary reads storage during parallel reading
static const size_t READS_BUFFER_SIZE = 536870912; // 512 MB in bytes

typedef ::KMerSplitter<runtime_k::RtSeq> RtSeqKMerSplitter;

template<class Read>
class DeBruijnKMerSplitter : public RtSeqKMerSplitter {
  unsigned K_;
  io::ReadStreamVector<io::IReader<Read>> &streams_;
  SingleReadStream *contigs_;

 private:
  typedef std::vector< std::vector<runtime_k::RtSeq> > KMerBuffer;

  template<class ReadStream>
  std::pair<size_t, size_t>
  FillBufferFromStream(ReadStream& stream,
                       KMerBuffer &tmp_entries,
                       unsigned num_files, size_t cell_size) const;
  void DumpBuffers(size_t num_files, size_t nthreads,
                   std::vector<KMerBuffer> &buffers,
                   FILE **ostreams) const;

 public:
  DeBruijnKMerSplitter(const std::string &work_dir,
                       unsigned K,
                       io::ReadStreamVector< io::IReader<Read> >& streams,
                       SingleReadStream* contigs_stream = 0)
      : RtSeqKMerSplitter(work_dir), K_(K), streams_(streams), contigs_(contigs_stream) {
  }

  virtual path::files_t Split(size_t num_files);
};

template<class Read> template<class ReadStream>
std::pair<size_t, size_t>
DeBruijnKMerSplitter<Read>::FillBufferFromStream(ReadStream &stream,
                                                 KMerBuffer &buffer,
                                                 unsigned num_files, size_t cell_size) const {
  typename ReadStream::read_type r;
  unsigned k_plus_one = K_ + 1;
  size_t reads = 0, kmers = 0, rl = 0;

  while (!stream.eof() && kmers < num_files * cell_size) {
    stream >> r;
    rl = std::max(rl, r.size());
    reads += 1;

    const Sequence &seq = r.sequence();
    if (seq.size() < k_plus_one)
      continue;

    runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq::max_size>(k_plus_one);
    buffer[this->GetFileNumForSeq(kmer, num_files)].push_back(kmer);
    for (size_t j = k_plus_one; j < seq.size(); ++j) {
      kmer <<= seq[j];
      buffer[this->GetFileNumForSeq(kmer, num_files)].push_back(kmer);
    }

    kmers += ((seq.size() - k_plus_one + 1) + 1);
  }

  return std::make_pair(reads, rl);
}

template<class Read>
void DeBruijnKMerSplitter<Read>::DumpBuffers(size_t num_files, size_t nthreads,
                                             std::vector<KMerBuffer> &buffers,
                                             FILE **ostreams) const {
  unsigned k_plus_one = K_ + 1;
  size_t item_size = sizeof(runtime_k::RtSeq::DataType),
             items = runtime_k::RtSeq::GetDataSize(k_plus_one);

# pragma omp parallel for shared(items, item_size, k_plus_one)
  for (unsigned k = 0; k < num_files; ++k) {
    size_t sz = 0;
    for (size_t i = 0; i < nthreads; ++i)
      sz += buffers[i][k].size();

    std::vector<runtime_k::RtSeq> SortBuffer;
    SortBuffer.reserve(sz);
    for (size_t i = 0; i < nthreads; ++i) {
      KMerBuffer &entry = buffers[i];
      for (size_t j = 0; j < entry[k].size(); ++j)
        SortBuffer.push_back(entry[k][j]);
    }
    std::sort(SortBuffer.begin(), SortBuffer.end(), runtime_k::RtSeq::less2_fast());
    auto it = std::unique(SortBuffer.begin(), SortBuffer.end());

#   pragma omp critical
    {
      for (auto I = SortBuffer.begin(), E = it; I != E; ++I)
        fwrite(I->data(), item_size, items, ostreams[k]);
    }
  }

  for (unsigned i = 0; i < nthreads; ++i) {
    for (unsigned j = 0; j < num_files; ++j) {
      buffers[i][j].clear();
    }
  }
}

template<class Read>
path::files_t DeBruijnKMerSplitter<Read>::Split(size_t num_files) {
  unsigned nthreads = streams_.size();

  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  // Determine the set of output files
  path::files_t out;
  for (unsigned i = 0; i < num_files; ++i)
    out.push_back(this->GetRawKMersFname(i));

  FILE** ostreams = new FILE*[num_files];
  for (unsigned i = 0; i < num_files; ++i)
    ostreams[i] = fopen(out[i].c_str(), "wb");

  size_t cell_size = READS_BUFFER_SIZE /
                     (nthreads * num_files * sizeof(runtime_k::RtSeq));
  INFO("Using cell size of " << cell_size);

  std::vector<KMerBuffer> tmp_entries(nthreads);
  for (unsigned i = 0; i < nthreads; ++i) {
    KMerBuffer &entry = tmp_entries[i];
    entry.resize(num_files);
    for (unsigned j = 0; j < num_files; ++j) {
      entry[j].reserve(1.25 * cell_size);
    }
  }

  size_t counter = 0, rl = 0;
  streams_.reset();
  while (!streams_.eof()) {
#   pragma omp parallel for num_threads(nthreads) reduction(+ : counter) shared(rl)
    for (size_t i = 0; i < nthreads; ++i) {
      std::pair<size_t, size_t> stats = FillBufferFromStream(streams_[i], tmp_entries[i], num_files, cell_size);
      counter += stats.first;

      // There is no max reduction in C/C++ OpenMP... Only in FORTRAN :(
#     pragma omp flush(rl)
      if (stats.second > rl)
#     pragma omp critical
      {
        rl = std::max(rl, stats.second);
      }
    }

    DumpBuffers(num_files, nthreads, tmp_entries, ostreams);
  }

  if (contigs_) {
    INFO("Adding contigs from previous K");
    size_t cnt = 0;
    contigs_->reset();
    while (!contigs_->eof()) {
      FillBufferFromStream(*contigs_, tmp_entries[cnt], num_files, cell_size);
      DumpBuffers(num_files, nthreads, tmp_entries, ostreams);
      if (++cnt >= nthreads)
        cnt = 0;
    }
  }

  for (unsigned i = 0; i < num_files; ++i)
    fclose(ostreams[i]);

  delete[] ostreams;

  INFO("Used " << counter << " reads. Maximum read length " << rl);

  return out;
}

template<class IdType>
class DeBruijnKMerIndex {
  unsigned K_;
  KMerIndex<runtime_k::RtSeq> index_;
  std::vector<EdgeInfo<IdType> > data_;

 public:
  DeBruijnKMerIndex(unsigned K) 
      : K_(K), index_(K + 1) {}

  unsigned K() const { return K_; }
  
  template<class Read>
  friend class DeBruijnKMerIndexBuilder;
};

template<class Read>
class DeBruijnKMerIndexBuilder {
  template<class ReadStream, class IdType>
  size_t FillCoverageFromStream(ReadStream &stream,
                                DeBruijnKMerIndex<IdType> &index) const;

 public:
  template<class IdType>
  size_t BuildIndex(DeBruijnKMerIndex<IdType> &index, unsigned k,
                    io::ReadStreamVector<io::IReader<Read> > &streams,
                    SingleReadStream* contigs_stream = 0);
};

template<class Read> template<class ReadStream, class IdType>
size_t
DeBruijnKMerIndexBuilder<Read>::FillCoverageFromStream(ReadStream &stream,
                                                       DeBruijnKMerIndex<IdType> &index) const {
  typename ReadStream::read_type r;
  unsigned k_plus_one = index.K() + 1;
  size_t rl = 0;

  while (!stream.eof()) {
    stream >> r;
    rl = std::max(rl, r.size());

    const Sequence &seq = r.sequence();
    if (seq.size() < k_plus_one)
      continue;

    runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq::max_size>(k_plus_one);
    EdgeInfo<IdType> &entry = index.data_[index.index_.seq_idx(kmer)];
#   pragma omp atomic
    entry.count_ += 1;
    for (size_t j = k_plus_one; j < seq.size(); ++j) {
      kmer <<= seq[j];
      EdgeInfo<IdType> &entry = index.data_[index.index_.seq_idx(kmer)];
#     pragma omp atomic
      entry.count_ += 1;
    }
  }

  return rl;
}

template<class Read> template<class IdType>
size_t
DeBruijnKMerIndexBuilder<Read>::BuildIndex(DeBruijnKMerIndex<IdType> &index, unsigned k,
                                           io::ReadStreamVector<io::IReader<Read> > &streams,
                                           SingleReadStream* contigs_stream) {
  unsigned nthreads = streams.size();
  DeBruijnKMerSplitter<Read> splitter(cfg::get().output_dir, k, streams, contigs_stream);
  KMerIndexBuilder<runtime_k::RtSeq> builder(cfg::get().output_dir, 16, streams.size());
  size_t sz = builder.BuildIndex(index.index_, splitter);

  // Now use the index to fill the coverage and EdgeId's
  INFO("Collecting edges information, this takes a while.");
  index.data_.resize(sz);

  size_t rl = 0;
  streams.reset();
# pragma omp parallel for num_threads(nthreads) shared(rl)
  for (size_t i = 0; i < nthreads; ++i) {
    size_t crl = FillCoverageFromStream(streams[i], index);

    // There is no max reduction in C/C++ OpenMP... Only in FORTRAN :(
#   pragma omp flush(rl)
    if (crl > rl)
#     pragma omp critical
    {
      rl = std::max(rl, crl);
    }
  }

  if (contigs_stream) {
    contigs_stream->reset();
    FillCoverageFromStream(*contigs_stream, index);
  }

  return rl;
}
 
}

#endif // DEBRUIJN_KMER_INDEX_HPP
