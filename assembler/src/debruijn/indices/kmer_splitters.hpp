#pragma once
/*
 * kmer_splitters.hpp
 *
 *  Created on: May 24, 2013
 *      Author: anton
 */

#include "io/read_stream_vector.hpp"

namespace debruijn_graph {
// used for temporary reads storage during parallel reading
static const size_t READS_BUFFER_SIZE = 536870912; // 512 MB in bytes

typedef ::KMerSplitter<runtime_k::RtSeq> RtSeqKMerSplitter;

class DeBruijnKMerSplitter : public RtSeqKMerSplitter {
 protected:
  typedef KMerVector<runtime_k::RtSeq> RtSeqKMerVector;
  typedef std::vector<RtSeqKMerVector> KMerBuffer;

  size_t FillBufferFromSequence(const Sequence &seq,
                                KMerBuffer &tmp_entries, unsigned num_files) const;
  void DumpBuffers(size_t num_files, size_t nthreads,
                   std::vector<KMerBuffer> &buffers,
                   FILE **ostreams) const;

 public:
  DeBruijnKMerSplitter(const std::string &work_dir,
                       unsigned K) : RtSeqKMerSplitter(work_dir, K) {
  }
};



size_t
DeBruijnKMerSplitter::FillBufferFromSequence(const Sequence &seq,
                                             KMerBuffer &buffer, unsigned num_files) const {
  size_t kmers = 0;

  if (seq.size() < K_)
    return kmers;

  runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(K_);
  buffer[this->GetFileNumForSeq(kmer, num_files)].push_back(kmer);
  for (size_t j = K_; j < seq.size(); ++j) {
    kmer <<= seq[j];
    buffer[this->GetFileNumForSeq(kmer, num_files)].push_back(kmer);
  }

  kmers += ((seq.size() - K_ + 1) + 1);

  return kmers;
}

void DeBruijnKMerSplitter::DumpBuffers(size_t num_files, size_t nthreads,
                                       std::vector<KMerBuffer> &buffers,
                                       FILE **ostreams) const {
# pragma omp parallel for
  for (unsigned k = 0; k < num_files; ++k) {
    size_t sz = 0;
    for (size_t i = 0; i < nthreads; ++i)
      sz += buffers[i][k].size();

    KMerVector<runtime_k::RtSeq> SortBuffer(K_, sz);
    for (size_t i = 0; i < nthreads; ++i) {
      KMerBuffer &entry = buffers[i];
      for (size_t j = 0; j < entry[k].size(); ++j)
        SortBuffer.push_back(entry[k][j]);
    }
    libcxx::sort(SortBuffer.begin(), SortBuffer.end(), KMerVector<runtime_k::RtSeq>::less2_fast());
    auto it = std::unique(SortBuffer.begin(), SortBuffer.end(), KMerVector<runtime_k::RtSeq>::equal_to());

#   pragma omp critical
    {
      fwrite(SortBuffer.data(), SortBuffer.el_data_size(), it - SortBuffer.begin(), ostreams[k]);
    }
  }

  for (unsigned i = 0; i < nthreads; ++i) {
    for (unsigned j = 0; j < num_files; ++j) {
      buffers[i][j].clear();
    }
  }
}

template<class Read>
class DeBruijnReadKMerSplitter : public DeBruijnKMerSplitter {
  io::ReadStreamVector<io::IReader<Read>> &streams_;
  SingleReadStream *contigs_;

  template<class ReadStream>
  std::pair<size_t, size_t>
  FillBufferFromStream(ReadStream& stream,
                       KMerBuffer &tmp_entries,
                       unsigned num_files, size_t cell_size) const;

 public:
  DeBruijnReadKMerSplitter(const std::string &work_dir,
                           unsigned K,
                           io::ReadStreamVector< io::IReader<Read> >& streams,
                           SingleReadStream* contigs_stream = 0)
      : DeBruijnKMerSplitter(work_dir, K),
        streams_(streams), contigs_(contigs_stream) {
  }

  size_t recommended_thread_num() const {
      return streams_.size();
  }

  virtual path::files_t Split(size_t num_files);
};

template<class Read> template<class ReadStream>
std::pair<size_t, size_t>
DeBruijnReadKMerSplitter<Read>::FillBufferFromStream(ReadStream &stream,
                                                     KMerBuffer &buffer,
                                                     unsigned num_files, size_t cell_size) const {
  typename ReadStream::read_type r;
  size_t reads = 0, kmers = 0, rl = 0;

  while (!stream.eof() && kmers < num_files * cell_size) {
    stream >> r;
    rl = std::max(rl, r.size());
    reads += 1;

    kmers += FillBufferFromSequence(r.sequence(), buffer, num_files);
  }

  return std::make_pair(reads, rl);
}

template<class Read>
path::files_t DeBruijnReadKMerSplitter<Read>::Split(size_t num_files) {
  unsigned nthreads = streams_.size();

  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  // Determine the set of output files
  path::files_t out;
  for (unsigned i = 0; i < num_files; ++i)
    out.push_back(this->GetRawKMersFname(i));

  FILE** ostreams = new FILE*[num_files];
  for (unsigned i = 0; i < num_files; ++i) {
    ostreams[i] = fopen(out[i].c_str(), "wb");
    VERIFY_MSG(ostreams[i], "Cannot open temporary file to write");
  }

  size_t cell_size = READS_BUFFER_SIZE /
                     (nthreads * num_files * runtime_k::RtSeq::GetDataSize(K_) * sizeof(runtime_k::RtSeq::DataType));
  // Set sane minimum cell size
  if (cell_size < 16384)
    cell_size = 16384;
  INFO("Using cell size of " << cell_size);

  std::vector<KMerBuffer> tmp_entries(nthreads);
  for (unsigned i = 0; i < nthreads; ++i) {
    KMerBuffer &entry = tmp_entries[i];
    entry.resize(num_files, RtSeqKMerVector(K_, 1.25 * cell_size));
  }

  size_t counter = 0, rl = 0, n = 15;
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

    if (counter >> n) {
      INFO("Processed " << counter << " reads");
      n += 1;
    }
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

template<class Graph>
class DeBruijnGraphKMerSplitter : public DeBruijnKMerSplitter {
  typedef typename Graph::SmartEdgeIt EdgeIt;
  typedef typename Graph::EdgeId EdgeId;

  const Graph &g_;

  size_t FillBufferFromEdges(EdgeIt &edge,
                             KMerBuffer &tmp_entries,
                             unsigned num_files, size_t cell_size) const;

 public:
  DeBruijnGraphKMerSplitter(const std::string &work_dir,
                            unsigned K, const Graph &g)
      : DeBruijnKMerSplitter(work_dir, K), g_(g) {}

  size_t recommended_thread_num() const {
      return 1;
  }

  virtual path::files_t Split(size_t num_files);
};

template<class Graph>
size_t
DeBruijnGraphKMerSplitter<Graph>::FillBufferFromEdges(EdgeIt &edge,
                                                      KMerBuffer &buffer,
                                                      unsigned num_files, size_t cell_size) const {
  size_t seqs = 0;
  for (size_t kmers = 0; !edge.IsEnd() && kmers < num_files * cell_size; ++edge) {
    const Sequence &nucls = g_.EdgeNucls(*edge);

    kmers += FillBufferFromSequence(nucls, buffer, num_files);
    seqs += 1;
  }

  return seqs;
}

template<class Graph>
path::files_t DeBruijnGraphKMerSplitter<Graph>::Split(size_t num_files) {
  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  // Determine the set of output files
  path::files_t out;
  for (unsigned i = 0; i < num_files; ++i)
    out.push_back(this->GetRawKMersFname(i));

  FILE** ostreams = new FILE*[num_files];
  for (unsigned i = 0; i < num_files; ++i) {
    ostreams[i] = fopen(out[i].c_str(), "wb");
    VERIFY_MSG(ostreams[i], "Cannot open temporary file to write");
  }
  size_t cell_size = READS_BUFFER_SIZE /
                     (num_files * runtime_k::RtSeq::GetDataSize(K_) * sizeof(runtime_k::RtSeq::DataType));
  INFO("Using cell size of " << cell_size);

  std::vector<KMerBuffer> tmp_entries(1);
  KMerBuffer &entry = tmp_entries[0];
  entry.resize(num_files, RtSeqKMerVector(K_, 1.25 * cell_size));

  size_t counter = 0, n = 10;
  for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ) {
    counter += FillBufferFromEdges(it, tmp_entries[0], num_files, cell_size);

    DumpBuffers(num_files, 1, tmp_entries, ostreams);

    if (counter >> n) {
      INFO("Processed " << counter << " edges");
      n += 1;
    }
  }

  for (unsigned i = 0; i < num_files; ++i)
    fclose(ostreams[i]);

  delete[] ostreams;

  INFO("Used " << counter << " sequences.");

  return out;
}

}
