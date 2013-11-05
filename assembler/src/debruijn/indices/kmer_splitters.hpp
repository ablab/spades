#pragma once
/*
 * kmer_splitters.hpp
 *
 *  Created on: May 24, 2013
 *      Author: anton
 */

#include "io/io_helper.hpp"

#include "file_limit.hpp"

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
                       unsigned K, uint32_t seed = 0)
      : RtSeqKMerSplitter(work_dir, K, seed) {
  }
};

inline size_t
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

inline
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

template<class ReadType>
class DeBruijnReadKMerSplitter : public DeBruijnKMerSplitter {
  io::ReadStreamList<ReadType> &streams_;
  io::SingleStream *contigs_;

  template<class ReadStream>
  std::pair<size_t, size_t>
  FillBufferFromStream(ReadStream& stream,
                       KMerBuffer &tmp_entries,
                       unsigned num_files, size_t cell_size) const;

  size_t rl_;

 public:
  DeBruijnReadKMerSplitter(const std::string &work_dir,
                           unsigned K, uint32_t seed,
                           io::ReadStreamList<ReadType>& streams,
                           io::SingleStream* contigs_stream = 0)
      : DeBruijnKMerSplitter(work_dir, K, seed),
        streams_(streams), contigs_(contigs_stream), rl_(0) {
  }

  virtual path::files_t Split(size_t num_files);

  size_t read_length() const { return rl_; }
};

template<class ReadType> template<class ReadStream>
std::pair<size_t, size_t>
DeBruijnReadKMerSplitter<ReadType>::FillBufferFromStream(ReadStream &stream,
                                                     KMerBuffer &buffer,
                                                     unsigned num_files, size_t cell_size) const {
  typename ReadStream::ReadT r;
  size_t reads = 0, kmers = 0, rl = 0;

  while (!stream.eof() && kmers < num_files * cell_size) {
    stream >> r;
    rl = std::max(rl, r.size());
    reads += 1;

    kmers += FillBufferFromSequence(r.sequence(), buffer, num_files);
  }

  return std::make_pair(reads, rl);
}

template<class ReadType>
path::files_t DeBruijnReadKMerSplitter<ReadType>::Split(size_t num_files) {
  unsigned nthreads = (unsigned) streams_.size();

  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  // Determine the set of output files
  path::files_t out;
  for (unsigned i = 0; i < num_files; ++i)
    out.push_back(this->GetRawKMersFname(i));

  size_t file_limit = num_files + 2*nthreads;
  size_t res = limit_file(file_limit);
  if (res < file_limit) {
    WARN("Failed to setup necessary limit for number of open files. The process might crash later on.");
    WARN("Do 'ulimit -n " << file_limit << "' in the console to overcome the limit");
  }

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
    entry.resize(num_files, RtSeqKMerVector(K_, (size_t) (1.25 * (double) cell_size)));
  }

  size_t counter = 0, rl = 0, n = 15;
  streams_.reset();
  while (!streams_.eof()) {
#   pragma omp parallel for num_threads(nthreads) reduction(+ : counter) shared(rl)
    for (size_t i = 0; i < nthreads; ++i) {
      std::pair<size_t, size_t> stats = FillBufferFromStream(streams_[i], tmp_entries[i], (unsigned) num_files, cell_size);
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
      FillBufferFromStream(*contigs_, tmp_entries[cnt], (unsigned) num_files, cell_size);
      DumpBuffers(num_files, nthreads, tmp_entries, ostreams);
      if (++cnt >= nthreads)
        cnt = 0;
    }
  }

  for (unsigned i = 0; i < num_files; ++i)
    fclose(ostreams[i]);

  delete[] ostreams;

  INFO("Used " << counter << " reads. Maximum read length " << rl);
  rl_ = rl;

  return out;
}

template<class Graph>
class DeBruijnGraphKMerSplitter : public DeBruijnKMerSplitter {
  typedef typename Graph::ConstEdgeIt EdgeIt;
  typedef typename Graph::EdgeId EdgeId;

  const Graph &g_;

  size_t FillBufferFromEdges(EdgeIt &edge,
                             KMerBuffer &tmp_entries,
                             unsigned num_files, size_t cell_size) const;

 public:
  DeBruijnGraphKMerSplitter(const std::string &work_dir,
                            unsigned K, const Graph &g)
      : DeBruijnKMerSplitter(work_dir, K), g_(g) {}

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

  size_t file_limit = num_files + 2*16;
  size_t res = limit_file(file_limit);
  if (res < file_limit) {
    WARN("Failed to setup necessary limit for number of open files. The process might crash later on.");
    WARN("Do 'ulimit -n " << file_limit << "' in the console to overcome the limit");
  }

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
  entry.resize(num_files, RtSeqKMerVector(K_, (size_t) (1.25 * (double) cell_size)));

  size_t counter = 0, n = 10;
  for (auto it = g_.ConstEdgeBegin(); !it.IsEnd(); ) {
    counter += FillBufferFromEdges(it, tmp_entries[0], (unsigned) num_files, cell_size);

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

class DeBruijnKMerKMerSplitter : public DeBruijnKMerSplitter {
  typedef MMappedFileRecordArrayIterator<runtime_k::RtSeq::DataType> kmer_iterator;

  unsigned K_source_;
  std::vector<std::string> kmers_;

  size_t FillBufferFromKMers(kmer_iterator &kmer,
                             KMerBuffer &tmp_entries,
                             unsigned num_files, size_t cell_size) const;

 public:
  DeBruijnKMerKMerSplitter(const std::string &work_dir,
                           unsigned K_target, unsigned K_source)
      : DeBruijnKMerSplitter(work_dir, K_target), K_source_(K_source) {}

  void AddKMers(const std::string &file) {
    kmers_.push_back(file);
  }

  virtual path::files_t Split(size_t num_files);
};

inline
size_t DeBruijnKMerKMerSplitter::FillBufferFromKMers(kmer_iterator &kmer,
                                                     KMerBuffer &buffer,
                                                     unsigned num_files, size_t cell_size) const {
  size_t seqs = 0;
  for (size_t kmers = 0; kmer.good() && kmers < num_files * cell_size; ++kmer) {
    Sequence nucls(runtime_k::RtSeq(K_source_, *kmer));
    kmers += FillBufferFromSequence(nucls, buffer, num_files);
    seqs += 1;
  }

  return seqs;
}

inline
path::files_t DeBruijnKMerKMerSplitter::Split(size_t num_files) {
  unsigned nthreads = (unsigned) kmers_.size();
  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  // Determine the set of output files
  path::files_t out;
  for (unsigned i = 0; i < num_files; ++i)
    out.push_back(this->GetRawKMersFname(i));

  size_t file_limit = num_files + 2*nthreads;
  size_t res = limit_file(file_limit);
  if (res < file_limit) {
    WARN("Failed to setup necessary limit for number of open files. The process might crash later on.");
    WARN("Do 'ulimit -n " << file_limit << "' in the console to overcome the limit");
  }

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
    entry.resize(num_files, RtSeqKMerVector(K_, (size_t) (1.25 * (double) cell_size)));
  }

  size_t counter = 0, n = 10;
  std::vector<kmer_iterator> its;
  its.reserve(nthreads);
  for (auto it = kmers_.begin(), et = kmers_.end(); it != et; ++it)
    its.emplace_back(*it, runtime_k::RtSeq::GetDataSize(K_source_));

  bool anygood = false;
  do {
#   pragma omp parallel for num_threads(nthreads) reduction(+ : counter)
    for (size_t i = 0; i < nthreads; ++i)
      counter += FillBufferFromKMers(its[i], tmp_entries[i], (unsigned) num_files, cell_size);

    DumpBuffers(num_files, nthreads, tmp_entries, ostreams);

    if (counter >> n) {
      INFO("Processed " << counter << " kmers");
      n += 1;
    }

    anygood = false;
    for (auto it = its.begin(), et = its.end(); it != et; ++it)
      anygood |= it->good();
  } while (anygood);

  for (unsigned i = 0; i < num_files; ++i)
    fclose(ostreams[i]);

  delete[] ostreams;

  INFO("Used " << counter << " kmers.");

  return out;
}


}
