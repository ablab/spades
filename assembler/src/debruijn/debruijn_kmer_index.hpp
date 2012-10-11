//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef DEBRUIJN_KMER_INDEX_HPP
#define DEBRUIJN_KMER_INDEX_HPP

#include "openmp_wrapper.h"
#include "standard.hpp"

#include "io/multifile_reader.hpp"

#include "mph_index/kmer_index.hpp"
#include "kmer_vector.hpp"

#include <vector>
#include <cstdlib>
#include <cstdio>

namespace debruijn_graph {

// Aux struct to count kmers during graph construction.
template<class IdType>
struct EdgeInfo {
  IdType edgeId_;
  int offset_;
  int count_;

  EdgeInfo() :
      edgeId_(), offset_(-1), count_(0) { }

  EdgeInfo(IdType edgeId, size_t offset) :
      edgeId_(edgeId), offset_(offset), count_(1) { }
};

// used for temporary reads storage during parallel reading
static const size_t READS_BUFFER_SIZE = 536870912; // 512 MB in bytes

typedef ::KMerSplitter<runtime_k::RtSeq> RtSeqKMerSplitter;

class DeBruijnKMerSplitter : public RtSeqKMerSplitter {
 protected:
  unsigned K_;

  typedef KMerVector<runtime_k::RtSeq> RtSeqKMerVector;
  typedef std::vector<RtSeqKMerVector> KMerBuffer;

  size_t FillBufferFromSequence(const Sequence &seq,
                                KMerBuffer &tmp_entries, unsigned num_files) const;
  void DumpBuffers(size_t num_files, size_t nthreads,
                   std::vector<KMerBuffer> &buffers,
                   FILE **ostreams) const;

 public:
  DeBruijnKMerSplitter(const std::string &work_dir,
                       unsigned K) : RtSeqKMerSplitter(work_dir), K_(K) {
  }
};

size_t
DeBruijnKMerSplitter::FillBufferFromSequence(const Sequence &seq,
                                             KMerBuffer &buffer, unsigned num_files) const {
  size_t kmers = 0;

  if (seq.size() < K_)
    return kmers;

  runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq::max_size>(K_);
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
    std::sort(SortBuffer.begin(), SortBuffer.end(), KMerVector<runtime_k::RtSeq>::less2_fast());
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
  INFO("Using cell size of " << cell_size);

  std::vector<KMerBuffer> tmp_entries(nthreads);
  for (unsigned i = 0; i < nthreads; ++i) {
    KMerBuffer &entry = tmp_entries[i];
    entry.resize(num_files, RtSeqKMerVector(K_, 1.25 * cell_size));
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

  size_t counter = 0;
  for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ) {
    counter += FillBufferFromEdges(it, tmp_entries[0], num_files, cell_size);

    DumpBuffers(num_files, 1, tmp_entries, ostreams);
  }

  for (unsigned i = 0; i < num_files; ++i)
    fclose(ostreams[i]);

  delete[] ostreams;

  INFO("Used " << counter << " sequences.");

  return out;
}


template<class IdType>
class DeBruijnKMerIndex {
 public:
  typedef runtime_k::RtSeq KMer;

 private:
  typedef EdgeInfo<IdType> KMerIndexValueType;
  typedef std::vector<KMerIndexValueType> KMerIndexStorageType;

  unsigned K_;
  std::string workdir_;
  KMerIndex<KMer> index_;
  KMerIndexStorageType data_;
  MMappedRecordArrayReader<typename KMer::DataType> *kmers;

 public:
  typedef typename KMerIndexStorageType::iterator value_iterator;
  typedef typename KMerIndexStorageType::const_iterator const_value_iterator;
  typedef typename MMappedRecordArrayReader<typename KMer::DataType>::iterator kmer_iterator;
  typedef typename MMappedRecordArrayReader<typename KMer::DataType>::const_iterator const_kmer_iterator;
  typedef size_t KMerIdx;
  static const size_t InvalidKMerIdx = SIZE_MAX;

  DeBruijnKMerIndex(unsigned K, const std::string &workdir)
      : K_(K), index_(K), kmers(NULL) {
    workdir_ = path::make_temp_dir(workdir, "kmeridx");
  }
  ~DeBruijnKMerIndex() {
    delete kmers;
    path::remove_dir(workdir_);
  }

  void clear() {
    index_.clear();
    data_.clear();
    KMerIndexStorageType().swap(data_);
    delete kmers;
    kmers = NULL;
  }

  unsigned K() const { return K_; }

  const KMerIndexValueType &operator[](KMerIdx idx) const {
    return data_[idx];
  }
  KMerIndexValueType &operator[](const KMer &s) {
    return operator[](index_.seq_idx(s));
  }
  const KMerIndexValueType &operator[](const KMer &s) const {
    return operator[](index_.seq_idx(s));
  }
  KMerIndexValueType &operator[](KMerIdx idx) {
    return data_[idx];
  }
  KMerIdx seq_idx(const KMer &s) const {
    KMerIdx idx = index_.seq_idx(s);
    return (contains(idx, s) ? idx : InvalidKMerIdx);
  }

  size_t size() const { return data_.size(); }

  value_iterator value_begin() {
    return data_.begin();
  }
  const_value_iterator value_begin() const {
    return data_.begin();
  }
  const_value_iterator value_cbegin() const {
    return data_.cbegin();
  }
  value_iterator value_end() {
    return data_.end();
  }
  const_value_iterator value_end() const {
    return data_.end();
  }
  const_value_iterator value_cend() const {
    return data_.cend();
  }

  kmer_iterator kmer_begin() {
    return kmers->begin();
  }
  const_kmer_iterator kmer_begin() const {
    return kmers->cbegin();
  }
  kmer_iterator kmer_end() {
    return kmers->end();
  }
  const_kmer_iterator kmer_end() const {
    return kmers->cend();
  }

  KMerIdx kmer_idx_begin() const {
    return 0;
  }

  KMerIdx kmer_idx_end() const {
    return data_.size();
  }

  bool contains(const KMer &k) const {
    KMerIdx idx = seq_idx(k);

    return idx != InvalidKMerIdx;
  }

  bool contains(KMerIdx idx) const {
    return idx < data_.size();
  }

  KMer kmer(KMerIdx idx) const {
    VERIFY(contains(idx));
    auto it = kmers->begin() + idx;
    return KMer(K_, (*it).data());
  }

  /**
   * Number of edges coming into param edge's end
   */
  unsigned RivalEdgeCount(const KMer &kmer) const {
    KMer kmer2 = kmer << 'A';
    unsigned res = 0;
    for (char c = 0; c < 4; ++c)
      if (contains(kmer2 >> c))
        res += 1;

    return res;
  }

  unsigned RivalEdgeCount(KMerIdx idx) const {
    KMer kmer2 = kmer(idx) << 'A';
    unsigned res = 0;
    for (char c = 0; c < 4; ++c)
      if (contains(kmer2 >> c))
        res += 1;

    return res;
  }

  /**
   * Number of edges going out of the param edge's end
   */
  unsigned NextEdgeCount(const KMer &kmer) const {
    unsigned res = 0;
    for (char c = 0; c < 4; ++c) {
      if (contains(kmer << c))
        res += 1;
    }

    return res;
  }

  unsigned NextEdgeCount(KMerIdx idx) const {
    KMer kmer = this->kmer(idx);

    unsigned res = 0;
    for (char c = 0; c < 4; ++c) {
      if (contains(kmer << c))
        res += 1;
    }

    return res;
  }

  KMer NextEdge(const KMer &kmer) const { // returns any next edge
    for (char c = 0; c < 4; ++c) {
      KMer s = kmer << c;
      if (contains(s))
        return s;
    }

    VERIFY_MSG(false, "Couldn't find requested edge!");
    return KMer(K_);
    // no next edges (we should request one here).
  }

  KMer NextEdge(KMerIdx idx) const { // returns any next edge
    KMer kmer = this->kmer(idx);

    for (char c = 0; c < 4; ++c) {
      KMer s = kmer << c;
      if (contains(s))
        return s;
    }

    VERIFY_MSG(false, "Couldn't find requested edge!");
    return KMer(K_);
    // no next edges (we should request one here).
  }

  std::pair<IdType, size_t> get(const KMer &kmer) const {
    KMerIdx idx = seq_idx(kmer);
    VERIFY(idx != InvalidKMerIdx);

    const KMerIndexValueType &entry = operator[](idx);
    return std::make_pair(entry.edgeId_, (size_t)entry.offset_);
  }

  std::pair<IdType, size_t> get(KMerIdx idx) const {
    const KMerIndexValueType &entry = operator[](idx);
    return std::make_pair(entry.edgeId_, (size_t)entry.offset_);
  }

  bool ContainsInIndex(KMerIdx idx) const {
    const KMerIndexValueType &entry = operator[](idx);
    return (entry.offset_ != -1);
  }

  bool ContainsInIndex(const KMer& kmer) const {
    KMerIdx idx = seq_idx(kmer);

    // Early exit if kmer has not been seen at all
    if (idx == InvalidKMerIdx)
      return false;

    // Otherwise, check, whether it's attached to any edge
    const KMerIndexValueType &entry = operator[](idx);
    return (entry.offset_ != -1);
  }

  bool DeleteIfEqual(const KMer &kmer, IdType id) {
    KMerIdx idx = seq_idx(kmer);

    // Early exit if kmer has not been seen at all
    if (idx == InvalidKMerIdx)
      return false;

    // Now we know that idx is in range. Check the edge id.
    KMerIndexValueType &entry = operator[](idx);

    if (entry.edgeId_ == id) {
      entry.offset_ = -1;
      return true;
    }

    return false;
  }

  void RenewKMers(const Sequence &nucls, IdType id, bool ignore_new_kmers = false) {
    VERIFY(nucls.size() >= K_);
    KMer kmer(K_, nucls);

    PutInIndex(kmer, id, 0, ignore_new_kmers);
    for (size_t i = K_, n = nucls.size(); i < n; ++i) {
      kmer <<= nucls[i];
      PutInIndex(kmer, id, i - K_ + 1, ignore_new_kmers);
    }
  }

  void DeleteKMers(const Sequence &nucls, IdType id) {
    VERIFY(nucls.size() >= K_);
    KMer kmer(K_, nucls);
    DeleteIfEqual(kmer, id);
    for (size_t i = K_, n = nucls.size(); i < n; ++i) {
      kmer <<= nucls[i];
      DeleteIfEqual(kmer, id);
    }
  }

  friend class DeBruijnKMerIndexBuilder;
 private:
  const std::string &workdir() const {
    return workdir_;
  }

  bool contains(KMerIdx idx, const KMer &k) const {
    // Sanity check
    if (idx >= data_.size())
      return false;

    auto it = kmers->begin() + idx;

    return (0 == memcmp(k.data(), (*it).data(), (*it).data_size()));
  }

  size_t raw_seq_idx(const typename KMerIndex<KMer>::KMerRawData &s) const {
    return index_.raw_seq_idx(s);
  }

  void PutInIndex(const KMer &kmer, IdType id, int offset, bool ignore_new_kmer = false) {
    size_t idx = seq_idx(kmer);
    if (contains(idx, kmer)) {
        KMerIndexValueType &entry = operator[](idx);
        entry.edgeId_ = id;
        entry.offset_ = offset;
    } else {
      VERIFY(ignore_new_kmer);
    }
  }
};

class DeBruijnKMerIndexBuilder {
  template<class ReadStream, class IdType>
  size_t FillCoverageFromStream(ReadStream &stream,
                                DeBruijnKMerIndex<IdType> &index) const;
  template<class IdType>
  void SortUniqueKMers(const KMerIndexBuilder<typename DeBruijnKMerIndex<IdType>::KMer> &builder,
                       DeBruijnKMerIndex<IdType> &index) const;

 public:
  template<class Read, class IdType>
  size_t BuildIndexFromStream(DeBruijnKMerIndex<IdType> &index,
                              io::ReadStreamVector<io::IReader<Read> > &streams,
                              SingleReadStream* contigs_stream = 0) const;
  template<class Graph>
  void BuildIndexFromGraph(DeBruijnKMerIndex<typename Graph::EdgeId> &index,
                           const Graph &g) const;


 private:
  DECL_LOGGER("K-mer Index Building");
};

template<class ReadStream, class IdType>
size_t
DeBruijnKMerIndexBuilder::FillCoverageFromStream(ReadStream &stream,
                                                 DeBruijnKMerIndex<IdType> &index) const {
  unsigned K = index.K();
  size_t rl = 0;

  while (!stream.eof()) {
    typename ReadStream::read_type r;
    stream >> r;
    rl = std::max(rl, r.size());

    const Sequence &seq = r.sequence();
    if (seq.size() < K)
      continue;

    runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq::max_size>(K);

    size_t idx = index.seq_idx(kmer);
    VERIFY(index.contains(idx, kmer));
#   pragma omp atomic
    index.data_[idx].count_ += 1;
    for (size_t j = K; j < seq.size(); ++j) {
      kmer <<= seq[j];
      idx = index.seq_idx(kmer);
      VERIFY(index.contains(idx, kmer));

#     pragma omp atomic
      index.data_[idx].count_ += 1;
    }
  }

  return rl;
}

template<class IdType>
void
DeBruijnKMerIndexBuilder::SortUniqueKMers(const KMerIndexBuilder<typename DeBruijnKMerIndex<IdType>::KMer> &builder,
                                          DeBruijnKMerIndex<IdType> &index) const {
  typedef typename DeBruijnKMerIndex<IdType>::KMer KMer;
  unsigned K = index.K();

  if (!index.kmers)
    index.kmers =
        new MMappedRecordArrayReader<typename KMer::DataType>(builder.GetFinalKMersFname(), KMer::GetDataSize(K), /* unlink */ true);

  size_t swaps = 0;
  INFO("Arranging kmers in hash map order");
  for (auto I = index.kmers->begin(), E = index.kmers->end(); I != E; ++I) {
    size_t cidx = I - index.kmers->begin();
    size_t kidx = index.raw_seq_idx(*I);
    while (cidx != kidx) {
      auto J = index.kmers->begin() + kidx;
      using std::swap;
      swap(*I, *J);
      swaps += 1;

      kidx = index.raw_seq_idx(*I);
    }
  }
  INFO("Done. Total swaps: " << swaps);
}

template<class Read, class IdType>
size_t
DeBruijnKMerIndexBuilder::BuildIndexFromStream(DeBruijnKMerIndex<IdType> &index,
                                               io::ReadStreamVector<io::IReader<Read> > &streams,
                                               SingleReadStream* contigs_stream) const {
  unsigned nthreads = streams.size();
  DeBruijnReadKMerSplitter<Read> splitter(index.workdir(), index.K(), streams, contigs_stream);
  KMerIndexBuilder<typename DeBruijnKMerIndex<IdType>::KMer> builder(index.workdir(), 16, streams.size());
  size_t sz = builder.BuildIndex(index.index_, splitter, /* save final */ true);

  SortUniqueKMers(builder, index);

  // Now use the index to fill the coverage and EdgeId's
  INFO("Collecting k-mer coverage information, this takes a while.");
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

template<class Graph>
void
DeBruijnKMerIndexBuilder::BuildIndexFromGraph(DeBruijnKMerIndex<typename Graph::EdgeId> &index,
                                              const Graph &g) const {
  DeBruijnGraphKMerSplitter<Graph> splitter(index.workdir(), index.K(), g);
  KMerIndexBuilder<typename DeBruijnKMerIndex<typename Graph::EdgeId>::KMer> builder(index.workdir(), 16, 1);
  size_t sz = builder.BuildIndex(index.index_, splitter, /* save final */ true);

  SortUniqueKMers(builder, index);

  // Now use the index to fill the coverage and EdgeId's
  INFO("Collecting k-mer coverage information, this takes a while.");
  index.data_.resize(sz);

  for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
    typename Graph::EdgeId edge = *it;
    index.RenewKMers(g.EdgeNucls(edge), edge);
  }
}

}

#endif // DEBRUIJN_KMER_INDEX_HPP
