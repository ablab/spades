#pragma once
/*
 * edge_multi_index.hpp
 *
 *  Created on: May 24, 2013
 *      Author: anton
 */
#include "debruijn_kmer_index.hpp"
#include "kmer_splitters.hpp"

namespace debruijn_graph {
template <class Seq>
class DeBruijnEdgeMultiIndexBuilder;
template<class IdType, class Seq = runtime_k::RtSeq,
    class traits = kmer_index_traits<Seq> >

class DeBruijnEdgeMultiIndex : public EditableDeBruijnKMerIndex<vector<EdgeInfo<IdType> >, Seq, traits> {
  typedef EditableDeBruijnKMerIndex<vector<EdgeInfo<IdType> >, Seq, traits> base;
 public:
  typedef Seq                      KMer;
  typedef KMerIndex<KMer, traits>  KMerIndexT;

  DeBruijnEdgeMultiIndex(unsigned K, const std::string &workdir)
      : base(K, workdir) {}

  ~DeBruijnEdgeMultiIndex() {}

  bool ContainsInIndex(typename base::KMerIdx idx) const {
    const typename base::KMerIndexValueType &entry = base::operator[](idx);
    bool res = false;
    for (auto iter = entry.begin(); iter != entry.end(); ++iter)
      if (iter->offset_ != -1) {
        res = true;
        break;
      }
    return res;
  }

  bool ContainsInIndex(const KMer& kmer) const {
    typename base::KMerIdx idx = base::seq_idx(kmer);

    // Early exit if kmer has not been seen at all
    if (idx == base::InvalidKMerIdx)
      return false;

    // Otherwise, check, whether it's attached to any edge
    const typename base::KMerIndexValueType &entry = base::operator[](idx);
    bool res = false;
    for (auto iter = entry.begin(); iter != entry.end(); ++iter)
      if (iter->offset_ != -1) {
        res = true;
        break;
      }
    return res;
  }

  std::vector<EdgePosition> get(const KMer &kmer) const {
    typename base::KMerIdx idx = base::seq_idx(kmer);
    VERIFY(idx != base::InvalidKMerIdx);

    const typename base::KMerIndexValueType &entry = base::operator[](idx);
    return entry;
  }

  std::vector<EdgePosition> get(typename base::KMerIdx idx) const {
    const typename base::KMerIndexValueType &entry = base::operator[](idx);
    return entry;
  }

  bool DeleteIfEqual(const KMer &kmer, IdType id) {
    typename base::KMerIdx idx = base::seq_idx(kmer);

    // Early exit if kmer has not been seen at all
    if (idx == base::InvalidKMerIdx)
      return false;

    // Now we know that idx is in range. Check the edge id.
    typename base::KMerIndexValueType &entry = base::operator[](idx);

    if (entry.edgeId_ == id) {
      entry.offset_ = -1;
      return true;
    }

    return false;
  }

  void RenewKMers(const Sequence &nucls, IdType id, bool ignore_new_kmers = false) {
    VERIFY(nucls.size() >= base::K());
    KMer kmer(base::K(), nucls);

    PutInIndex(kmer, id, 0, ignore_new_kmers);
    for (size_t i = base::K(), n = nucls.size(); i < n; ++i) {
      kmer <<= nucls[i];
      PutInIndex(kmer, id, i - base::K() + 1, ignore_new_kmers);
    }
  }

  void DeleteKMers(const Sequence &nucls, IdType id) {
    VERIFY(nucls.size() >= base::K());
    KMer kmer(base::K(), nucls);
    DeleteIfEqual(kmer, id);
    for (size_t i = base::K(), n = nucls.size(); i < n; ++i) {
      kmer <<= nucls[i];
      DeleteIfEqual(kmer, id);
    }
  }

  friend class DeBruijnEdgeMultiIndexBuilder<Seq>;

private:
  void PutInIndex(const KMer &kmer, IdType id, int offset, bool ignore_new_kmer = false) {
    size_t idx = base::seq_idx(kmer);
    if (base::contains(idx, kmer)) {
      std::vector<EdgeInfo<IdType> > &entry = base::operator[](idx);
      EdgeInfo<IdType> new_entry;
      new_entry.edgeId_ = id;
      new_entry.offset_ = offset;
      entry.push_back(new_entry);
    } else {
      VERIFY(ignore_new_kmer);
    }
  }
};

template <class Seq>
class DeBruijnEdgeMultiIndexBuilder : public EditableDeBruijnKMerIndexBuilder<Seq> {
  template <class ReadStream, class IdType>
  size_t FillCoverageFromStream(ReadStream &stream,
                                DeBruijnEdgeMultiIndex<IdType, Seq> &index) const;

 public:
  template <class IdType, class Read>
  size_t BuildIndexFromStream(DeBruijnEdgeMultiIndex<IdType, Seq> &index,
                              io::ReadStreamVector<io::IReader<Read> > &streams,
                              SingleReadStream* contigs_stream = 0) const;

  template <class IdType, class Graph>
  void BuildIndexFromGraph(DeBruijnEdgeMultiIndex<IdType, Seq> &index,
                           const Graph &g) const;

  template <class IdType, class Graph>
  void UpdateIndexFromGraph(DeBruijnEdgeMultiIndex<IdType, Seq> &index,
                            const Graph &g) const;


 protected:
  template <class KMerCounter, class Index>
  void SortUniqueKMers(KMerCounter &counter, Index &index) const;

 protected:
  DECL_LOGGER("Edge MultiIndex Building");
};

template <>
class DeBruijnEdgeMultiIndexBuilder<runtime_k::RtSeq> :
      public DeBruijnKMerIndexBuilder<runtime_k::RtSeq> {
  typedef DeBruijnKMerIndexBuilder<runtime_k::RtSeq> base;

  template <class ReadStream, class IdType>
  size_t FillCoverageFromStream(ReadStream &stream,
                                DeBruijnEdgeMultiIndex<IdType, runtime_k::RtSeq> &index) const {
    size_t rl = 0;
//TODO: not needed now, implement later
//
//    while (!stream.eof()) {
//      typename ReadStream::read_type r;
//      stream >> r;
//      rl = std::max(rl, r.size());
//
//      const Sequence &seq = r.sequence();
//      if (seq.size() < K)
//        continue;
//
//      runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(K);
//
//      size_t idx = index.seq_idx(kmer);
//      VERIFY(index.contains(idx, kmer));
//#   pragma omp atomic
//      index.data_[idx].count_ += 1;
//      for (size_t j = K; j < seq.size(); ++j) {
//        kmer <<= seq[j];
//        idx = index.seq_idx(kmer);
//        VERIFY(index.contains(idx, kmer));
//
//#     pragma omp atomic
//        index.data_[idx].count_ += 1;
//      }
//    }

    return rl;
  }

 public:
  template <class IdType, class Read>
  size_t BuildIndexFromStream(DeBruijnEdgeMultiIndex<IdType, runtime_k::RtSeq> &index,
                              io::ReadStreamVector<io::IReader<Read> > &streams,
                              SingleReadStream* contigs_stream = 0) const {
    unsigned nthreads = streams.size();

    base::BuildIndexFromStream(index, streams, contigs_stream);

    // Now use the index to fill the coverage and EdgeId's
    INFO("Collecting k-mer coverage information from reads, this takes a while.");

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

  template <class IdType, class Graph>
  void BuildIndexFromGraph(DeBruijnEdgeMultiIndex<IdType, runtime_k::RtSeq> &index,
                           const Graph &g) const {
    base::BuildIndexFromGraph(index, g);

    // Now use the index to fill the coverage and EdgeId's
    INFO("Collecting k-mer coverage information from graph, this takes a while.");

    UpdateIndexFromGraph(index, g);
  }

  template <class IdType, class Graph>
  void UpdateIndexFromGraph(DeBruijnEdgeMultiIndex<IdType, runtime_k::RtSeq> &index,
                            const Graph &g) const {
    for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
      typename Graph::EdgeId edge = *it;
      index.RenewKMers(g.EdgeNucls(edge), edge);
    }
  }

 protected:
  DECL_LOGGER("Edge Index Building");
};

}
