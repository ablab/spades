#pragma once
/*
 * edge_index.hpp
 *
 *  Created on: May 24, 2013
 *      Author: anton
 */

#include "debruijn_kmer_index.hpp"
#include "kmer_splitters.hpp"

namespace debruijn_graph {

template <class Seq>
class DeBruijnEdgeIndexBuilder;

// Aux struct to count kmers during graph construction.
template<class IdType>
struct EdgeInfo {
    IdType edgeId_;
    int offset_;
    int count_;

    EdgeInfo(IdType edgeId = IdType(), int offset = -1, int count = 0) :
            edgeId_(edgeId), offset_(offset), count_(count) { }
};

template<class Graph, class Seq = runtime_k::RtSeq, class traits = kmer_index_traits<Seq>>
class DeBruijnEdgeIndex : public DeBruijnKMerIndex<EdgeInfo<typename Graph::EdgeId>, traits> {
    typedef typename Graph::EdgeId IdType;
    typedef DeBruijnKMerIndex<EdgeInfo<IdType>, traits> base;
    const Graph &graph_;
  public:
    typedef typename traits::SeqType KMer;
    typedef KMerIndex<traits>        KMerIndexT;

    DeBruijnEdgeIndex(unsigned K, const Graph &graph, const std::string &workdir)
            : base(K, workdir), graph_(graph) {}

    ~DeBruijnEdgeIndex() {}

    bool contains(const KMer& kmer) const {
        typename base::KMerIdx idx = base::seq_idx(kmer);
        return contains(idx, kmer);
    }

    KMer kmer(typename base::KMerIdx idx) const {
        VERIFY(idx < this->size());
        const typename base::KMerIndexValueType &entry = base::operator[](idx);
        VERIFY(entry.offset_ != -1);
        return KMer(this->K_, graph_.EdgeNucls(entry.edgeId_), entry.offset_);
    }

    bool contains(typename base::KMerIdx idx, const KMer &k) const {
        // Sanity check
        if (idx == base::InvalidKMerIdx || idx >= this->size())
            return false;

        const typename base::KMerIndexValueType &entry = base::operator[](idx);

        if (entry.offset_ == -1)
            return false;

        return k == KMer(this->K_, graph_.EdgeNucls(entry.edgeId_), entry.offset_);
    }

  public:
    /**
     * Number of edges coming into param edge's end
     */
    unsigned RivalEdgeCount(const KMer &kmer) const {
        KMer kmer2 = kmer << 'A';
        unsigned res = 0;
        for (char c = 0; c < 4; ++c)
            if (RandomAccessContains(kmer2 >> c))
                res += 1;

        return res;
    }
  private:
    bool RandomAccessContains(const KMer &kmer) const {
        VERIFY(false);
        return false;
    }

    KMer RandomAccesskmer(typename base::KMerIdx idx) const {
        VERIFY(false);
        return KMer();
    }

  public:
    unsigned RivalEdgeCount(typename base::KMerIdx idx) const {
        KMer kmer2 = RandomAccesskmer(idx) << 'A';
        unsigned res = 0;
        for (char c = 0; c < 4; ++c)
            if (RandomAccessContains(kmer2 >> c))
                res += 1;

        return res;
    }

    /**
     * Number of edges going out of the param edge's end
     */
    unsigned NextEdgeCount(const KMer &kmer) const {
        unsigned res = 0;
        for (char c = 0; c < 4; ++c)
            if (RandomAccessContains(kmer << c))
                res += 1;

        return res;
    }

    unsigned NextEdgeCount(typename base::KMerIdx idx) const {
        KMer kmer = RandomAccesskmer(idx);

        unsigned res = 0;
        for (char c = 0; c < 4; ++c)
            if (RandomAccessContains(kmer << c))
                res += 1;

        return res;
    }

    KMer NextEdge(const KMer &kmer) const { // returns any next edge
        for (char c = 0; c < 4; ++c) {
            KMer s = kmer << c;
            typename base::KMerIdx idx = base::seq_idx(s);
            if (RandomAccessContains(idx))
                return RandomAccesskmer(idx);
        }

        VERIFY_MSG(false, "Couldn't find requested edge!");
        return KMer(base::K());
        // no next edges (we should request one here).
    }

    KMer NextEdge(typename base::KMerIdx idx) const { // returns any next edge
        KMer kmer = RandomAccesskmer(idx);

        for (char c = 0; c < 4; ++c) {
            KMer s = kmer << c;
            if (RandomAccessContains(s))
                return s;
        }

        VERIFY_MSG(false, "Couldn't find requested edge!");
        return KMer(base::K());
        // no next edges (we should request one here).
    }

    std::pair<IdType, size_t> get(const KMer &kmer) const {
        typename base::KMerIdx idx = base::seq_idx(kmer);
        VERIFY(this->contains(idx, kmer));

        const EdgeInfo<IdType> &entry = base::operator[](idx);
        return std::make_pair(entry.edgeId_, (size_t)entry.offset_);
    }

    std::pair<IdType, size_t> get(typename base::KMerIdx idx) const {
        const EdgeInfo<IdType> &entry = base::operator[](idx);
        return std::make_pair(entry.edgeId_, (size_t)entry.offset_);
    }

    bool DeleteIfEqual(const KMer &kmer, IdType id) {
        typename base::KMerIdx idx = base::seq_idx(kmer);
        if (!contains(idx, kmer))
            return false;

        EdgeInfo<IdType> &entry = base::operator[](idx);
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

    template<class Writer>
    void BinWrite(Writer &writer) const {
        base::index_.serialize(writer);
        size_t sz = base::data_.size();
        writer.write((char*)&sz, sizeof(sz));
        for (size_t i = 0; i < sz; ++i)
            writer.write((char*)&(base::data_[i].count_), sizeof(base::data_[0].count_));
        traits::raw_serialize(writer, base::kmers);
    }

    template<class Reader>
    void BinRead(Reader &reader, const std::string &FileName) {
        base::clear();
        base::index_.deserialize(reader);
        size_t sz = 0;
        reader.read((char*)&sz, sizeof(sz));
        base::data_.resize(sz);
        for (size_t i = 0; i < sz; ++i)
            reader.read((char*)&(base::data_[i].count_), sizeof(base::data_[0].count_));
        base::kmers = traits::raw_deserialize(reader, FileName);
    }

    friend class DeBruijnEdgeIndexBuilder<KMer>;

  private:
    void PutInIndex(const KMer &kmer, IdType id, int offset, bool ignore_new_kmer = false) {
        size_t idx = base::seq_idx(kmer);
        if (contains(idx, kmer)) {
            EdgeInfo<IdType> &entry = base::operator[](idx);
            entry.edgeId_ = id;
            entry.offset_ = offset;
        }
        //    } else {
        //      VERIFY(ignore_new_kmer);
        //      idx = base::insert(kmer, EdgeInfo<IdType>(id, offset, 1));
        //
        //      VERIFY(idx != base::InvalidKMerIdx);
        //    }
    }
};

template <class Seq>
class DeBruijnEdgeIndexBuilder : public DeBruijnKMerIndexBuilder<kmer_index_traits<Seq>> {
    template <class ReadStream, class Graph>
    size_t FillCoverageFromStream(ReadStream &stream,
                                  DeBruijnEdgeIndex<Graph, Seq> &index) const;

  public:
    template <class Graph, class Read>
            size_t BuildIndexFromStream(DeBruijnEdgeIndex<Graph, Seq> &index,
                                io::ReadStreamVector<io::IReader<Read> > &streams,
                                SingleReadStream* contigs_stream = 0) const;

    template <class Graph>
            void BuildIndexFromGraph(DeBruijnEdgeIndex<Graph, Seq> &index,
                             const Graph &g) const;

    template <class Graph>
            void UpdateIndexFromGraph(DeBruijnEdgeIndex<Graph, Seq> &index,
                              const Graph &g) const;


  protected:
    DECL_LOGGER("Edge Index Building");
};

template <>
class DeBruijnEdgeIndexBuilder<runtime_k::RtSeq> :
            public DeBruijnKMerIndexBuilder<kmer_index_traits<runtime_k::RtSeq> > {
    typedef DeBruijnKMerIndexBuilder<kmer_index_traits<runtime_k::RtSeq> > base;

    template <class ReadStream, class Graph>
    size_t FillCoverageFromStream(ReadStream &stream,
                                  DeBruijnEdgeIndex<Graph, runtime_k::RtSeq> &index) const {
        unsigned K = index.K();
        size_t rl = 0;

        while (!stream.eof()) {
            typename ReadStream::read_type r;
            stream >> r;
            rl = std::max(rl, r.size());

            const Sequence &seq = r.sequence();
            if (seq.size() < K)
                continue;

            runtime_k::RtSeq kmer = seq.start<runtime_k::RtSeq>(K);

            size_t idx = index.seq_idx(kmer);
            if (index.contains(idx, kmer)) {
#   pragma omp atomic
                index.data_[idx].count_ += 1;
            }
            for (size_t j = K; j < seq.size(); ++j) {
                kmer <<= seq[j];
                idx = index.seq_idx(kmer);
                if (index.contains(idx, kmer)) {
#     pragma omp atomic
                    index.data_[idx].count_ += 1;
                }
            }
        }

        return rl;
    }

  public:

    template<class Graph, class Read>
    size_t ParallelFillCoverage(
        DeBruijnEdgeIndex<Graph, runtime_k::RtSeq> &index,
        io::ReadStreamVector<io::IReader<Read> > &streams,
        SingleReadStream* contigs_stream = 0) const {
        INFO("Collecting k-mer coverage information from reads, this takes a while.");

        unsigned nthreads = streams.size();
        size_t rl = 0;
        streams.reset();
#pragma omp parallel for num_threads(nthreads) shared(rl)
        for (size_t i = 0; i < nthreads; ++i) {
            size_t crl = FillCoverageFromStream(streams[i], index);

            // There is no max reduction in C/C++ OpenMP... Only in FORTRAN :(
#pragma omp flush(rl)
            if (crl > rl)
#pragma omp critical
            {
                rl = std::max(rl, crl);
            }
        }

        // Contigs have zero coverage!
#if 0
        if (contigs_stream) {
            contigs_stream->reset();
            FillCoverageFromStream(*contigs_stream, index);
        }
#endif

#ifndef NDEBUG
        for (auto idx = index.kmer_idx_begin(), eidx = index.kmer_idx_end();
             idx != eidx; ++idx) {
            runtime_k::RtSeq k = index.kmer(idx);

            VERIFY(index[k].count_ == index[!k].count_);
        }
#endif
        return rl;
    }

    template<class Graph, class Read>
    size_t BuildIndexFromStream(
        DeBruijnEdgeIndex<Graph, runtime_k::RtSeq> &index,
        io::ReadStreamVector<io::IReader<Read> > &streams,
        SingleReadStream* contigs_stream = 0) const {
        base::BuildIndexFromStream(index, streams, contigs_stream);

        // Now use the index to fill the coverage and EdgeId's
        return ParallelFillCoverage(index, streams, contigs_stream);
    }

    template<class Graph, class Read>
    size_t BuildIndexWithCoverageFromGraph(Graph &graph,
                                           DeBruijnEdgeIndex<Graph, runtime_k::RtSeq> &index,
                                           io::ReadStreamVector<io::IReader<Read> > &streams,
                                           SingleReadStream* contigs_stream = 0) const {
        BuildIndexFromGraph(index, graph);

        // Now use the index to fill the coverage and EdgeId's
        return ParallelFillCoverage(index, streams, contigs_stream);
    }

    template <class Graph>
    void BuildIndexFromGraph(DeBruijnEdgeIndex<Graph, runtime_k::RtSeq> &index,
                             const Graph &g) const {
        base::BuildIndexFromGraph(index, g);

        // Now use the index to fill the coverage and EdgeId's
        INFO("Collecting k-mer coverage information from graph, this takes a while.");

        UpdateIndexFromGraph(index, g);
    }

    template <class Graph>
    void UpdateIndexFromGraph(DeBruijnEdgeIndex<Graph, runtime_k::RtSeq> &index,
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
