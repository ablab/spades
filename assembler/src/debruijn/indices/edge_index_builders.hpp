#pragma once

#include "debruijn_edge_index.hpp"

namespace debruijn_graph {

template <class Builder>
class GraphPositionFillingIndexBuilder : public Builder {
    typedef Builder base;
public:
    typedef typename Builder::IndexT IndexT;
    typedef typename IndexT::KMer Kmer;
//    typedef typename IndexT::GraphT GraphT;

    template<class Graph>
    void BuildIndexFromGraph(IndexT &index,
                             const Graph/*T*/ &g) const {
        base::BuildIndexFromGraph(index, g);

        // Now use the index to fill the coverage and EdgeId's
        INFO("Collecting k-mer coverage information from graph, this takes a while.");
        EdgeInfoUpdater<IndexT, Graph> updater(g, index);
        updater.UpdateAll();
    }

};

template <class Builder>
class CoverageFillingEdgeIndexBuilder : public Builder {
    typedef Builder base;
 public:
    typedef typename Builder::IndexT IndexT;
    typedef typename IndexT::KMer Kmer;
    typedef typename IndexT::KMerIdx KmerIdx;

 private:

    bool ContainsWrap(bool check_contains, IndexT& index, KmerIdx idx, const Kmer& kmer, key_storing_tag) const {
        return !check_contains || index.contains(idx, kmer);
    }

    bool ContainsWrap(bool check_contains, IndexT& index, KmerIdx idx, const Kmer& kmer, key_free_tag) const {
        VERIFY(!check_contains);
        return true;
    }

    template<class ReadStream>
    size_t FillCoverageFromStream(ReadStream &stream,
                                  IndexT &index, bool check_contains) const {
        unsigned k = index.k();
        size_t rl = 0;

        while (!stream.eof()) {
            typename ReadStream::read_type r;
            stream >> r;
            rl = std::max(rl, r.size());

            const Sequence &seq = r.sequence();
            if (seq.size() < k)
                continue;

            Kmer kmer = seq.start<Kmer>(k);
            kmer >>= 'A';
            for (size_t j = k - 1; j < seq.size(); ++j) {
                kmer <<= seq[j];
                KmerIdx idx = index.seq_idx(kmer);
                typedef typename IndexT::key_storing_policy_tag key_storing_policy;
                //contains is not used since index might be still empty here
                if (index.valid_idx(idx) && ContainsWrap(check_contains, index, idx, kmer, key_storing_policy())) {
#     pragma omp atomic
                    index[idx].count += 1;
                }
            }
        }

        return rl;
    }

 public:

    template<class Streams>
    size_t ParallelFillCoverage(IndexT &index,
                                Streams &streams,
                                bool check_contains = true) const {
        INFO("Collecting k-mer coverage information from reads, this takes a while.");
        unsigned nthreads = (unsigned) streams.size();
        size_t rl = 0;
        streams.reset();
#pragma omp parallel for num_threads(nthreads) shared(rl)
        for (size_t i = 0; i < nthreads; ++i) {
            size_t crl = FillCoverageFromStream(streams[i], index, check_contains);

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
            FillCoverageFromStream(*contigs_stream, index, check_contains);
        }
#endif

//todo if this verify is neede, put it outside
//#ifndef NDEBUG
//        for (auto idx = index.kmer_idx_begin(), eidx = index.kmer_idx_end();
//             idx != eidx; ++idx) {
//
//            Kmer k = index.kmer(idx);
//
//            VERIFY(index[k].count == index[!k].count);
//        }
//#endif

        return rl;
    }

    template<class Streams>
    size_t BuildIndexFromStream(IndexT &index,
                                Streams &streams,
                                SingleReadStream* contigs_stream = 0) const {
        base::BuildIndexFromStream(index, streams, contigs_stream);

        return ParallelFillCoverage(index, streams, false);
    }

//    template<class Streams>
//    size_t BuildIndexWithCoverageFromGraph(
//            GraphT &graph, IndexT &index,
//            Streams &streams,
//            SingleReadStream* contigs_stream = 0) const {
//        this->BuildIndexFromGraph(index, graph);
//
//        return ParallelFillCoverage(index, streams, contigs_stream, true);
//    }
};

template<class Index>
struct EdgeIndexHelper {
    typedef Index IndexT;
    typedef typename IndexT::KMer Kmer;
    typedef typename IndexT::KMerIdx KMerIdx;
    typedef typename IndexT::traits_t traits_t;
//    typedef typename IndexT::IdType IdType;
    typedef DeBruijnStreamKMerIndexBuilder<Kmer, IndexT> DeBruijnStreamKMerIndexBuilderT;
    typedef CoverageFillingEdgeIndexBuilder<DeBruijnStreamKMerIndexBuilderT> CoverageFillingEdgeIndexBuilderT;
    typedef DeBruijnGraphKMerIndexBuilder<IndexT> DeBruijnGraphKMerIndexBuilderT;
    typedef GraphPositionFillingIndexBuilder<DeBruijnGraphKMerIndexBuilderT> GraphPositionFillingIndexBuilderT;
    typedef CoverageFillingEdgeIndexBuilder<GraphPositionFillingIndexBuilderT> CoverageAndGraphPositionFillingIndexBuilderT;
};

}
