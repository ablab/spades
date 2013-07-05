#pragma once

#include "debruijn_edge_index.hpp"

namespace debruijn_graph {

template <class Builder>
class GraphPositionFillingIndexBuilder : public Builder {
    typedef Builder base;
public:
    typedef typename Builder::IndexT IndexT;
    typedef typename IndexT::KMer Kmer;
    typedef typename IndexT::GraphT GraphT;

    void BuildIndexFromGraph(IndexT &index,
                             const GraphT &g) const {
        base::BuildIndexFromGraph(index, g);

        // Now use the index to fill the coverage and EdgeId's
        INFO("Collecting k-mer coverage information from graph, this takes a while.");
        EdgeInfoUpdater<IndexT> updater(g, index);
        updater.UpdateAll();
    }

};

template <class Builder>
class CoverageFillingEdgeIndexBuilder : public Builder {
    typedef Builder base;
 public:
    typedef typename Builder::IndexT IndexT;
    typedef typename IndexT::KMer Kmer;
    typedef typename IndexT::GraphT GraphT;

    template<class ReadStream>
    size_t FillCoverageFromStream(ReadStream &stream,
                                  IndexT &index) const {
        unsigned K = index.K();
        size_t rl = 0;

        while (!stream.eof()) {
            typename ReadStream::read_type r;
            stream >> r;
            rl = std::max(rl, r.size());

            const Sequence &seq = r.sequence();
            if (seq.size() < K)
                continue;

            Kmer kmer = seq.start<Kmer>(K);

            size_t idx = index.seq_idx(kmer);
            if (index.contains(idx, kmer)) {
#   pragma omp atomic
                index[idx].count += 1;
            }
            for (size_t j = K; j < seq.size(); ++j) {
                kmer <<= seq[j];
                idx = index.seq_idx(kmer);
                if (index.contains(idx, kmer)) {
#     pragma omp atomic
                    index[idx].count += 1;
                }
            }
        }

        return rl;
    }

    template<class Streams>
    size_t ParallelFillCoverage(IndexT &index,
                                Streams &streams,
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

            VERIFY(index[k].count == index[!k].count);
        }
#endif
        return rl;
    }

 public:

    template<class Streams>
    size_t BuildIndexFromStream(IndexT &index,
                                /*io::ReadStreamVector<io::IReader<Read> >*/Streams &streams,
                                SingleReadStream* contigs_stream = 0) const {
        base::BuildIndexFromStream(index, streams, contigs_stream);

        return ParallelFillCoverage(index, streams, contigs_stream);
    }

    template<class Streams>
    size_t BuildIndexWithCoverageFromGraph(
            GraphT &graph, IndexT &index,
            /*io::ReadStreamVector<io::IReader<Read> >*/Streams &streams,
            SingleReadStream* contigs_stream = 0) const {
        BuildIndexFromGraph(index, graph);

        return ParallelFillCoverage(index, streams, contigs_stream);
    }
};

template<class Index>
struct EdgeIndexHelper {
    typedef Index IndexT;
//    typedef typename IndexT::GraphT GraphT;
//    typedef typename IndexT::KMer Kmer;
//    typedef typename IndexT::KMerIdx KMerIdx;
//    typedef typename Index::InnerIndexT InnerIndexT;
    typedef typename IndexT::KMer Kmer;
    typedef typename IndexT::KMerIdx KMerIdx;
    typedef typename IndexT::traits_t traits_t;
    typedef typename IndexT::IdType IdType;
    typedef typename IndexT::BuilderT BuilderT;
    typedef DeBruijnKMerIndexBuilder<Kmer, BuilderT> DeBruijnKMerIndexBuilderT;
    typedef DeBruijnGraphKMerIndexBuilder<DeBruijnKMerIndexBuilderT> DeBruijnGraphKMerIndexBuilderT;
    typedef GraphPositionFillingIndexBuilder<DeBruijnGraphKMerIndexBuilderT> GraphPositionFillingIndexBuilderT;
    typedef CoverageFillingEdgeIndexBuilder<DeBruijnGraphKMerIndexBuilderT> CoverageFillingEdgeIndexBuilderT;
    typedef CoverageFillingEdgeIndexBuilder<GraphPositionFillingIndexBuilderT> CoverageAndGraphPositionFillingIndexBuilderT;
};

}
