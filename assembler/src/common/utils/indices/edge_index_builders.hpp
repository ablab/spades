//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "edge_info_updater.hpp"
#include "perfect_hash_map_builder.hpp"

namespace debruijn_graph {

template<class Index>
class GraphPositionFillingIndexBuilder {
public:
    typedef Index IndexT;
    typedef typename Index::KMer Kmer;

    template<class Graph>
    void BuildIndexFromGraph(Index &index,
                             const Graph/*T*/ &g, size_t read_buffer_size = 0) const {
        debruijn_graph::BuildIndexFromGraph(index, g, read_buffer_size);

        // Now use the index to fill the coverage and EdgeId's
        INFO("Collecting k-mer coverage information from graph, this takes a while.");
        EdgeInfoUpdater<Index, Graph> updater(g, index);
        updater.UpdateAll();
    }

};

template<typename> struct Void { typedef void type; };

template<typename T, typename Sfinae = void>
struct has_contains: std::false_type {};

template<typename T>
struct has_contains<
    T
    , typename Void<
        //decltype( std::declval<T&>().contains(typename T::KMerIdx(0), typename T::KMer()) )
        decltype( ((T*)(0))->contains(*((typename T::KeyWithHash*)(0))) )
    >::type
>: std::true_type {};

template <class Builder>
class CoverageFillingEdgeIndexBuilder : public Builder {
    typedef Builder base;
 public:
    typedef typename Builder::IndexT IndexT;
    typedef typename IndexT::KMer Kmer;
    typedef typename IndexT::KMerIdx KmerIdx;
    typedef typename IndexT::KeyWithHash KeyWithHash;

 private:


    bool ContainsWrap(bool check_contains, IndexT& index, const KeyWithHash &kwh, std::true_type) const {
        return !check_contains || index.contains(kwh);
    }

    bool ContainsWrap(bool /*check_contains*/, IndexT&/* index*/, const KeyWithHash &/*kwh*/, std::false_type) const {
        VERIFY(false);
//        VERIFY(!check_contains);
        return true;
    }

    template<class ReadStream>
    size_t FillCoverageFromStream(ReadStream &stream,
                                  IndexT &index, bool check_contains) const {
        unsigned k = index.k();
        size_t rl = 0;

        while (!stream.eof()) {
            typename ReadStream::ReadT r;
            stream >> r;
            rl = std::max(rl, r.size());

            const Sequence &seq = r.sequence();
            if (seq.size() < k)
                continue;

            KeyWithHash kwh = index.ConstructKWH(seq.start<Kmer>(k) >> 'A');
            for (size_t j = k - 1; j < seq.size(); ++j) {
                kwh <<= seq[j];
                //contains is not used since index might be still empty here
                if (kwh.is_minimal() && index.valid(kwh) && ContainsWrap(check_contains, index, kwh, has_contains<IndexT>())) {
#     pragma omp atomic
                    index.get_raw_value_reference(kwh).count += 1;
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
                                io::SingleStream* contigs_stream = 0) const {
        debruijn_graph::BuildIndexFromStream(index, streams, contigs_stream);

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
    typedef typename Index::KMer Kmer;
    typedef typename Index::KMerIdx KMerIdx;
    typedef typename Index::traits_t traits_t;
    typedef CoverageFillingEdgeIndexBuilder<Index> CoverageFillingEdgeIndexBuilderT;
    typedef GraphPositionFillingIndexBuilder<Index> GraphPositionFillingIndexBuilderT;
    typedef CoverageFillingEdgeIndexBuilder<GraphPositionFillingIndexBuilderT> CoverageAndGraphPositionFillingIndexBuilderT;
};

}
