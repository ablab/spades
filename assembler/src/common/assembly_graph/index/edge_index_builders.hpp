//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "edge_info_updater.hpp"
#include "utils/ph_map/perfect_hash_map_builder.hpp"

namespace debruijn_graph {

template<class Graph, class KmerFilter>
class DeBruijnGraphKMerSplitter : public utils::DeBruijnKMerSplitter<KmerFilter> {
    typedef typename omnigraph::GraphEdgeIterator<Graph> EdgeIt;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename adt::iterator_range<EdgeIt> EdgeRange;
    using typename utils::DeBruijnKMerSplitter<KmerFilter>::RawKMers;

    const Graph &g_;

    size_t FillBufferFromEdges(EdgeRange &r, unsigned thread_id);

public:
    DeBruijnGraphKMerSplitter(fs::TmpDir work_dir,
                              unsigned K, const Graph &g,
                              size_t read_buffer_size = 0)
            : utils::DeBruijnKMerSplitter<KmerFilter>(work_dir, K, KmerFilter(), read_buffer_size),
              g_(g) {}

    RawKMers Split(size_t num_files, unsigned nthreads) override;
};

template<class Graph, class KmerFilter>
size_t
DeBruijnGraphKMerSplitter<Graph, KmerFilter>::FillBufferFromEdges(EdgeRange &r,
                                                                  unsigned thread_id) {
    size_t seqs = 0;
    for (auto &it = r.begin(); it != r.end(); ++it) {
        const Sequence &nucls = g_.EdgeNucls(*it);

        seqs += 1;
        if (this->FillBufferFromSequence(nucls, thread_id))
            break;
    }

    return seqs;
}

template<class Graph, class KmerFilter>
typename DeBruijnGraphKMerSplitter<Graph, KmerFilter>::RawKMers
DeBruijnGraphKMerSplitter<Graph, KmerFilter>::Split(size_t num_files, unsigned nthreads) {
    auto out = this->PrepareBuffers(num_files, nthreads, this->read_buffer_size_);

    omnigraph::IterationHelper<Graph, EdgeId> edges(g_);
    auto its = edges.Chunks(nthreads);

    // Turn chunks into iterator ranges
    std::vector<EdgeRange> ranges;
    for (size_t i = 0; i < its.size() - 1; ++i)
        ranges.emplace_back(its[i], its[i+1]);

    VERIFY(ranges.size() <= nthreads);

    size_t counter = 0, n = 10;
    while (!std::all_of(ranges.begin(), ranges.end(),
                        [](const EdgeRange &r) { return r.begin() == r.end(); })) {
#       pragma omp parallel for num_threads(nthreads) reduction(+ : counter)
        for (size_t i = 0; i < ranges.size(); ++i)
            counter += FillBufferFromEdges(ranges[i], omp_get_thread_num());

        this->DumpBuffers(out);

        if (counter >> n) {
            INFO("Processed " << counter << " edges");
            n += 1;
        }
    }

    INFO("Used " << counter << " sequences.");

    this->ClearBuffers();

    return out;
}

template<class Index>
class GraphPositionFillingIndexBuilder {
public:
    typedef Index IndexT;
    typedef typename Index::KMer Kmer;

    template<class Graph>
    void BuildIndexFromGraph(Index &index, const Graph &g,
                             fs::TmpDir workdir, size_t read_buffer_size = 0) const {
        unsigned nthreads = omp_get_max_threads();

        DeBruijnGraphKMerSplitter<Graph,
                                  utils::StoringTypeFilter<typename Index::storing_type>>
                splitter(workdir, index.k(), g, read_buffer_size);
        utils::KMerDiskCounter<RtSeq> counter(workdir, splitter);
        BuildIndex(index, counter, 16, nthreads);

        // Now use the index to fill the coverage and EdgeId's
        INFO("Collecting edge information from graph, this takes a while.");
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
    void FillCoverageFromStream(ReadStream &stream,
                                IndexT &index, bool check_contains) const {
        unsigned k = index.k();

        while (!stream.eof()) {
            typename ReadStream::ReadT r;
            stream >> r;

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
    }

 public:

    template<class Streams>
    void ParallelFillCoverage(IndexT &index,
                              Streams &streams,
                              bool check_contains = true) const {
        INFO("Collecting k-mer coverage information from reads, this takes a while.");
        unsigned nthreads = (unsigned) streams.size();
        streams.reset();
#pragma omp parallel for num_threads(nthreads)
        for (size_t i = 0; i < nthreads; ++i) {
            FillCoverageFromStream(streams[i], index, check_contains);
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
    }
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
