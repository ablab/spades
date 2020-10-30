//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "edge_info_updater.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/core/kmer_iterator.hpp"

#include "utils/kmer_mph/kmer_splitters.hpp"
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
    auto ranges = edges.Ranges(nthreads);

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

template<class Graph, class KmerFilter>
class DeBruijnEdgeKMerSplitter : public utils::DeBruijnKMerSplitter<KmerFilter> {
    typedef typename Graph::EdgeId EdgeId;
    using typename utils::DeBruijnKMerSplitter<KmerFilter>::RawKMers;

    const Graph &g_;
    const std::vector<EdgeId> &edges_;

 public:
    DeBruijnEdgeKMerSplitter(fs::TmpDir work_dir,
                             unsigned K, const Graph &g, const std::vector<EdgeId> &edges,
                             size_t read_buffer_size = 0)
        : utils::DeBruijnKMerSplitter<KmerFilter>(work_dir, K, KmerFilter(), read_buffer_size),
          g_(g), edges_(edges) {}

    RawKMers Split(size_t num_files, unsigned nthreads) override;
};

template<class Graph, class KmerFilter>
typename DeBruijnEdgeKMerSplitter<Graph, KmerFilter>::RawKMers
DeBruijnEdgeKMerSplitter<Graph, KmerFilter>::Split(size_t num_files, unsigned nthreads) {
    auto out = this->PrepareBuffers(num_files, nthreads, this->read_buffer_size_);

    std::vector<std::pair<size_t, size_t>> ranges(nthreads);
    size_t chunk_size = std::max(size_t(1), edges_.size()/(10*nthreads));
    size_t range_begin = 0, range_end = 0;
    while (range_end < edges_.size()) {
        range_begin = range_end;
        range_end += chunk_size;
        if (range_end > edges_.size()) {
            range_end = edges_.size();
        }

        ranges.emplace_back(range_begin, range_end);
    }

    size_t counter = 0, n = 10;
    while (!std::all_of(ranges.begin(), ranges.end(),
                        [](const std::pair<size_t, size_t> &r) { return r.first == r.second; })) {
#   pragma omp parallel for num_threads(nthreads) reduction(+ : counter)
        for (size_t chunk_id = 0; chunk_id < ranges.size(); ++chunk_id) {
            while (ranges[chunk_id].first < ranges[chunk_id].second) {
                const Sequence &nucls = g_.EdgeNucls(edges_[ranges[chunk_id].first]);
                ranges[chunk_id].first += 1;
                counter += 1;

                if (this->FillBufferFromSequence(nucls, omp_get_thread_num())) {
                    break;
                }
            }
        }
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

    template<class EdgeRange>
    class KMerGraphStorage {
      public:
        KMerGraphStorage(const Graph &g, unsigned k,
                         std::vector<EdgeRange> buckets)
                : k_(k), g_(g), buckets_(std::move(buckets)) {
            segment_policy_.reset(1);

            sizes_.resize(buckets_.size());
#           pragma omp parallel for
            for (size_t i = 0; i < buckets_.size(); ++i) {
                size_t sz = std::distance(bucket_begin(i), bucket_end(i));
                sizes_[i] = sz;
            }
            for (size_t sz : sizes_)
                total_kmers_ += sz;
        }

        unsigned k() const { return k_; }
        size_t total_kmers() const {  return total_kmers_;  }

        size_t bucket_size(size_t i) const { return sizes_.at(i); }

        auto bucket_begin(size_t i) const {
            return boost::make_filter_iterator(utils::StoringTypeFilter<typename Index::storing_type>(),
                                               kmer_begin(buckets_.at(i).begin(), buckets_.at(i).end(), g_, k_),
                                               kmer_end(buckets_.at(i).end(), g_));
        }

        auto bucket_end(size_t i) const {
            return boost::make_filter_iterator(utils::StoringTypeFilter<typename Index::storing_type>(),
                                               kmer_end(buckets_.at(i).end(), g_),
                                               kmer_end(buckets_.at(i).end(), g_));
        }

        size_t num_buckets() const { return buckets_.size(); }
        kmer::KMerSegmentPolicy<Kmer> segment_policy() const { return segment_policy_; }

      private:
        unsigned k_;
        const Graph &g_;
        size_t total_kmers_ = 0;
        std::vector<EdgeRange> buckets_;
        std::vector<size_t> sizes_;
        kmer::KMerSegmentPolicy<Kmer> segment_policy_;
    };

    class KMerFullGraphStorage : public KMerGraphStorage<omnigraph::IterationHelper<Graph, EdgeId>::EdgeRange> {
        using IterationHelper = omnigraph::IterationHelper<Graph, EdgeId>;
        using base = KMerGraphStorage<IterationHelper::EdgeRange>;

      public:
        KMerFullGraphStorage(const Graph &g, unsigned k, unsigned num_buckets)
                : base(g, k, IterationHelper(g).Ranges(num_buckets)) {}
    };

    class KMerPartialGraphStorage : public KMerGraphStorage<adt::iterator_range<std::vector<EdgeId>::const_iterator>> {
        using EdgeRange = adt::iterator_range<std::vector<EdgeId>::const_iterator>;
        using base = KMerGraphStorage<EdgeRange>;

        std::vector<EdgeRange> ranges(const std::vector<EdgeId> &edges, unsigned num_buckets) const {
            if (num_buckets == 1)
                return { adt::make_range(edges.begin(), edges.end()) };

            std::vector<EdgeRange> res;
            size_t chunk_size = std::max(size_t(1), edges.size() / num_buckets);
            size_t range_begin = 0, range_end = 0;
            while (range_end < edges.size()) {
                range_begin = range_end;
                range_end += chunk_size;
                if (range_end > edges.size())
                    range_end = edges.size();

                res.emplace_back(edges.begin() + range_begin, edges.begin() + range_end);
            }

            return res;
        }

      public:
        KMerPartialGraphStorage(const Graph &g, unsigned k,
                                const std::vector<EdgeId> &edges, unsigned num_buckets)
                : base(g, k, ranges(edges, num_buckets)) {}
    };

    void BuildIndexFromGraph(Index &index, const Graph &g) const {
        unsigned nthreads = omp_get_max_threads();

        KMerFullGraphStorage storage(g, index.k(), 10 * nthreads);

        BuildIndex(index, storage, nthreads);

        // Now use the index to fill the coverage and EdgeId's
        INFO("Collecting edge information from graph, this takes a while.");
        EdgeInfoUpdater<Graph>().UpdateAll(g, index);
    }

    void BuildIndexFromGraph(Index &index, const Graph &g,
                             const std::vector<EdgeId> &edges) const {
        unsigned nthreads = omp_get_max_threads();

        KMerPartialGraphStorage storage(g, index.k(), edges, 10 * nthreads);

        BuildIndex(index, storage, nthreads);

        // Now use the index to fill the coverage and EdgeId's
        INFO("Collecting edge information from graph, this takes a while.");
        EdgeInfoUpdater<Graph>().Update(g, index, edges);
    }

    void BuildIndexFromGraph(Index &index, const Graph &g,
                             fs::TmpDir workdir, size_t read_buffer_size = 0) const {
        unsigned nthreads = omp_get_max_threads();

        using Splitter = DeBruijnGraphKMerSplitter<Graph,
                                                   utils::StoringTypeFilter<typename Index::storing_type>>;

        kmers::KMerDiskCounter<RtSeq> counter(workdir,
                                              Splitter(workdir, index.k(), g, read_buffer_size));

        BuildIndex(index, counter, 10 * nthreads, nthreads);

        // Now use the index to fill the coverage and EdgeId's
        INFO("Collecting edge information from graph, this takes a while.");
        EdgeInfoUpdater<Graph>().UpdateAll(g, index);
    }

    void BuildIndexFromGraph(Index &index,
                             const Graph &g, const std::vector<typename Graph::EdgeId> &edges,
                             fs::TmpDir workdir, size_t read_buffer_size = 0) const {
        unsigned nthreads = omp_get_max_threads();

        using Splitter = DeBruijnEdgeKMerSplitter<Graph,
                                                  utils::StoringTypeFilter<typename Index::storing_type>>;

        kmers::KMerDiskCounter<RtSeq> counter(workdir,
                                              Splitter(workdir, index.k(), g, edges, read_buffer_size));

        BuildIndex(index, counter, 10 * nthreads, nthreads);

        // Now use the index to fill the coverage and EdgeId's
        INFO("Update edge information.");
        EdgeInfoUpdater<Graph>().Update(g, index, edges);
    }

};

}
