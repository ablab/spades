//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************


#pragma once

#include "assembly_graph/core/coverage.hpp"
#include "assembly_graph/graph_support/detail_coverage.hpp"
#include "utils/parallel/openmp_wrapper.h"

namespace debruijn_graph {

template<class Graph, class PHM>
class GraphCoverageFiller {
    typedef typename omnigraph::GraphEdgeIterator<Graph> EdgeIt;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename adt::iterator_range<EdgeIt> EdgeRange;

    const Graph& g_;
    const PHM& phm_;
    omnigraph::FlankingCoverage<Graph>& flanking_coverage_;
    omnigraph::CoverageIndex<Graph>& coverage_index_;
    unsigned k_;
    size_t avg_range_;
  public:
    GraphCoverageFiller(const Graph& g, unsigned k, const PHM& phm,
                        omnigraph::FlankingCoverage<Graph>& flanking_coverage,
                        omnigraph::CoverageIndex<Graph>& coverage_index)
            : g_(g),
              phm_(phm),
              flanking_coverage_(flanking_coverage),
              coverage_index_(coverage_index),
              k_(k), avg_range_(flanking_coverage_.averaging_range()) {}

    void inc_coverage(EdgeId edge_id, size_t offset, uint32_t value) {
        coverage_index_.IncRawCoverage(edge_id, value);
        if (offset < avg_range_)
            flanking_coverage_.IncRawCoverage(edge_id, value);
    }

    size_t FillCoverageFromEdges(EdgeRange &r) {
        size_t seqs = 0;
        for (auto &it = r.begin(); it != r.end() && seqs < 100000; ++it) {
            EdgeId e = *it;
            const Sequence &seq = g_.EdgeNucls(e);
                
            seqs += 1;
            RtSeq kmer = seq.start<RtSeq>(this->k_) >> 'A';
            for (size_t j = this->k_ - 1; j < seq.size(); ++j) {
                kmer <<= seq[j];

                auto kwh = phm_.ConstructKWH(kmer);
                uint32_t cov = phm_.get_value(kwh, utils::InvertableStoring::trivial_inverter());
                inc_coverage(e, j - this->k_ + 1, cov);
            }
        }
        
        return seqs;
    }


    void Fill(unsigned nthreads) {
        omnigraph::IterationHelper<Graph, EdgeId> edges(g_);
        auto its = edges.Chunks(10*nthreads);

        // Turn chunks into iterator ranges
        std::vector<EdgeRange> ranges;
        for (size_t i = 0; i < its.size() - 1; ++i)
            ranges.emplace_back(its[i], its[i+1]);

        size_t counter = 0, n = 10;
        while (!std::all_of(ranges.begin(), ranges.end(),
                            [](const EdgeRange &r) { return r.begin() == r.end(); })) {
#               pragma omp parallel for num_threads(nthreads) reduction(+ : counter) schedule(guided)
            for (size_t i = 0; i < ranges.size(); ++i)
                counter += FillCoverageFromEdges(ranges[i]);

            if (counter >> n) {
                INFO("Processed " << counter << " edges");
                n += 1;
            }
        }
    }
};

template<class Graph, class PHM>
void FillCoverageAndFlankingFromPHM(const PHM& phm, Graph& g,
                                    omnigraph::FlankingCoverage<Graph>& flanking_coverage) {
    GraphCoverageFiller<Graph, PHM>(g,
                                    unsigned(g.k() + 1), phm,
                                    flanking_coverage, g.coverage_index()).Fill(omp_get_max_threads());
}

}
