#pragma once

#include "assembly_graph/core/coverage.hpp"
#include "assembly_graph/graph_support/detail_coverage.hpp"

namespace debruijn_graph {

template<class StoringType>
struct SimultaneousCoverageCollector {
};

template<>
struct SimultaneousCoverageCollector<SimpleStoring> {
    template<class SimultaneousCoverageFiller, class Info>
    static void CollectCoverage(SimultaneousCoverageFiller& filler, const Info &edge_info) {
        filler.inc_coverage(edge_info);
    }
};

template<>
struct SimultaneousCoverageCollector<InvertableStoring> {
    template<class SimultaneousCoverageFiller, class Info>
    static void CollectCoverage(SimultaneousCoverageFiller& filler, const Info &edge_info) {
        filler.inc_coverage(edge_info);
        filler.inc_coverage(edge_info.conjugate(filler.k()));
    }
};

template<class Graph, class CountIndex>
class SimultaneousCoverageFiller {
    const Graph& g_;
    const CountIndex& count_index_;
    omnigraph::FlankingCoverage<Graph>& flanking_coverage_;
    omnigraph::CoverageIndex<Graph>& coverage_index_;
    typedef typename CountIndex::KmerPos Value;
public:
    SimultaneousCoverageFiller(const Graph& g, const CountIndex& count_index,
                               omnigraph::FlankingCoverage<Graph>& flanking_coverage,
                               omnigraph::CoverageIndex<Graph>& coverage_index) :
            g_(g),
            count_index_(count_index),
            flanking_coverage_(flanking_coverage),
            coverage_index_(coverage_index) {
    }

    size_t k() const {
        return count_index_.k();
    }

    void inc_coverage(const Value &edge_info) {
        coverage_index_.IncRawCoverage(edge_info.edge_id, edge_info.count);
        if (edge_info.offset < flanking_coverage_.averaging_range()) {
            flanking_coverage_.IncRawCoverage(edge_info.edge_id, edge_info.count);
        }
    }

    void Fill() {
        for (auto I = count_index_.value_cbegin(), E = count_index_.value_cend();
             I != E; ++I) {
            const auto& edge_info = *I;
            //VERIFY(edge_info.valid());
            if (edge_info.valid()) {
                VERIFY(edge_info.edge_id.get() != NULL);
                SimultaneousCoverageCollector<typename CountIndex::storing_type>::CollectCoverage(*this, edge_info);
            } else {
                VERIFY(edge_info.removed());
                WARN("Duplicating k+1-mers in graph (known bug in construction)");
            }
        }
    }
};

template<class Graph, class CountIndex>
void FillCoverageAndFlanking(const CountIndex& count_index, Graph& g,
                             FlankingCoverage<Graph>& flanking_coverage) {
    SimultaneousCoverageFiller<Graph, CountIndex> filler(g, count_index, flanking_coverage, g.coverage_index());
    filler.Fill();
}

}