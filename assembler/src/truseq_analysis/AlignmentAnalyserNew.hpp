#pragma once

#include <io/single_read.hpp>
#include "consistent_mapping.h"

namespace alignment_analysis {

    class AlignmentAnalyserNew {
    private:
        typedef debruijn_graph::DeBruijnGraph Graph;
        typedef Graph::EdgeId EdgeId;
        typedef Graph::VertexId VertexId;
    public:
        AlignmentAnalyserNew(Graph const &graph, size_t step) : graph_(graph), step_(step) { }
        vector <ConsistentMapping> Analyse(const MappingPath<EdgeId> &path) const;
    private:
        void Cut(vector<ConsistentMapping> &path, VertexId start) const;

        size_t StepBack(const vector<ConsistentMapping> &path) const;

        const Graph &graph_;
        size_t step_;
        DECL_LOGGER("AlignmentAnalyserNew")
    };

}
