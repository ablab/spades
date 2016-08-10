//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include "consistent_mapping.h"

namespace alignment_analysis {

    class AlignmentAnalyserNew {
    private:
        typedef debruijn_graph::DeBruijnGraph Graph;
        typedef Graph::EdgeId EdgeId;
        typedef Graph::VertexId VertexId;
    public:
        AlignmentAnalyserNew(Graph const &graph, size_t step) : graph_(graph), step_(step) { }
        vector <ConsistentMapping> Analyse(const omnigraph::MappingPath<EdgeId> &path) const;
    private:
        void Cut(vector<ConsistentMapping> &path, VertexId start) const;

        size_t StepBack(const vector<ConsistentMapping> &path) const;

        const Graph &graph_;
        size_t step_;
        DECL_LOGGER("AlignmentAnalyserNew")
    };

}
