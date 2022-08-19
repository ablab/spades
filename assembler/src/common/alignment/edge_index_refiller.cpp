//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#include "edge_index_refiller.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/index/edge_index_builders.hpp"

#include <vector>

namespace debruijn_graph {

using EdgeIndex = KmerFreeEdgeIndex<ConjugateDeBruijnGraph>;
using EdgeIndex64 = KmerFreeEdgeIndex<ConjugateDeBruijnGraph, uint64_t>;
using EdgeIndex32 = KmerFreeEdgeIndex<ConjugateDeBruijnGraph, uint32_t>;

EdgeIndexRefiller::EdgeIndexRefiller(const std::filesystem::path &workdir)
    : workdir_(workdir)
{}

template<class EdgeIndex>
void EdgeIndexRefiller::Refill(EdgeIndex &index, const Graph &g,
                               bool count) {
    typedef GraphPositionFillingIndexBuilder<EdgeIndex> IndexBuilder;
    if (count) {
        IndexBuilder().BuildIndexFromGraph(index, g, fs::tmp::make_temp_dir(workdir_, "edge_index"));
    } else {
        IndexBuilder().BuildIndexFromGraph(index, g);
    }
}

template
void EdgeIndexRefiller::Refill(EdgeIndex &index, const Graph &g, bool);

template
void EdgeIndexRefiller::Refill(EdgeIndex64 &index, const Graph &g, bool);

template
void EdgeIndexRefiller::Refill(EdgeIndex32 &index, const Graph &g, bool);


template<class EdgeIndex>
void EdgeIndexRefiller::Refill(EdgeIndex &index,
                               const Graph &g,
                               const std::vector<EdgeId> &edges,
                               bool count) {
    typedef GraphPositionFillingIndexBuilder<EdgeIndex> IndexBuilder;
    if (count) {
        IndexBuilder().BuildIndexFromGraph(index, g, edges, fs::tmp::make_temp_dir(workdir_, "edge_index"));
    } else {
        IndexBuilder().BuildIndexFromGraph(index, g, edges);
    }
}

template
void EdgeIndexRefiller::Refill(EdgeIndex &index,
                               const Graph &g,
                               const std::vector<typename ConjugateDeBruijnGraph::EdgeId> &edges,
                               bool);

template
void EdgeIndexRefiller::Refill(EdgeIndex64 &index,
                               const ConjugateDeBruijnGraph &g,
                               const std::vector<typename ConjugateDeBruijnGraph::EdgeId> &edges,
                               bool);

template
void EdgeIndexRefiller::Refill(EdgeIndex32 &index,
                               const ConjugateDeBruijnGraph &g,
                               const std::vector<typename ConjugateDeBruijnGraph::EdgeId> &edges,
                               bool);
}
