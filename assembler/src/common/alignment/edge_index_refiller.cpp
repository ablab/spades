//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/core/kmer_iterator.hpp"
#include "assembly_graph/index/edge_index_builders.hpp"
#include "utils/filesystem/temporary.hpp"

#include "edge_index_refiller.hpp"

namespace debruijn_graph {

using EdgeIndex = KmerFreeEdgeIndex<ConjugateDeBruijnGraph>;
using EdgeIndex64 = KmerFreeEdgeIndex<ConjugateDeBruijnGraph, uint64_t>;
using EdgeIndex32 = KmerFreeEdgeIndex<ConjugateDeBruijnGraph, uint32_t>;

EdgeIndexRefiller::EdgeIndexRefiller(const std::string &workdir)
    : workdir_(workdir)
{}

template<class EdgeIndex>
void EdgeIndexRefiller::Refill(EdgeIndex &index, const Graph &g) {
    typedef GraphPositionFillingIndexBuilder<EdgeIndex> IndexBuilder;
    IndexBuilder().BuildIndexFromGraph(index, g);
}

template
void EdgeIndexRefiller::Refill(EdgeIndex &index, const Graph &g);

template
void EdgeIndexRefiller::Refill(EdgeIndex64 &index, const Graph &g);

template
void EdgeIndexRefiller::Refill(EdgeIndex32 &index, const Graph &g);


template<class EdgeIndex>
void EdgeIndexRefiller::Refill(EdgeIndex &index,
                               const Graph &g,
                               const std::vector<EdgeId> &edges) {
    typedef GraphPositionFillingIndexBuilder<EdgeIndex> IndexBuilder;
    IndexBuilder().BuildIndexFromGraph(index, g, edges);
}

template
void EdgeIndexRefiller::Refill(EdgeIndex &index,
                               const Graph &g,
                               const std::vector<typename ConjugateDeBruijnGraph::EdgeId> &edges);

template
void EdgeIndexRefiller::Refill(EdgeIndex64 &index,
                               const ConjugateDeBruijnGraph &g,
                               const std::vector<typename ConjugateDeBruijnGraph::EdgeId> &edges);

template
void EdgeIndexRefiller::Refill(EdgeIndex32 &index,
                               const ConjugateDeBruijnGraph &g,
                               const std::vector<typename ConjugateDeBruijnGraph::EdgeId> &edges);
}
