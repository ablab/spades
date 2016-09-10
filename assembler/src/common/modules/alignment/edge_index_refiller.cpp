//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/indices/edge_index_builders.hpp"
#include "utils/indices/edge_multi_index.hpp"
#include "core/graph.hpp"

#include "edge_index_refiller.hpp"

namespace debruijn_graph {

using EdgeIndex = KmerFreeEdgeIndex<ConjugateDeBruijnGraph>;

template<>
void EdgeIndexRefiller::Refill(EdgeIndex &index,
                               const ConjugateDeBruijnGraph &g) {
    typedef typename EdgeIndexHelper<EdgeIndex>::GraphPositionFillingIndexBuilderT IndexBuilder;
    IndexBuilder().BuildIndexFromGraph(index, g);
}

using PacIndex = DeBruijnEdgeMultiIndex<ConjugateDeBruijnGraph::EdgeId>;

template<>
void EdgeIndexRefiller::Refill(PacIndex &index,
                               const ConjugateDeBruijnGraph &g) {
    typedef typename debruijn_graph::EdgeIndexHelper<PacIndex>::GraphPositionFillingIndexBuilderT Builder;
    Builder().BuildIndexFromGraph(index, g);
}

}
