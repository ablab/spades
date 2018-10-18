//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "index/edge_index_builders.hpp"
#include "index/edge_multi_index.hpp"
#include "core/graph.hpp"
#include "utils/filesystem/temporary.hpp"

#include "edge_index_refiller.hpp"

namespace debruijn_graph {

using EdgeIndex = KmerFreeEdgeIndex<ConjugateDeBruijnGraph>;

template<>
void EdgeIndexRefiller::Refill(EdgeIndex &index,
                               const ConjugateDeBruijnGraph &g) {
    auto workdir = fs::tmp::make_temp_dir(workdir_, "edge_index");

    typedef typename EdgeIndexHelper<EdgeIndex>::GraphPositionFillingIndexBuilderT IndexBuilder;
    IndexBuilder().BuildIndexFromGraph(index, g, workdir);
}

using PacIndex = DeBruijnEdgeMultiIndex<ConjugateDeBruijnGraph::EdgeId>;

template<>
void EdgeIndexRefiller::Refill(PacIndex &index,
                               const ConjugateDeBruijnGraph &g) {
    auto workdir = fs::tmp::make_temp_dir(workdir_, "edge_index");

    typedef typename debruijn_graph::EdgeIndexHelper<PacIndex>::GraphPositionFillingIndexBuilderT Builder;
    Builder().BuildIndexFromGraph(index, g, workdir);
}

}
