//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/index/edge_index_builders.hpp"
#include "assembly_graph/core/graph.hpp"
#include "utils/filesystem/temporary.hpp"

#include "edge_index_refiller.hpp"

namespace debruijn_graph {

using EdgeIndex = KmerFreeEdgeIndex<ConjugateDeBruijnGraph>;
using EdgeIndex64 = KmerFreeEdgeIndex<ConjugateDeBruijnGraph, uint64_t>;
using EdgeIndex32 = KmerFreeEdgeIndex<ConjugateDeBruijnGraph, uint32_t>;

EdgeIndexRefiller::EdgeIndexRefiller(const std::string &workdir)
    : workdir_(workdir)
{}

template<>
void EdgeIndexRefiller::Refill(EdgeIndex &index,
                               const ConjugateDeBruijnGraph &g) {
    auto workdir = fs::tmp::make_temp_dir(workdir_, "edge_index");

    typedef GraphPositionFillingIndexBuilder<EdgeIndex> IndexBuilder;
    IndexBuilder().BuildIndexFromGraph(index, g, workdir);
}

template<>
void EdgeIndexRefiller::Refill(EdgeIndex64 &index,
                               const ConjugateDeBruijnGraph &g) {
    auto workdir = fs::tmp::make_temp_dir(workdir_, "edge_index");

    typedef GraphPositionFillingIndexBuilder<EdgeIndex64> IndexBuilder;
    IndexBuilder().BuildIndexFromGraph(index, g, workdir);
}

template<>
void EdgeIndexRefiller::Refill(EdgeIndex32 &index,
                               const ConjugateDeBruijnGraph &g) {
    auto workdir = fs::tmp::make_temp_dir(workdir_, "edge_index");

    typedef GraphPositionFillingIndexBuilder<EdgeIndex32> IndexBuilder;
    IndexBuilder().BuildIndexFromGraph(index, g, workdir);
}

}
