//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "contracted_graph_builder.hpp"

#include "modules/path_extend/scaff_supplementary.hpp"

namespace contracted_graph {

class ContractedGraphFactoryHelper {
 public:
    typedef path_extend::ScaffoldingUniqueEdgeStorage UniqueStorage;
    typedef debruijn_graph::VertexId VertexId;
    typedef debruijn_graph::Graph Graph;

    explicit ContractedGraphFactoryHelper(const debruijn_graph::Graph &g) : g_(g) {}

    std::shared_ptr<ContractedGraph> ConstructFromUniqueStorage(const UniqueStorage &unique_storage) const;

    std::shared_ptr<ContractedGraph> ExtractContractedSubgraph(const ContractedGraph &other,
                                                               const std::unordered_set<VertexId> &vertices) const;

  private:
    const Graph& g_;
};
}