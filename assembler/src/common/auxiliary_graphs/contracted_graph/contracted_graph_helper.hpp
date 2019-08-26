//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "contracted_graph_builder.hpp"

#include "assembly_graph/graph_support/scaff_supplementary.hpp"

namespace contracted_graph {

class ContractedGraphFactoryHelper {
 public:
    typedef path_extend::ScaffoldingUniqueEdgeStorage UniqueStorage;
    typedef path_extend::read_cloud::SimpleGraph<scaffold_graph::ScaffoldVertex> SimpleGraph;

    explicit ContractedGraphFactoryHelper(const Graph &g) : g_(g) {}

    std::shared_ptr<ContractedGraph> ConstructFromUniqueStorage(const UniqueStorage &unique_storage) const;

    std::shared_ptr<ContractedGraph> ExtractContractedSubgraph(const ContractedGraph &other,
                                                               const std::unordered_set<VertexId> &vertices) const;

  private:
    const Graph& g_;
};
}