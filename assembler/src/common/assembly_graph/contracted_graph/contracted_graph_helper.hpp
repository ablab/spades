//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/assembly_graph/contracted_graph/contracted_graph_builder.hpp"
#include "common/assembly_graph/graph_support/scaff_supplementary.hpp"

namespace contracted_graph {

class ContractedGraphFactoryHelper {
 public:
    typedef path_extend::ScaffoldingUniqueEdgeStorage UniqueStorage;
    typedef path_extend::read_cloud::cluster_storage::Cluster::InternalGraph InternalGraph;

    explicit ContractedGraphFactoryHelper(const Graph& g_) : g_(g_) {}

    ContractedGraph ConstructFromUniqueStorage(const UniqueStorage& unique_storage) const;

    //todo move this somewhere (or better yet, get rid of SimpleGraphT)
    ContractedGraph ConstructFromInternalGraph(const InternalGraph& internal_graph) const;

    ContractedGraph ExtractContractedSubgraph(const ContractedGraph& other,
                                              const std::unordered_set<VertexId>& vertices) const;

  private:
    const Graph& g_;
};
}