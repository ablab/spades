//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "contracted_graph_helper.hpp"

namespace contracted_graph {

std::shared_ptr<ContractedGraph> ContractedGraphFactoryHelper::ConstructFromUniqueStorage(
        const UniqueStorage &unique_storage) const {
    std::function<bool(debruijn_graph::EdgeId)> edge_predicate = [&unique_storage](debruijn_graph::EdgeId edge) {
      return unique_storage.IsUnique(edge);
    };
    DBGContractedGraphFactory factory(g_, edge_predicate);
    factory.Construct();
    return factory.GetGraph();
}
std::shared_ptr<ContractedGraph> ContractedGraphFactoryHelper::ExtractContractedSubgraph(
        const ContractedGraph &other,
        const std::unordered_set<VertexId> &vertices) const {
    SubgraphContractedGraphFactory factory(other, vertices);
    factory.Construct();
    return factory.GetGraph();
}
}