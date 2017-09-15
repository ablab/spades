#pragma once

#include "common/assembly_graph/contracted_graph/contracted_graph_builder.hpp"

namespace contracted_graph {

class ContractedGraphFactoryHelper {
 public:
    typedef path_extend::ScaffoldingUniqueEdgeStorage UniqueStorage;
 private:
    const Graph& g_;
 public:
    explicit ContractedGraphFactoryHelper(const Graph& g_) : g_(g_) {}

    ContractedGraph ConstructFromUniqueStorage(const UniqueStorage& unique_storage) const;

    ContractedGraph ConstructFromInternalGraph(const cluster_storage::Cluster::InternalGraph& internal_graph) const;

    ContractedGraph TransposeContractedGraph(const ContractedGraph& other) const;

    ContractedGraph ExtractContractedSubgraph(const ContractedGraph& other,
                                              const std::unordered_set<VertexId>& vertices) const;
};
}