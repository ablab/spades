#include "contracted_graph_helper.hpp"

namespace contracted_graph {

ContractedGraph ContractedGraphFactoryHelper::ConstructFromUniqueStorage(const ContractedGraphFactoryHelper::UniqueStorage& unique_storage) const {
    std::function<bool(EdgeId)> edge_predicate = [&unique_storage](EdgeId edge) {
      return unique_storage.IsUnique(edge);
    };
    DBGContractedGraphFactory factory(g_, edge_predicate);
    factory.Construct();
    return *(factory.GetGraph());
}
ContractedGraph ContractedGraphFactoryHelper::ConstructFromInternalGraph(
    const cluster_storage::Cluster::InternalGraph& internal_graph) const {
    SimpleContractedGraphFactory factory(g_, internal_graph);
    factory.Construct();
    return *(factory.GetGraph());
}
ContractedGraph ContractedGraphFactoryHelper::ExtractContractedSubgraph(const ContractedGraph& other,
                                                                        const std::unordered_set<VertexId>& vertices) const {
    SubgraphContractedGraphFactory factory(other, vertices);
    factory.Construct();
    return *(factory.GetGraph());
}
}