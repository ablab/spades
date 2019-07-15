#pragma once
#include "common/assembly_graph/contracted_graph/contracted_graph.hpp"
#include "adt/concurrent_dsu.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/cluster_storage.hpp"

namespace contracted_graph {

class ContractedGraphFactory {
 public:
    ContractedGraphFactory(const Graph &g) :
        g_(g), graph_ptr_(make_shared<ContractedGraph>(g)) {}
    virtual ~ContractedGraphFactory() = default;
    virtual void Construct() = 0;
    shared_ptr<ContractedGraph> GetGraph() {
        return graph_ptr_;
    }
 protected:
    const Graph &g_;
    shared_ptr<ContractedGraph> graph_ptr_;
};

class PartsBasedContractedFactory : public ContractedGraphFactory{

 protected:
    using ContractedGraphFactory::graph_ptr_;
    using ContractedGraphFactory::g_;
    typedef path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

    struct ContractedGraphParts {
      vector <ScaffoldVertex> long_edges_;
      unordered_set <VertexId> long_edge_ends_;
      unordered_map <VertexId, size_t> vertex_to_capacity_;
      unordered_map <VertexId, VertexId> vertex_to_root_;
    };

    virtual ContractedGraphParts ConstructParts() const = 0;

    void ConstructFromParts(ContractedGraphParts&& parts);
    DECL_LOGGER("DSUBasedContractedGraphFactory");
 public:
    PartsBasedContractedFactory(const Graph &g): ContractedGraphFactory(g) {}
    virtual ~PartsBasedContractedFactory() {}

    void Construct() override;
};

class DBGContractedGraphFactory : public PartsBasedContractedFactory {
    using PartsBasedContractedFactory::g_;
    using PartsBasedContractedFactory::graph_ptr_;
    using PartsBasedContractedFactory::ContractedGraphParts;
    typedef dsu::ConcurrentDSU contracted_dsu_t;

    const std::function<bool(EdgeId)> edge_predicate_;

 public:

    DBGContractedGraphFactory(const Graph& g, const std::function<bool(EdgeId)>& edge_predicate) :
        PartsBasedContractedFactory(g), edge_predicate_(edge_predicate) {}

 private:
    ContractedGraphParts ConstructParts() const override;

    void ProcessEdge(contracted_dsu_t& graph_dsu, ContractedGraphParts& parts,
                     const std::unordered_map<VertexId, size_t>& vertex_to_id,
                     const std::unordered_map<size_t, VertexId>& id_to_vertex, const EdgeId& edge) const;

    DECL_LOGGER("DBGContractedGraphFactory");
};

class SubgraphContractedGraphFactory: public ContractedGraphFactory {
    const ContractedGraph& other_;
    const std::unordered_set<VertexId>& vertices_;
 public:
    SubgraphContractedGraphFactory(const ContractedGraph& other, const std::unordered_set<VertexId>& vertices) :
        ContractedGraphFactory(other.GetAssemblyGraph()), other_(other), vertices_(vertices) {}

    void Construct() override;

 private:
    void ExtractSubgraphFromContractedGraph(const ContractedGraph& other, const std::unordered_set<VertexId>& vertices);
};

class SimpleContractedGraphFactory: public PartsBasedContractedFactory {
    using PartsBasedContractedFactory::g_;
    using PartsBasedContractedFactory::graph_ptr_;
    using PartsBasedContractedFactory::ContractedGraphParts;
    typedef dsu::ConcurrentDSU contracted_dsu_t;
    typedef path_extend::read_cloud::cluster_storage::Cluster::InternalGraph InternalGraph;
    typedef path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
    const InternalGraph& internal_graph_;

 public:
    SimpleContractedGraphFactory(const Graph& assembly_graph_, const InternalGraph& internal_graph_)
        : PartsBasedContractedFactory(assembly_graph_), internal_graph_(internal_graph_) {}

 private:
    ContractedGraphParts ConstructParts() const override;
    DECL_LOGGER("SimpleContractedGraphFactory");
};

} //contracted_graph