//***************************************************************************
//* Copyright (c) 2023-2025 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "contracted_graph.hpp"

#include "adt/concurrent_dsu.hpp"

#include <memory>

namespace contracted_graph {

class ContractedGraphFactory {
 public:
    using Graph = debruijn_graph::Graph;

    ContractedGraphFactory(const Graph &g) :
        g_(g), graph_ptr_(std::make_shared<ContractedGraph>(g)) {}
    virtual ~ContractedGraphFactory() = default;
    virtual void Construct() = 0;
    std::shared_ptr<ContractedGraph> GetGraph() {
        return graph_ptr_;
    }
 protected:
    const Graph &g_;
    std::shared_ptr<ContractedGraph> graph_ptr_;
};

class PartsBasedContractedFactory : public ContractedGraphFactory {
 public:
    using VertexId = debruijn_graph::VertexId;
    using ContractedGraphFactory::Graph;

    PartsBasedContractedFactory(const Graph &g): ContractedGraphFactory(g) {}
    virtual ~PartsBasedContractedFactory() {}

    void Construct() override;
  protected:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

    struct ContractedGraphParts {
      std::vector <ScaffoldVertex> long_edges_;
      std::unordered_set <VertexId> long_edge_ends_;
      std::unordered_map <VertexId, size_t> vertex_to_capacity_;
      std::unordered_map <VertexId, VertexId> vertex_to_root_;
    };

    virtual ContractedGraphParts ConstructParts() const = 0;
    void ConstructFromParts(ContractedGraphParts &&parts);

    using ContractedGraphFactory::graph_ptr_;
    using ContractedGraphFactory::g_;
    DECL_LOGGER("DSUBasedContractedGraphFactory");
};

class DBGContractedGraphFactory : public PartsBasedContractedFactory {
  public:
    using contracted_dsu_t = dsu::ConcurrentDSU;
    using EdgeId = debruijn_graph::EdgeId;
    using PartsBasedContractedFactory::Graph;

    DBGContractedGraphFactory(const Graph &g, const std::function<bool(EdgeId)> &edge_predicate) :
        PartsBasedContractedFactory(g), edge_predicate_(edge_predicate) {}

 private:
    ContractedGraphParts ConstructParts() const override;

    void ProcessEdge(contracted_dsu_t &graph_dsu, ContractedGraphParts &parts,
                     const std::unordered_map<VertexId, size_t> &vertex_to_id,
                     const std::unordered_map<size_t, VertexId> &id_to_vertex, EdgeId edge) const;

    using PartsBasedContractedFactory::g_;
    using PartsBasedContractedFactory::graph_ptr_;
    using PartsBasedContractedFactory::ContractedGraphParts;
    const std::function<bool(EdgeId)> edge_predicate_;

    DECL_LOGGER("DBGContractedGraphFactory");
};

class SubgraphContractedGraphFactory: public ContractedGraphFactory {
 public:
    using VertexId = debruijn_graph::VertexId;

    SubgraphContractedGraphFactory(const ContractedGraph &other, const std::unordered_set<VertexId> &vertices) :
        ContractedGraphFactory(other.GetAssemblyGraph()), other_(other), vertices_(vertices) {}

    void Construct() override;

 private:
    void ExtractSubgraphFromContractedGraph(const ContractedGraph &other, const std::unordered_set<VertexId> &vertices);

    const ContractedGraph &other_;
    const std::unordered_set<VertexId> &vertices_;
};

} //contracted_graph
