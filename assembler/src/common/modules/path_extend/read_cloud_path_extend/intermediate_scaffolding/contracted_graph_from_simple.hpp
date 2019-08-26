//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "simple_graph.hpp"
#include "auxiliary_graphs/contracted_graph/contracted_graph_builder.hpp"

namespace path_extend {
namespace read_cloud {

class SimpleContractedGraphFactory : public contracted_graph::PartsBasedContractedFactory {
  public:
    typedef dsu::ConcurrentDSU contracted_dsu_t;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef path_extend::read_cloud::SimpleGraph<ScaffoldVertex> SimpleGraph;

    SimpleContractedGraphFactory(const Graph &assembly_graph, const SimpleGraph &simple_graph)
        : PartsBasedContractedFactory(assembly_graph), simple_graph_(simple_graph) {}

  private:
    ContractedGraphParts ConstructParts() const override;

    using PartsBasedContractedFactory::g_;
    using PartsBasedContractedFactory::graph_ptr_;
    using PartsBasedContractedFactory::ContractedGraphParts;
    const SimpleGraph &simple_graph_;

    DECL_LOGGER("SimpleContractedGraphFactory");
};

class ContractedGraphFromSimpleHelper {
  public:
    typedef SimpleGraph<scaffold_graph::ScaffoldVertex> ScaffoldVertexSimpleGraph;

    explicit ContractedGraphFromSimpleHelper(const Graph &g) : g_(g) {}

    std::shared_ptr<contracted_graph::ContractedGraph> ConstructFromSimpleGraph(
        const ScaffoldVertexSimpleGraph &simple_graph) const;
  private:
    const Graph& g_;
};

}

}