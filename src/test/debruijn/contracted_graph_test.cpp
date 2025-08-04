//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2023-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "auxiliary_graphs/contracted_graph/contracted_graph_helper.hpp"
#include "graphio.hpp"

#include <filesystem>
#include <gtest/gtest.h>

using namespace debruijn_graph;

std::shared_ptr<contracted_graph::ContractedGraph> CreateContractedGraph(const Graph& g, size_t length_threshold) {
    auto length_predicate = [length_threshold, &g](EdgeId edge) {
      return g.length(edge) >= length_threshold;
    };

    contracted_graph::DBGContractedGraphFactory graph_factory(g, length_predicate);
    graph_factory.Construct();
    auto contracted_graph = graph_factory.GetGraph();

    return contracted_graph;
}

TEST(ContractedGraph, Construction) {
    Graph g(55);
    std::string save_path = "src/test/debruijn/graph_fragments/contracted_graph/simple_bulge";
    graphio::ScanBasicGraph(save_path, g);
    EXPECT_EQ(g.size(), 8);
    EXPECT_EQ(g.k(), 55);

    size_t contraction_threshold = 80;
    auto contracted_graph = CreateContractedGraph(g, contraction_threshold);

    EXPECT_EQ(contracted_graph->size(), g.size());
    EXPECT_EQ(contracted_graph->CountEdges(), 10);
    omnigraph::IterationHelper<Graph, VertexId> vertex_it_helper(g);
    for (const auto& vertex: vertex_it_helper) {
        EXPECT_TRUE(contracted_graph->ContainsVertex(vertex));
    }

    const size_t first_id = 132238743;
    const size_t second_id = 102357497;
    const size_t third_id = 243148631;
    const size_t fourth_id = 170924699;

    std::map<size_t, VertexId> simple_vertex_map;
    for (const auto& vertex: vertex_it_helper) {
        simple_vertex_map.insert({vertex.int_id(), vertex});
    }
    EXPECT_EQ(contracted_graph->GetOutDegree(simple_vertex_map.at(first_id)), 2);
    EXPECT_EQ(contracted_graph->GetInDegree(simple_vertex_map.at(first_id)), 0);
    EXPECT_EQ(contracted_graph->GetOutDegree(simple_vertex_map.at(second_id)), 1);
    EXPECT_EQ(contracted_graph->GetInDegree(simple_vertex_map.at(second_id)), 2);
    EXPECT_EQ(contracted_graph->GetOutDegree(simple_vertex_map.at(third_id)), 2);
    EXPECT_EQ(contracted_graph->GetInDegree(simple_vertex_map.at(third_id)), 1);
    EXPECT_EQ(contracted_graph->GetOutDegree(simple_vertex_map.at(fourth_id)), 0);
    EXPECT_EQ(contracted_graph->GetInDegree(simple_vertex_map.at(fourth_id)), 2);

    contraction_threshold = 200;
    contracted_graph = CreateContractedGraph(g, contraction_threshold);
    EXPECT_EQ(contracted_graph->size(), 6);
    EXPECT_EQ(contracted_graph->CountEdges(), 8);
    EXPECT_EQ(contracted_graph->GetOutDegree(simple_vertex_map.at(first_id)), 2);
    EXPECT_EQ(contracted_graph->GetInDegree(simple_vertex_map.at(first_id)), 0);
    EXPECT_EQ(contracted_graph->GetOutDegree(simple_vertex_map.at(fourth_id)), 0);
    EXPECT_EQ(contracted_graph->GetInDegree(simple_vertex_map.at(fourth_id)), 2);
    std::vector<scaffold_graph::ScaffoldVertex> out_edges;
    for (const auto &edge: contracted_graph->OutgoingEdges(simple_vertex_map.at(first_id))) {
        out_edges.emplace_back(edge);
    }
    EXPECT_EQ(out_edges.size(), 2);
    out_edges.clear();
    for (const auto &edge: contracted_graph->OutgoingEdges(simple_vertex_map.at(fourth_id))) {
        out_edges.emplace_back(edge);
    }
    EXPECT_EQ(out_edges.size(), 0);

    std::vector<scaffold_graph::ScaffoldVertex> in_edges;
    for (const auto &edge: contracted_graph->IncomingEdges(simple_vertex_map.at(first_id))) {
        in_edges.emplace_back(edge);
    }
    EXPECT_EQ(in_edges.size(), 0);
    in_edges.clear();
    for (const auto &edge: contracted_graph->IncomingEdges(simple_vertex_map.at(fourth_id))) {
        in_edges.emplace_back(edge);
    }
    EXPECT_EQ(in_edges.size(), 2);

    contraction_threshold = 5000;
    contracted_graph = CreateContractedGraph(g, contraction_threshold);
    EXPECT_EQ(contracted_graph->size(), 4);
    EXPECT_EQ(contracted_graph->CountEdges(), 4);
    EXPECT_EQ(contracted_graph->GetOutDegree(simple_vertex_map.at(first_id)), 2);
    EXPECT_EQ(contracted_graph->GetInDegree(simple_vertex_map.at(first_id)), 0);

    contraction_threshold = 25000;
    contracted_graph = std::move(CreateContractedGraph(g, contraction_threshold));
    EXPECT_EQ(contracted_graph->size(), 2);
    EXPECT_EQ(contracted_graph->CountEdges(), 0);
}

TEST(ContractedGraph, Subgraph) {
    typedef contracted_graph::ContractedGraphFactoryHelper ContractedGraphFactoryHelper;
    Graph g(55);
    std::string save_path = "./src/test/debruijn/graph_fragments/contracted_graph/simple_bulge";
    graphio::ScanBasicGraph(save_path, g);

    size_t contraction_threshold = 80;
    auto contracted_graph = CreateContractedGraph(g, contraction_threshold);

    const size_t first_id = 132238743;
    const size_t second_id = 102357497;
    const size_t third_id = 243148631;
    const size_t fourth_id = 170924699;

    omnigraph::IterationHelper<Graph, VertexId> vertex_it_helper(g);

    std::map<size_t, VertexId> vertex_map;
    for (const auto& vertex: vertex_it_helper) {
        vertex_map.insert({vertex.int_id(), vertex});
    }
    ContractedGraphFactoryHelper helper(g);

    std::unordered_set<VertexId> no_fourth_subset({vertex_map.at(first_id),
                                                   vertex_map.at(second_id),
                                                   vertex_map.at(third_id)});
    auto first_subgraph = helper.ExtractContractedSubgraph(*contracted_graph, no_fourth_subset);
    EXPECT_EQ(first_subgraph->size(), 3);
    EXPECT_EQ(first_subgraph->CountEdges(), 3);

    std::unordered_set<VertexId> first_third_subset({vertex_map.at(first_id), vertex_map.at(third_id)});
    auto second_subgraph = helper.ExtractContractedSubgraph(*contracted_graph, first_third_subset);
    EXPECT_EQ(second_subgraph->size(), 2);
    EXPECT_EQ(second_subgraph->CountEdges(), 0);

    std::unordered_set<VertexId> no_third_subset({vertex_map.at(first_id),
                                                  vertex_map.at(second_id),
                                                  vertex_map.at(fourth_id)});
    auto third_subgraph = helper.ExtractContractedSubgraph(*contracted_graph, no_third_subset);
    EXPECT_EQ(third_subgraph->size(), 3);
    EXPECT_EQ(third_subgraph->CountEdges(), 2);
}
