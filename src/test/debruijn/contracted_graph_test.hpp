#pragma once

#include "test_utils.hpp"
#include "auxiliary_graphs/contracted_graph/contracted_graph_helper.hpp"

#include <boost/test/unit_test.hpp>

namespace debruijn_graph {

std::shared_ptr<contracted_graph::ContractedGraph> CreateContractedGraph(const Graph& g, size_t length_threshold) {
    auto length_predicate = [length_threshold, &g](EdgeId edge) {
      return g.length(edge) >= length_threshold;
    };

    contracted_graph::DBGContractedGraphFactory graph_factory(g, length_predicate);
    graph_factory.Construct();
    auto contracted_graph = graph_factory.GetGraph();

    return contracted_graph;
}

BOOST_FIXTURE_TEST_SUITE(contracted_graph_tests, fs::TmpFolderFixture);

BOOST_AUTO_TEST_CASE( ContractedGraphBuilder ) {
    Graph g(55);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/contracted_graph/simple_bulge", g);
    BOOST_CHECK_EQUAL(g.size(), 8);
    BOOST_CHECK_EQUAL(g.k(), 55);

    size_t contraction_threshold = 80;
    auto contracted_graph = CreateContractedGraph(g, contraction_threshold);

    BOOST_CHECK_EQUAL(contracted_graph->size(), g.size());
    BOOST_CHECK_EQUAL(contracted_graph->CountEdges(), 10);
    omnigraph::IterationHelper<Graph, VertexId> vertex_it_helper(g);
    for (const auto& vertex: vertex_it_helper) {
        BOOST_CHECK(contracted_graph->ContainsVertex(vertex));
    }

    const size_t first_id = 132238743;
    const size_t second_id = 102357497;
    const size_t third_id = 243148631;
    const size_t fourth_id = 170924699;

    std::map<size_t, VertexId> simple_vertex_map;
    for (const auto& vertex: vertex_it_helper) {
        simple_vertex_map.insert({vertex.int_id(), vertex});
    }
    BOOST_CHECK_EQUAL(contracted_graph->GetOutDegree(simple_vertex_map.at(first_id)), 2);
    BOOST_CHECK_EQUAL(contracted_graph->GetInDegree(simple_vertex_map.at(first_id)), 0);
    BOOST_CHECK_EQUAL(contracted_graph->GetOutDegree(simple_vertex_map.at(second_id)), 1);
    BOOST_CHECK_EQUAL(contracted_graph->GetInDegree(simple_vertex_map.at(second_id)), 2);
    BOOST_CHECK_EQUAL(contracted_graph->GetOutDegree(simple_vertex_map.at(third_id)), 2);
    BOOST_CHECK_EQUAL(contracted_graph->GetInDegree(simple_vertex_map.at(third_id)), 1);
    BOOST_CHECK_EQUAL(contracted_graph->GetOutDegree(simple_vertex_map.at(fourth_id)), 0);
    BOOST_CHECK_EQUAL(contracted_graph->GetInDegree(simple_vertex_map.at(fourth_id)), 2);

    contraction_threshold = 200;
    contracted_graph = CreateContractedGraph(g, contraction_threshold);
    BOOST_CHECK_EQUAL(contracted_graph->size(), 6);
    BOOST_CHECK_EQUAL(contracted_graph->CountEdges(), 8);
    BOOST_CHECK_EQUAL(contracted_graph->GetOutDegree(simple_vertex_map.at(first_id)), 2);
    BOOST_CHECK_EQUAL(contracted_graph->GetInDegree(simple_vertex_map.at(first_id)), 0);
    BOOST_CHECK_EQUAL(contracted_graph->GetOutDegree(simple_vertex_map.at(fourth_id)), 0);
    BOOST_CHECK_EQUAL(contracted_graph->GetInDegree(simple_vertex_map.at(fourth_id)), 2);
    std::vector<scaffold_graph::ScaffoldVertex> out_edges;
    for (const auto &edge: contracted_graph->OutcomingEdges(simple_vertex_map.at(first_id))) {
        out_edges.emplace_back(edge);
    }
    BOOST_CHECK_EQUAL(out_edges.size(), 2);
    out_edges.clear();
    for (const auto &edge: contracted_graph->OutcomingEdges(simple_vertex_map.at(fourth_id))) {
        out_edges.emplace_back(edge);
    }
    BOOST_CHECK_EQUAL(out_edges.size(), 0);

    std::vector<scaffold_graph::ScaffoldVertex> in_edges;
    for (const auto &edge: contracted_graph->IncomingEdges(simple_vertex_map.at(first_id))) {
        in_edges.emplace_back(edge);
    }
    BOOST_CHECK_EQUAL(in_edges.size(), 0);
    in_edges.clear();
    for (const auto &edge: contracted_graph->IncomingEdges(simple_vertex_map.at(fourth_id))) {
        in_edges.emplace_back(edge);
    }
    BOOST_CHECK_EQUAL(in_edges.size(), 2);

    contraction_threshold = 5000;
    contracted_graph = CreateContractedGraph(g, contraction_threshold);
    BOOST_CHECK_EQUAL(contracted_graph->size(), 4);
    BOOST_CHECK_EQUAL(contracted_graph->CountEdges(), 4);
    BOOST_CHECK_EQUAL(contracted_graph->GetOutDegree(simple_vertex_map.at(first_id)), 2);
    BOOST_CHECK_EQUAL(contracted_graph->GetInDegree(simple_vertex_map.at(first_id)), 0);

    contraction_threshold = 25000;
    contracted_graph = std::move(CreateContractedGraph(g, contraction_threshold));
    BOOST_CHECK_EQUAL(contracted_graph->size(), 2);
    BOOST_CHECK_EQUAL(contracted_graph->CountEdges(), 0);
}

BOOST_AUTO_TEST_CASE( ContractedSubgraph ) {
    typedef contracted_graph::ContractedGraphFactoryHelper ContractedGraphFactoryHelper;
    Graph g(55);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/contracted_graph/simple_bulge", g);

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
    BOOST_CHECK_EQUAL(first_subgraph->size(), 3);
    BOOST_CHECK_EQUAL(first_subgraph->CountEdges(), 3);

    std::unordered_set<VertexId> first_third_subset({vertex_map.at(first_id), vertex_map.at(third_id)});
    auto second_subgraph = helper.ExtractContractedSubgraph(*contracted_graph, first_third_subset);
    BOOST_CHECK_EQUAL(second_subgraph->size(), 2);
    BOOST_CHECK_EQUAL(second_subgraph->CountEdges(), 0);

    std::unordered_set<VertexId> no_third_subset({vertex_map.at(first_id),
                                                  vertex_map.at(second_id),
                                                  vertex_map.at(fourth_id)});
    auto third_subgraph = helper.ExtractContractedSubgraph(*contracted_graph, no_third_subset);
    BOOST_CHECK_EQUAL(third_subgraph->size(), 3);
    BOOST_CHECK_EQUAL(third_subgraph->CountEdges(), 2);
}

BOOST_AUTO_TEST_SUITE_END()
}
