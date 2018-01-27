//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include <boost/test/unit_test.hpp>
#include "test_utils.hpp"
#include "assembly_graph/handlers/id_track_handler.hpp"

namespace debruijn_graph {
template<class Graph>
class IteratorOrderChecker {
private:
    const Graph &graph1_;
    const Graph &graph2_;
public:
    IteratorOrderChecker(const Graph &graph1, const Graph &graph2) :graph1_(graph1), graph2_(graph2) {
    }

    template<typename iterator>
    bool CheckOrder(iterator it1, iterator it2) {
        while(!it1.IsEnd() && !it2.IsEnd()) {
//            cout << graph1_.int_id(*it1) << " " << graph2_.int_id(*it2) << endl;
            if(graph1_.int_id(*it1) != graph2_.int_id(*it2))
                return false;
            ++it1;
            ++it2;
        }
        return it1.IsEnd() && it2.IsEnd();
    }
};

template<class Graph>
class RandomGraphConstructor {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    size_t iteration_number_;
    size_t max_size_;
    unsigned rand_seed_;

    Sequence GenerateRandomSequence(size_t length) {
        string result(length, 'A');
        for(size_t i = 0; i < length; i++) {
            result[i] = nucl(rand() % 4);
        }
        return Sequence(result);
    }

    void AddRandomVertex(Graph &graph) {
        graph.AddVertex();
    }

    VertexId GetRandomVertex(const Graph &graph) {
        size_t num = rand() % graph.size();
        auto it = graph.SmartVertexBegin();
        for(; num > 0; num--) {
            ++it;
        }
        return *it;
    }

    EdgeId GetRandomEdge(const Graph &graph) {
        EdgeId result = *(graph.SmartEdgeBegin());
        size_t cur = 0;
        for(auto it = graph.SmartEdgeBegin(); !it.IsEnd(); ++it) {
            cur++;
            if(rand() % cur == 0) {
                result = *it;
            }
        }
        return result;
    }

    void AddRandomEdge(Graph &graph) {
        graph.AddEdge(GetRandomVertex(graph), GetRandomVertex(graph), GenerateRandomSequence(rand() % 1000 + graph.k() + 1));
    }

    void RemoveRandomVertex(Graph &graph) {
        graph.ForceDeleteVertex(GetRandomVertex(graph));
    }

    void RemoveRandomEdge(Graph &graph) {
        graph.DeleteEdge(GetRandomEdge(graph));
    }

    void PerformRandomOperation(Graph &graph) {
        if (graph.size() == 0) {
            AddRandomVertex(graph);
        } else if (graph.SmartEdgeBegin().IsEnd()) {
            if (rand() % 2 == 0) {
                AddRandomVertex(graph);
            } else {
                AddRandomEdge(graph);
            }
        } else if (graph.size() > 100) {
            RemoveRandomVertex(graph);
        } else {
            size_t tmp = rand() % 9;
            if (tmp == 0) {
                AddRandomVertex(graph);
            } else if (tmp <= 6) {
                AddRandomEdge(graph);
            } else {
                RemoveRandomEdge(graph);
            }
        }
    }

public:
    RandomGraphConstructor(size_t iteration_number, size_t max_size, unsigned rand_seed = 100) :
            iteration_number_(iteration_number), max_size_(max_size), rand_seed_(rand_seed) {
    }

    void Generate(Graph &graph) {
        srand(rand_seed_);
        for (size_t i = 0; i < 1000; i++) {
            PerformRandomOperation(graph);
        }
    }
};

BOOST_FIXTURE_TEST_SUITE(robust_order_tests, fs::TmpFolderFixture)

BOOST_AUTO_TEST_CASE( OrderTest ) {
    string file_name = "src/test/debruijn/graph_fragments/saves/test_save";
    Graph graph(55);
    RandomGraphConstructor<Graph>(1000, 100, 100).Generate(graph);
    graphio::ConjugateDataPrinter<Graph> printer(graph);
    printer.SaveGraph(file_name);
    printer.SaveEdgeSequences(file_name);
    Graph new_graph(55);
    graphio::ConjugateDataScanner<Graph> scanner(new_graph);
    scanner.LoadGraph(file_name);
    IteratorOrderChecker<Graph> checker(graph, new_graph);
    BOOST_CHECK(checker.CheckOrder(graph.SmartVertexBegin(), new_graph.SmartVertexBegin()));
    BOOST_CHECK(checker.CheckOrder(graph.SmartEdgeBegin(), new_graph.SmartEdgeBegin()));
//    BOOST_CHECK(checker.CheckOrder(graph.SmartVertexBegin(), new_graph.SmartVertexBegin()));
//    BOOST_CHECK(checker.CheckOrder(graph.SmartEdgeBegin(), new_graph.SmartEdgeBegin()));
}
BOOST_AUTO_TEST_SUITE_END()
}
