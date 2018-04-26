//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include <boost/test/unit_test.hpp>
#include "test_utils.hpp"
#include "random_graph.hpp"
#include "assembly_graph/handlers/id_track_handler.hpp"
#include "io/binary/graph.hpp"

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

BOOST_FIXTURE_TEST_SUITE(robust_order_tests, fs::TmpFolderFixture)

BOOST_AUTO_TEST_CASE( OrderTest ) {
    string file_name = "src/test/debruijn/graph_fragments/saves/test_save";
    Graph graph(55);
    RandomGraphConstructor<Graph>(graph, /*max_size*/100).Generate(/*iterations*/1000);

    io::GraphIO<Graph> io;
    io.Save(file_name, graph);

    Graph new_graph(55);
    io.Load(file_name, new_graph);

    IteratorOrderChecker<Graph> checker(graph, new_graph);
    BOOST_CHECK(checker.CheckOrder(graph.SmartVertexBegin(), new_graph.SmartVertexBegin()));
    BOOST_CHECK(checker.CheckOrder(graph.SmartEdgeBegin(), new_graph.SmartEdgeBegin()));
//    BOOST_CHECK(checker.CheckOrder(graph.SmartVertexBegin(), new_graph.SmartVertexBegin()));
//    BOOST_CHECK(checker.CheckOrder(graph.SmartEdgeBegin(), new_graph.SmartEdgeBegin()));
}
BOOST_AUTO_TEST_SUITE_END()
}
