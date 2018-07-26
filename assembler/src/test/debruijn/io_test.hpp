//***************************************************************************
//* Copyright (c) 2015-2018 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include <boost/test/unit_test.hpp>
#include "test_utils.hpp"
#include "random_graph.hpp"
#include "assembly_graph/handlers/id_track_handler.hpp"
#include "io/binary/graphio.hpp"

namespace debruijn_graph {

template<typename T>
void CompareContainers(const T &lhs, const T &rhs) {
    auto i = lhs.begin();
    auto j = rhs.begin();
    for (; i != lhs.end() && j != rhs.end(); ++i, ++j) {
        BOOST_CHECK_EQUAL(*i, *j);
    }
    BOOST_CHECK(i == lhs.end());
    BOOST_CHECK(j == lhs.end());
}

template<typename I>
void CompareGraphIterators(I it1, I it2) {
    while(!it1.IsEnd() && !it2.IsEnd()) {
        BOOST_CHECK_EQUAL((*it1).int_id(), (*it2).int_id());
        ++it1;
        ++it2;
    }
    BOOST_CHECK(it1.IsEnd());
    BOOST_CHECK(it2.IsEnd());
}

const Graph &CommonGraph() {
    static Graph graph(55);
    if (!graph.size())
        RandomGraph<Graph>(graph, /*max_size*/100).Generate(/*iterations*/1000);
    return graph;
}

BOOST_FIXTURE_TEST_SUITE(robust_order_tests, fs::TmpFolderFixture)

using namespace io::binary;
const char *file_name = "src/test/debruijn/graph_fragments/saves/test_save";

BOOST_AUTO_TEST_CASE( OrderTest ) {
    const auto &graph = CommonGraph();

    GraphIO<Graph> io;
    io.Save(file_name, graph);

    Graph new_graph(graph.k());
    io.Load(file_name, new_graph);

    CompareGraphIterators(graph.SmartVertexBegin(), new_graph.SmartVertexBegin());
    CompareGraphIterators(graph.SmartEdgeBegin(), new_graph.SmartEdgeBegin());
}

BOOST_AUTO_TEST_CASE(TestPairedInfoIO) {
    using Index = UnclusteredPairedInfoIndexT<Graph>;
    const auto &graph = CommonGraph();

    Index pi(graph);
    RandomPairedIndex<Index>(pi, 100).Generate(100);

    PairedIndexIO<Index> io;
    io.Save(file_name, pi);

    io::IdMapper<Graph::EdgeId> mapper;
    for (auto i = graph.SmartEdgeBegin(); !i.IsEnd(); ++i)
        mapper[(*i).int_id()] = *i;
    Index ni(graph);
    io.Load(file_name, ni, mapper);

    BOOST_CHECK_EQUAL(pi.size(), ni.size());
    for (auto pit = omnigraph::de::pair_begin(pi), nit = omnigraph::de::pair_begin(ni);
         pit != omnigraph::de::pair_end(pi); ++pit, ++nit) {
        BOOST_CHECK_EQUAL(pit.first(), nit.first());
        BOOST_CHECK_EQUAL(pit.second(), nit.second());
        BOOST_CHECK_EQUAL(pit->size(), nit->size());

        for (auto ppit = pit->begin(), npit = nit->begin();
             ppit != pit->end(); ++ppit, ++npit) {
            BOOST_CHECK_EQUAL(ppit->weight, npit->weight);
            BOOST_CHECK_EQUAL(ppit->d, npit->d);
        }
    }
}

BOOST_AUTO_TEST_CASE(TestKmerMapperIO) {
    const auto &graph = CommonGraph();

    KmerMapper<Graph> kmer_mapper(graph);
    RandomKmerMapper<Graph>(kmer_mapper).Generate(100);

    KmerMapperIO<Graph> io;
    io.Save(file_name, kmer_mapper);

    KmerMapper<Graph> new_mapper(graph);
    io.Load(file_name, new_mapper);

    CompareContainers(kmer_mapper, new_mapper);
}

BOOST_AUTO_TEST_SUITE_END()
}
