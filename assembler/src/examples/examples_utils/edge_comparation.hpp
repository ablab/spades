//***************************************************************************
//* Copyright (c) 2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/assembly_graph/core/graph.hpp"

#include <string>
#include <unordered_set>

#include <gtest/gtest.h>

namespace spades_example {

void CompareEdges(const debruijn_graph::Graph& g1, const debruijn_graph::Graph& g2) {
    /*
    * Edges of de Bruijn graph present as reads of k-1 length in a form of std::string. Thus, edges of the graph
    * can be compared as a set of std::string.
    *
    * There are several iterators to go through edges in debruijn_graph::Graph. They are declared in
    * common/assembly_graph/core/observable_graph.hpp file.
    */

    INFO("Asserting edges");
    std::unordered_set <std::string> edges1;
    for (auto it = g1.ConstEdgeBegin(); !it.IsEnd(); ++it) {
        edges1.insert(g1.EdgeNucls(*it).str());
    }
    std::unordered_set <std::string> edges2;
    for (auto it = g2.ConstEdgeBegin(); !it.IsEnd(); ++it) {
        edges2.insert(g2.EdgeNucls(*it).str());
    }

    ASSERT_EQ(edges1.size(), edges2.size()) << "The number of edges in graphs are not equal";
    for (auto it = edges1.begin(); it != edges1.end(); ++it) {
        ASSERT_TRUE(edges2.count(*it) > 0) << "An edge contained in one graph is not contained in the second";
    }

    INFO("Edges of the graphs are equal");
}

}
