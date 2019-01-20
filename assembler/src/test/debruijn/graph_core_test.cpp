//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"

#include <vector>
#include <set>
#include <string>

#include <gtest/gtest.h>

using namespace debruijn_graph;

TEST( GraphCore, Empty ) {
    Graph g(11);
    ASSERT_EQ(11u, g.k());
    ASSERT_EQ(0u, g.size());
}

TEST( GraphCore, OneVertex ) {
    Graph g(11);
    g.AddVertex();
    ASSERT_EQ(2u, g.size());
    VertexId v = *(g.begin());
    VertexId rcv = g.conjugate(v);
    ASSERT_NE(v, rcv);
    ASSERT_EQ(v, g.conjugate(rcv));
}

std::pair<std::vector<VertexId>, std::vector<EdgeId>> createGraph(Graph &graph, int edgeNumber) {
    std::vector<VertexId> v;
    std::vector<EdgeId> e;
    v.push_back(graph.AddVertex());
    for (int i = 0; i < edgeNumber; i++) {
        v.push_back(graph.AddVertex());
        e.push_back(
                graph.AddEdge(v[v.size() - 2], v[v.size() - 1],
                        Sequence("AAAAAAAAAAAAAAAAA")));
    }
    return make_pair(v, e);
}

TEST( GraphCore, OneEdge ) {
    Graph g(11);
    auto data = createGraph(g, 1);
    ASSERT_EQ(1u, g.OutgoingEdgeCount(data.first[0]));
    ASSERT_EQ(0u, g.OutgoingEdgeCount(data.first[1]));
    ASSERT_EQ(data.second[0], g.GetUniqueOutgoingEdge(data.first[0]));
    ASSERT_EQ(g.conjugate(data.second[0]),
            g.GetUniqueOutgoingEdge(g.conjugate(data.first[1])));
    ASSERT_EQ(data.second[0],
            g.conjugate(g.conjugate(data.second[0])));
    ASSERT_EQ(!(g.EdgeNucls(data.second[0])),
            g.EdgeNucls(g.conjugate(data.second[0])));
}

TEST( GraphCore, VertexMethods ) {
    Graph g(11);
    auto data = createGraph(g, 2);
    ASSERT_EQ(data.second[0], g.GetUniqueIncomingEdge(data.first[1]));
    ASSERT_EQ(data.second[0], g.GetUniqueOutgoingEdge(data.first[0]));
    ASSERT_EQ(false, g.CanCompressVertex(data.first[0]));
    ASSERT_EQ(true, g.CanCompressVertex(data.first[1]));
    ASSERT_EQ(false, g.CheckUniqueIncomingEdge(data.first[0]));
    ASSERT_EQ(true, g.CheckUniqueIncomingEdge(data.first[1]));
    ASSERT_EQ(false, g.CheckUniqueOutgoingEdge(data.first[2]));
    ASSERT_EQ(true, g.CheckUniqueOutgoingEdge(data.first[1]));
    ASSERT_EQ(true, g.IsDeadEnd(data.first[2]));
    ASSERT_EQ(false, g.IsDeadEnd(data.first[1]));
    ASSERT_EQ(true, g.IsDeadStart(data.first[0]));
    ASSERT_EQ(false, g.IsDeadStart(data.first[1]));
}

TEST( GraphCore, SmartIterator ) {
    Graph g(11);
    auto data = createGraph(g, 4);
    size_t num = 0;
    std::set<VertexId> visited;
//    std::less<VertexId> comp;
//    auto it = g.SmartVertexBegin(comp);
    for (auto it = g.SmartVertexBegin(); !it.IsEnd(); ++it) {
        num++;
        DEBUG( "with seq in vert" << g.VertexNucls(*it).str());
        visited.insert(*it);
    }
    ASSERT_EQ(num, data.first.size() * 2);
    for (size_t i = 0; i < data.first.size(); i++) {
        ASSERT_NE(visited.find(data.first[i]), visited.end());
        ASSERT_NE(visited.find(g.conjugate(data.first[i])), visited.end());
    }
}


TEST( GraphCore, SelfRCEdgeMerge ) {
    Graph g(5);
    VertexId v1 = g.AddVertex();
    VertexId v2 = g.AddVertex();
    EdgeId edge1 = g.AddEdge(v1, v2, Sequence("AACGCTATT"));
    EdgeId edge2 = g.AddEdge(v2, g.conjugate(v2), Sequence("CTATTCACGTGAATAG"));
    std::vector<EdgeId> path = {edge1, edge2};
    g.MergePath(path);
    ASSERT_EQ(2u, g.size());
    ASSERT_EQ(1u, g.OutgoingEdgeCount(v1));
    ASSERT_EQ(Sequence("AACGCTATTCACGTGAATAGCGTT"), g.EdgeNucls(g.GetUniqueOutgoingEdge(v1)));
}
