//***************************************************************************
//* Copyright (c) 2015-2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "random_graph.hpp"

#include "paired_info/index_point.hpp"
#include "paired_info/paired_info_helpers.hpp"
//#include "io/binary/paired_index.hpp"

#include <gtest/gtest.h>
#include <map>
#include <vector>

using namespace omnigraph::de;
using namespace debruijn_graph;

class MockGraph {

public:
    typedef int EdgeId;

    MockGraph() {
        std::pair<EdgeId, EdgeId> data[] = {{1,  2},
                                            {3,  4},
                                            {5,  7},
                                            {8,  9},
                                            {13, 14}};
        for (const auto &pair : data) {
            conjs_[pair.first] = pair.second;
            conjs_[pair.second] = pair.first;
            lengths_[pair.first] = pair.first;
            lengths_[pair.second] = pair.first;
        }
    }

    EdgeId conjugate(EdgeId id) const {
        return conjs_.find(id)->second;
    }

    int length(EdgeId id) const {
        return lengths_.find(id)->second;
    }

    int int_id(EdgeId id) const {
        return id;
    }

private:
    std::map<EdgeId, EdgeId> conjs_;
    std::map<EdgeId, int> lengths_;
};

struct EdgeData {
    MockGraph::EdgeId f, s;
    float d;
    EdgeData(MockGraph::EdgeId first, MockGraph::EdgeId second, float dist)
            : f(first), s(second), d(dist) {}
};

using MockIndex = UnclusteredPairedInfoIndexT<MockGraph>;
using MockClIndex = PairedInfoIndexT<MockGraph>;
using EdgeSet = std::set<MockIndex::EdgeId>;

template<typename Index>
EdgeSet GetNeighbours(const Index &pi, MockGraph::EdgeId e) {
    EdgeSet result;
    for (auto i : pi.Get(e))
        result.insert(i.first);
    return result;
}

EdgeSet GetHalfNeighbours(const MockIndex &pi, MockGraph::EdgeId e) {
    EdgeSet result;
    for (auto i : pi.GetHalf(e))
        result.insert(i.first);
    return result;
}

using EdgeDataSet = std::set<PairInfo<MockIndex::EdgeId>>;

EdgeDataSet GetEdgePairInfo(const MockIndex &pi) {
    EdgeDataSet result;
    for (auto i = pair_begin(pi); i != pair_end(pi); ++i) {
        for (auto j : *i) {
            result.emplace(i.first(), i.second(), j);
        }
    }
    return result;
}

EdgeDataSet GetNeighbourInfo(const MockIndex &pi, MockGraph::EdgeId e) {
    EdgeDataSet result;
    for (auto i : pi.Get(e))
        for (auto j : i.second)
            result.emplace(e, i.first, j);
    return result;
}

bool Contains(const MockIndex &pi, MockGraph::EdgeId e1, MockGraph::EdgeId e2, float distance) {
    for (auto p : pi.Get(e1, e2))
        if (math::eq(p.d, distance))
            return true;
    return false;
}

TEST(PairedInfo, Construct) {
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 8, RawPoint(1, 3));
    ASSERT_EQ(pi.size(), 2);
    //Check for normal info
    ASSERT_TRUE(pi.contains(1));
    ASSERT_TRUE(pi.contains(1, 8));
    ASSERT_TRUE(Contains(pi, 1, 8, 1));
    //Check for conjugate info
    ASSERT_TRUE(pi.contains(9));
    ASSERT_TRUE(pi.contains(9, 2));
    ASSERT_TRUE(Contains(pi, 9, 2, 8));
    //Check for more points
    pi.Add(1, 3, {2, 2});
    pi.Add(1, 3, {3, 1});
    ASSERT_TRUE(Contains(pi, 1, 8, 1));
    ASSERT_TRUE(Contains(pi, 1, 3, 2));
    ASSERT_TRUE(Contains(pi, 4, 2, 4));
    ASSERT_EQ(pi.size(), 6);
    //Check for loops
    pi.Add(1, 1, {0, 1});
    ASSERT_EQ(pi.size(), 8);
    ASSERT_TRUE(Contains(pi, 1, 1, 0));
    ASSERT_TRUE(Contains(pi, 2, 2, 0));
    //Check for self-conjugates
    //pi.Add(1, 2, {1, 1});
    //ASSERT_EQ(pi.size(), 9);
    //ASSERT_TRUE(Contains(pi, 1, 2, 1));
}

TEST(PairedInfo, Access) {
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 8, {1, 3});
    pi.Add(1, 3, {2, 2});
    pi.Add(1, 3, {3, 1});
    RawHistogram test0;
    RawHistogram test1;
    test1.insert({2, 1});
    test1.insert({3, 2});
    auto proxy1 = pi.Get(1);
    ASSERT_EQ(proxy1[1].Unwrap(), test0);
    ASSERT_EQ(proxy1[3].Unwrap(), test1);
    auto proxy3 = pi.Get(3);
    ASSERT_EQ(proxy3[7].Unwrap(), test0);
    ASSERT_EQ(proxy3[1].Unwrap(), test0);
    auto proxy2 = pi.Get(2);
    ASSERT_EQ(proxy2[1].Unwrap(), test0);
    ASSERT_EQ(proxy2[4].Unwrap(), test0);
    RawHistogram test4;
    test4.insert({4, 1});
    test4.insert({5, 2});
    auto proxy4 = pi.Get(4);
    ASSERT_EQ(proxy4[1].Unwrap(), test0);
    ASSERT_EQ(proxy4[2].Unwrap(), test4);
}

TEST(PairedInfo, HalfAccess) {
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 8, {1, 3});
    pi.Add(1, 3, {2, 2});
    pi.Add(1, 3, {3, 1});
    RawHistogram test0;
    RawHistogram test1;
    test1.insert({2, 1});
    test1.insert({3, 2});
    auto proxy1 = pi.GetHalf(1);
    ASSERT_EQ(proxy1[1].Unwrap(), test0);
    ASSERT_EQ(proxy1[3].Unwrap(), test1);
    auto proxy3 = pi.GetHalf(3);
    ASSERT_EQ(proxy3[7].Unwrap(), test0);
    ASSERT_EQ(proxy3[1].Unwrap(), test0);
    auto proxy2 = pi.GetHalf(2);
    ASSERT_EQ(proxy2[1].Unwrap(), test0);
    ASSERT_EQ(proxy2[4].Unwrap(), test0);
    RawHistogram test4;
    test4.insert({4, 2});
    test4.insert({5, 1});
    auto proxy4 = pi.GetHalf(4);
    ASSERT_EQ(proxy4[1].Unwrap(), test0);
    ASSERT_EQ(proxy4[2].Unwrap(), test0);
}

//Backwards info is currently unused
/*TEST(PairedInfo, PairedInfoBackData) {
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 8, {1, 3});
    pi.Add(1, 3, {2, 2});
    pi.Add(1, 3, {3, 1});
    RawHistogram test0;
    ASSERT_EQ(pi.GetBack(1, 3).Unwrap(), test0);
    RawHistogram test1;
    test1.insert({-2, 1});
    test1.insert({-3, 2});
    ASSERT_EQ(pi.GetBack(3, 1).Unwrap(), test1);
    RawHistogram test2;
    test2.insert({-4, 2});
    test2.insert({-5, 1});
    ASSERT_EQ(pi.GetBack(2, 4).Unwrap(), test2);
}*/

TEST(PairedInfo, Remove) {
    MockGraph graph;
    MockClIndex pi(graph); //Deleting currently works only with the clustered index
    pi.Add(1, 8, {1, 1, 0});
    pi.Add(1, 3, {2, 2, 0});
    pi.Add(1, 3, {3, 1, 0});
    ASSERT_EQ(pi.size(), 6);
    //Check for single pair remove
    ASSERT_EQ(pi.Remove(1, 3, {2, 1, 0}), 2);
    ASSERT_EQ(pi.size(), 4);
    HistogramWithWeight test1;
    test1.insert({3, 2, 0});
    ASSERT_EQ(pi.Get(1, 3).Unwrap(), test1);
    ASSERT_EQ(pi.Remove(1, 3, {3, 1, 0}), 2);
    ASSERT_FALSE(pi.contains(1, 3));
    ASSERT_EQ(pi.size(), 2);
    //Check for full remove
    pi.Add(1, 3, {2, 2, 0});
    pi.Add(1, 3, {3, 1, 0});
    pi.Remove(1, 3);
    ASSERT_FALSE(pi.contains(1, 3));
    ASSERT_FALSE(pi.contains(4, 2));
    ASSERT_EQ(pi.size(), 2);
    //Check for nonexisting remove
    ASSERT_EQ(pi.Remove(1, 3, {1, 1, 0}), 0);
    ASSERT_FALSE(pi.contains(1, 3));
    ASSERT_EQ(pi.Remove(1, 2, {2, 2, 0}), 0);
    ASSERT_TRUE(pi.contains(1, 8));
    ASSERT_EQ(pi.size(), 2);
    //Check for neighbours remove
    pi.Add(1, 3, {3, 1, 0});
    pi.Add(13, 2, {7, 1, 0});
    pi.Remove(1);
    ASSERT_FALSE(pi.contains(1, 3));
    ASSERT_FALSE(pi.contains(4, 2));
    ASSERT_FALSE(pi.contains(1, 14));
    ASSERT_FALSE(pi.contains(13, 2));
}

TEST(PairedInfo, UnexistingEdges) {
    //Check that accessing missing edges doesn't broke index
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 8, {1, 1});
    pi.Add(1, 3, {2, 2});
    EdgeSet empty;
    ASSERT_TRUE(pi.Get(14).empty());
    ASSERT_EQ(GetNeighbours(pi, 5), empty);
    ASSERT_EQ(GetNeighbours(pi, 8), empty);
    RawHistogram empty_hist;
    ASSERT_EQ(pi.Get(1, 5).Unwrap(), empty_hist);
    ASSERT_EQ(pi.Get(8, 1).Unwrap(), empty_hist);
}

TEST(PairedInfo, Prune) {
    MockGraph graph;
    MockClIndex pi(graph); //Deleting currently works only with the clustered index
    pi.Add(1, 8, {1, 1, 0});
    pi.Add(1, 3, {2, 2, 0});
    //Check for auto-prune
    ASSERT_EQ(pi.Remove(1, 8, {1, 1, 0}), 2);
    ASSERT_FALSE(pi.contains(1, 2));
    ASSERT_EQ(pi.Remove(1, 3, {2, 1, 0}), 2);
    ASSERT_FALSE(pi.contains(1, 3));
    ASSERT_FALSE(pi.contains(1));
    //Check for nonexisting remove
    ASSERT_EQ(pi.Remove(1, 2, {1, 1, 0}), 0);
    ASSERT_EQ(pi.Remove(1, 3, {1, 1, 0}), 0);
    ASSERT_EQ(pi.Remove(1), 0);
}

TEST(PairedInfo, Neighbours) {
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 3, RawPoint(1, 1));
    pi.Add(1, 9, RawPoint(2, 1));
    pi.Add(8, 14, RawPoint(3, 1));
    pi.Add(3, 2, RawPoint(4, 1));
    pi.Add(13, 4, RawPoint(5, 1));
    pi.Add(9, 2, RawPoint(6, 1));
    pi.Add(2, 13, RawPoint(7, 1));

    EdgeSet test1 = {3, 4, 8, 9};
    ASSERT_EQ(GetNeighbours(pi, 1), test1);

    EdgeSet test2 = {13};
    ASSERT_EQ(GetNeighbours(pi, 2), test2);

    EdgeSet test3 = {2, 14};
    ASSERT_EQ(GetNeighbours(pi, 3), test3);

    EdgeSet test14 = {1};
    ASSERT_EQ(GetNeighbours(pi, 14), test14);

    EdgeDataSet testF1 = {{1, 3, RawPoint(1, 1)}, {1, 4, RawPoint(2, 1)},
                          {1, 8, RawPoint(-1, 1)}, {1, 9, RawPoint(2, 1)}};
    ASSERT_EQ(GetNeighbourInfo(pi, 1), testF1);

    EdgeDataSet testF3 = {{3, 2, RawPoint(4, 1)}, {3, 14, RawPoint(-5, 1)}};
    ASSERT_EQ(GetNeighbourInfo(pi, 3), testF3);
}

TEST(PairedInfo, DoubledInfo) {
    MockGraph graph;
    MockIndex pi(graph);

    pi.Add(1, 3, RawPoint(1, 1));
    pi.Add(1, 8, RawPoint(1, 1));
    pi.Add(4, 2, RawPoint(1, 1));

    EdgeSet test0;
    EdgeSet neighbours;
    //Check that neighbours don't repeat
    for (auto ep : pi.Get(1)) {
        ASSERT_TRUE(neighbours.insert(ep.first).second);
        ASSERT_FALSE(ep.second.empty());
    }
    ASSERT_EQ(neighbours, EdgeSet({3, 8}));
    ASSERT_EQ(GetNeighbours(pi, 1), EdgeSet({3, 8}));
    //Check that the info is full
    EdgeDataSet testF1 = {{1, 3, RawPoint(1, 1)}, {1, 8, RawPoint(1, 1)}, {1, 3, RawPoint(-1, 1)}};
    ASSERT_EQ(GetNeighbourInfo(pi, 1), testF1);
}

/*TEST(PairedInfo, PairedInfoRawData) {
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 3, RawPoint(2, 1));
    pi.Add(2, 4, RawPoint(-3, 1));
    RawHistogram test1;
    test1.insert({1, 1});
    test1.insert({2, 1});
    ASSERT_EQ(pi.Get(1, 3).Unwrap(), test1);
    RawHistogram test2;
    test2.insert({2, 1});
    ASSERT_EQ(pi.RawGet(1, 3).Unwrap(), test2);
    RawHistogram test3;
    test3.insert({-3, 1});
    ASSERT_EQ(pi.RawGet(2, 4).Unwrap(), test3);
    RawHistogram test2b;
    test2b.insert({-1, 1});
    ASSERT_EQ(pi.RawGet(3, 1).Unwrap(), test2b);
    RawHistogram test3b;
    test3b.insert({4, 1});
    ASSERT_EQ(pi.RawGet(4, 2).Unwrap(), test3b);
}*/

TEST(PairedInfo, HalfNeighbours) {
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 3, RawPoint(1, 1));
    pi.Add(1, 9, RawPoint(2, 1));
    pi.Add(8, 14, RawPoint(3, 1));
    pi.Add(3, 2, RawPoint(4, 1));
    pi.Add(13, 4, RawPoint(5, 1));
    pi.Add(9, 2, RawPoint(6, 1));
    pi.Add(2, 13, RawPoint(7, 1));

    //PrintPi(pi);

    EdgeSet test1 = {3, 4, 8, 9};
    ASSERT_EQ(GetHalfNeighbours(pi, 1), test1);
    EdgeSet test2 = {13};
    ASSERT_EQ(GetHalfNeighbours(pi, 2), test2);
    EdgeSet test3 = {14};
    ASSERT_EQ(GetHalfNeighbours(pi, 3), test3);
    EdgeSet test4;
    ASSERT_EQ(GetHalfNeighbours(pi, 14), test4);
}

TEST(PairedInfo, MoreNeighbours) {
    MockGraph graph;
    MockIndex pi(graph);
    //Check that an empty index has an empty iterator range
    ASSERT_TRUE(pair_begin(pi) == pair_end(pi));
    EdgeDataSet empty;
    ASSERT_EQ(GetEdgePairInfo(pi), empty);
    RawPoint p0 = {0, 0}, p1 = {10, 1}, p2 = {20, 1};
    pi.Add(1, 3, p1);
    pi.Add(1, 9, p2);
    pi.Add(1, 1, p0);
    pi.Add(2, 4, p1);
    EdgeSet test1 = {1, 3, 9};
    ASSERT_EQ(GetNeighbours(pi, 1), test1);
    EdgeSet test2 = {2, 4};
    ASSERT_EQ(GetNeighbours(pi, 2), test2);
}


TEST(PairedInfo, Merge) {
    MockGraph graph;
    MockIndex pi(graph), spi(graph), tpi(graph);
    RawPoint p = {1, 1};
    tpi.Add(1, 1, p);
    pi.Merge(tpi);
    ASSERT_EQ(pi.size(), 2);
    spi.Add(1, 9, p);
    spi.Add(3, 13, p);
    ASSERT_EQ(spi.size(), 4);
    pi.Merge(spi);
    ASSERT_EQ(pi.size(), 6);
    ASSERT_TRUE(Contains(pi, 1, 9, 1));
    ASSERT_TRUE(Contains(pi, 3, 13, 1));
}


TEST(PairedInfo, PairTraverse) {
    MockGraph graph;
    MockIndex pi(graph);
    //Check that an empty index has an empty iterator range
    ASSERT_TRUE(pair_begin(pi) == pair_end(pi));
    EdgeDataSet empty;
    ASSERT_EQ(GetEdgePairInfo(pi), empty);
    RawPoint p0 = {0, 0}, p1 = {10, 1}, p2 = {20, 1}, pj1 = {12, 1}, pj2 = {27, 1};
    pi.Add(1, 3, p1);
    pi.Add(1, 9, p2);
    pi.Add(1, 1, p0);
    pi.Add(2, 4, p1);
    ASSERT_NE(omnigraph::de::pair_begin(pi), omnigraph::de::pair_end(pi));
    //PrintPI(pi);
    EdgeDataSet test1 = {{1, 1, p0}, {2, 2, p0},
                         {1, 3, p1}, {4, 2, pj1},
                         {1, 9, p2}, {8, 2, pj2},
                         {2, 4, p1}, {3, 1, pj1}};
    ASSERT_EQ(GetEdgePairInfo(pi), test1);
}

using TestIndex = UnclusteredPairedInfoIndexT<debruijn_graph::Graph>;

TEST(PairedInfo, RandomSymmetry) {
    debruijn_graph::Graph graph(55);
    debruijn_graph::RandomGraph<debruijn_graph::Graph>(graph, /*max_size*/100).Generate(/*iterations*/1000);

    TestIndex pi(graph);
    debruijn_graph::RandomPairedIndex<TestIndex>(pi, 100).Generate(20);

    for (auto it = pair_begin(pi); it != pair_end(pi); ++it) {
        auto info = *it;
        auto conj_info = pi.Get(graph.conjugate(it.second()), graph.conjugate(it.first()));
        ASSERT_EQ(info.size(), conj_info.size());
        auto offset = DEDistance(graph.length(it.first())) - DEDistance(graph.length(it.second()));
        for (auto i = info.begin(), ci = conj_info.begin(); i != info.end(); ++i, ++ci) {
            auto conj_point = *ci;
            conj_point.d += offset;
            ASSERT_EQ(*i, conj_point);
        }
    }
}
