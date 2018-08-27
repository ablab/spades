//***************************************************************************
//* Copyright (c) 2015-2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <boost/test/unit_test.hpp>
#include "paired_info/paired_info_helpers.hpp"
#include "random_graph.hpp"
#include "io/binary/paired_index.hpp"

namespace omnigraph {

namespace de {

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

BOOST_AUTO_TEST_SUITE(pair_info_tests)

BOOST_AUTO_TEST_CASE(PairedInfoConstruct) {
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 8, RawPoint(1, 3));
    BOOST_CHECK_EQUAL(pi.size(), 2);
    //Check for normal info
    BOOST_CHECK(pi.contains(1));
    BOOST_CHECK(pi.contains(1, 8));
    BOOST_CHECK(Contains(pi, 1, 8, 1));
    //Check for conjugate info
    BOOST_CHECK(pi.contains(9));
    BOOST_CHECK(pi.contains(9, 2));
    BOOST_CHECK(Contains(pi, 9, 2, 8));
    //Check for more points
    pi.Add(1, 3, {2, 2});
    pi.Add(1, 3, {3, 1});
    BOOST_CHECK(Contains(pi, 1, 8, 1));
    BOOST_CHECK(Contains(pi, 1, 3, 2));
    BOOST_CHECK(Contains(pi, 4, 2, 4));
    BOOST_CHECK_EQUAL(pi.size(), 6);
    //Check for loops
    pi.Add(1, 1, {0, 1});
    BOOST_CHECK_EQUAL(pi.size(), 8);
    BOOST_CHECK(Contains(pi, 1, 1, 0));
    BOOST_CHECK(Contains(pi, 2, 2, 0));
    //Check for self-conjugates
    //pi.Add(1, 2, {1, 1});
    //BOOST_CHECK_EQUAL(pi.size(), 9);
    //BOOST_CHECK(Contains(pi, 1, 2, 1));
}

BOOST_AUTO_TEST_CASE(PairedInfoAccess) {
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
    BOOST_CHECK_EQUAL(proxy1[1].Unwrap(), test0);
    BOOST_CHECK_EQUAL(proxy1[3].Unwrap(), test1);
    auto proxy3 = pi.Get(3);
    BOOST_CHECK_EQUAL(proxy3[7].Unwrap(), test0);
    BOOST_CHECK_EQUAL(proxy3[1].Unwrap(), test0);
    auto proxy2 = pi.Get(2);
    BOOST_CHECK_EQUAL(proxy2[1].Unwrap(), test0);
    BOOST_CHECK_EQUAL(proxy2[4].Unwrap(), test0);
    RawHistogram test4;
    test4.insert({4, 1});
    test4.insert({5, 2});
    auto proxy4 = pi.Get(4);
    BOOST_CHECK_EQUAL(proxy4[1].Unwrap(), test0);
    BOOST_CHECK_EQUAL(proxy4[2].Unwrap(), test4);
}

BOOST_AUTO_TEST_CASE(PairedInfoHalfAccess) {
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
    BOOST_CHECK_EQUAL(proxy1[1].Unwrap(), test0);
    BOOST_CHECK_EQUAL(proxy1[3].Unwrap(), test1);
    auto proxy3 = pi.GetHalf(3);
    BOOST_CHECK_EQUAL(proxy3[7].Unwrap(), test0);
    BOOST_CHECK_EQUAL(proxy3[1].Unwrap(), test0);
    auto proxy2 = pi.GetHalf(2);
    BOOST_CHECK_EQUAL(proxy2[1].Unwrap(), test0);
    BOOST_CHECK_EQUAL(proxy2[4].Unwrap(), test0);
    RawHistogram test4;
    test4.insert({4, 2});
    test4.insert({5, 1});
    auto proxy4 = pi.GetHalf(4);
    BOOST_CHECK_EQUAL(proxy4[1].Unwrap(), test0);
    BOOST_CHECK_EQUAL(proxy4[2].Unwrap(), test0);
}

//Backwards info is currently unused
/*BOOST_AUTO_TEST_CASE(PairedInfoBackData) {
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 8, {1, 3});
    pi.Add(1, 3, {2, 2});
    pi.Add(1, 3, {3, 1});
    RawHistogram test0;
    BOOST_CHECK_EQUAL(pi.GetBack(1, 3).Unwrap(), test0);
    RawHistogram test1;
    test1.insert({-2, 1});
    test1.insert({-3, 2});
    BOOST_CHECK_EQUAL(pi.GetBack(3, 1).Unwrap(), test1);
    RawHistogram test2;
    test2.insert({-4, 2});
    test2.insert({-5, 1});
    BOOST_CHECK_EQUAL(pi.GetBack(2, 4).Unwrap(), test2);
}*/

BOOST_AUTO_TEST_CASE(PairedInfoRemove) {
    MockGraph graph;
    MockClIndex pi(graph); //Deleting currently works only with the clustered index
    pi.Add(1, 8, {1, 1, 0});
    pi.Add(1, 3, {2, 2, 0});
    pi.Add(1, 3, {3, 1, 0});
    BOOST_CHECK_EQUAL(pi.size(), 6);
    //Check for single pair remove
    BOOST_CHECK_EQUAL(pi.Remove(1, 3, {2, 1, 0}), 2);
    BOOST_CHECK_EQUAL(pi.size(), 4);
    HistogramWithWeight test1;
    test1.insert({3, 2, 0});
    BOOST_CHECK_EQUAL(pi.Get(1, 3).Unwrap(), test1);
    BOOST_CHECK_EQUAL(pi.Remove(1, 3, {3, 1, 0}), 2);
    BOOST_CHECK(!pi.contains(1, 3));
    BOOST_CHECK_EQUAL(pi.size(), 2);
    //Check for full remove
    pi.Add(1, 3, {2, 2, 0});
    pi.Add(1, 3, {3, 1, 0});
    pi.Remove(1, 3);
    BOOST_CHECK(!pi.contains(1, 3));
    BOOST_CHECK(!pi.contains(4, 2));
    BOOST_CHECK_EQUAL(pi.size(), 2);
    //Check for nonexisting remove
    BOOST_CHECK_EQUAL(pi.Remove(1, 3, {1, 1, 0}), 0);
    BOOST_CHECK(!pi.contains(1, 3));
    BOOST_CHECK_EQUAL(pi.Remove(1, 2, {2, 2, 0}), 0);
    BOOST_CHECK(pi.contains(1, 8));
    BOOST_CHECK_EQUAL(pi.size(), 2);
    //Check for neighbours remove
    pi.Add(1, 3, {3, 1, 0});
    pi.Add(13, 2, {7, 1, 0});
    pi.Remove(1);
    BOOST_CHECK(!pi.contains(1, 3));
    BOOST_CHECK(!pi.contains(4, 2));
    BOOST_CHECK(!pi.contains(1, 14));
    BOOST_CHECK(!pi.contains(13, 2));
}

BOOST_AUTO_TEST_CASE(PairedInfoUnexistingEdges) {
    //Check that accessing missing edges doesn't broke index
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 8, {1, 1});
    pi.Add(1, 3, {2, 2});
    EdgeSet empty;
    BOOST_CHECK(pi.Get(14).empty());
    BOOST_CHECK_EQUAL(GetNeighbours(pi, 5), empty);
    BOOST_CHECK_EQUAL(GetNeighbours(pi, 8), empty);
    RawHistogram empty_hist;
    BOOST_CHECK_EQUAL(pi.Get(1, 5).Unwrap(), empty_hist);
    BOOST_CHECK_EQUAL(pi.Get(8, 1).Unwrap(), empty_hist);
}

BOOST_AUTO_TEST_CASE(PairedInfoPrune) {
    MockGraph graph;
    MockClIndex pi(graph); //Deleting currently works only with the clustered index
    pi.Add(1, 8, {1, 1, 0});
    pi.Add(1, 3, {2, 2, 0});
    //Check for auto-prune
    BOOST_CHECK_EQUAL(pi.Remove(1, 8, {1, 1, 0}), 2);
    BOOST_CHECK(!pi.contains(1, 2));
    BOOST_CHECK_EQUAL(pi.Remove(1, 3, {2, 1, 0}), 2);
    BOOST_CHECK(!pi.contains(1, 3));
    BOOST_CHECK(!pi.contains(1));
    //Check for nonexisting remove
    BOOST_CHECK_EQUAL(pi.Remove(1, 2, {1, 1, 0}), 0);
    BOOST_CHECK_EQUAL(pi.Remove(1, 3, {1, 1, 0}), 0);
    BOOST_CHECK_EQUAL(pi.Remove(1), 0);
}

BOOST_AUTO_TEST_CASE(PairedInfoNeighbours) {
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
    BOOST_CHECK_EQUAL(GetNeighbours(pi, 1), test1);

    EdgeSet test2 = {13};
    BOOST_CHECK_EQUAL(GetNeighbours(pi, 2), test2);

    EdgeSet test3 = {2, 14};
    BOOST_CHECK_EQUAL(GetNeighbours(pi, 3), test3);

    EdgeSet test14 = {1};
    BOOST_CHECK_EQUAL(GetNeighbours(pi, 14), test14);

    EdgeDataSet testF1 = {{1, 3, RawPoint(1, 1)}, {1, 4, RawPoint(2, 1)},
                          {1, 8, RawPoint(-1, 1)}, {1, 9, RawPoint(2, 1)}};
    BOOST_CHECK_EQUAL(GetNeighbourInfo(pi, 1), testF1);

    EdgeDataSet testF3 = {{3, 2, RawPoint(4, 1)}, {3, 14, RawPoint(-5, 1)}};
    BOOST_CHECK_EQUAL(GetNeighbourInfo(pi, 3), testF3);
}

BOOST_AUTO_TEST_CASE(PairedInfoDoubledInfo) {
    MockGraph graph;
    MockIndex pi(graph);

    pi.Add(1, 3, RawPoint(1, 1));
    pi.Add(1, 8, RawPoint(1, 1));
    pi.Add(4, 2, RawPoint(1, 1));

    EdgeSet test0;
    EdgeSet neighbours;
    //Check that neighbours don't repeat
    for (auto ep : pi.Get(1)) {
        BOOST_CHECK(neighbours.insert(ep.first).second);
        BOOST_CHECK(!ep.second.empty());
    }
    BOOST_CHECK_EQUAL(neighbours, EdgeSet({3, 8}));
    BOOST_CHECK_EQUAL(GetNeighbours(pi, 1), EdgeSet({3, 8}));
    //Check that the info is full
    EdgeDataSet testF1 = {{1, 3, RawPoint(1, 1)}, {1, 8, RawPoint(1, 1)}, {1, 3, RawPoint(-1, 1)}};
    BOOST_CHECK_EQUAL(GetNeighbourInfo(pi, 1), testF1);
}

/*BOOST_AUTO_TEST_CASE(PairedInfoRawData) {
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 3, RawPoint(2, 1));
    pi.Add(2, 4, RawPoint(-3, 1));
    RawHistogram test1;
    test1.insert({1, 1});
    test1.insert({2, 1});
    BOOST_CHECK_EQUAL(pi.Get(1, 3).Unwrap(), test1);
    RawHistogram test2;
    test2.insert({2, 1});
    BOOST_CHECK_EQUAL(pi.RawGet(1, 3).Unwrap(), test2);
    RawHistogram test3;
    test3.insert({-3, 1});
    BOOST_CHECK_EQUAL(pi.RawGet(2, 4).Unwrap(), test3);
    RawHistogram test2b;
    test2b.insert({-1, 1});
    BOOST_CHECK_EQUAL(pi.RawGet(3, 1).Unwrap(), test2b);
    RawHistogram test3b;
    test3b.insert({4, 1});
    BOOST_CHECK_EQUAL(pi.RawGet(4, 2).Unwrap(), test3b);
}*/

BOOST_AUTO_TEST_CASE(PairedInfoHalfNeighbours) {
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
    BOOST_CHECK_EQUAL(GetHalfNeighbours(pi, 1), test1);
    EdgeSet test2 = {13};
    BOOST_CHECK_EQUAL(GetHalfNeighbours(pi, 2), test2);
    EdgeSet test3 = {14};
    BOOST_CHECK_EQUAL(GetHalfNeighbours(pi, 3), test3);
    EdgeSet test4;
    BOOST_CHECK_EQUAL(GetHalfNeighbours(pi, 14), test4);
}

BOOST_AUTO_TEST_CASE(PairedInfoMoreNeighbours) {
    MockGraph graph;
    MockIndex pi(graph);
    //Check that an empty index has an empty iterator range
    BOOST_CHECK(pair_begin(pi) == pair_end(pi));
    EdgeDataSet empty;
    BOOST_CHECK_EQUAL(GetEdgePairInfo(pi), empty);
    RawPoint p0 = {0, 0}, p1 = {10, 1}, p2 = {20, 1};
    pi.Add(1, 3, p1);
    pi.Add(1, 9, p2);
    pi.Add(1, 1, p0);
    pi.Add(2, 4, p1);
    EdgeSet test1 = {1, 3, 9};
    BOOST_CHECK_EQUAL(GetNeighbours(pi, 1), test1);
    EdgeSet test2 = {2, 4};
    BOOST_CHECK_EQUAL(GetNeighbours(pi, 2), test2);
}


BOOST_AUTO_TEST_CASE(PairedInfoMerge) {
    MockGraph graph;
    MockIndex pi(graph), spi(graph), tpi(graph);
    RawPoint p = {1, 1};
    tpi.Add(1, 1, p);
    pi.Merge(tpi);
    BOOST_CHECK_EQUAL(pi.size(), 2);
    spi.Add(1, 9, p);
    spi.Add(3, 13, p);
    BOOST_CHECK_EQUAL(spi.size(), 4);
    pi.Merge(spi);
    BOOST_CHECK_EQUAL(pi.size(), 6);
    BOOST_CHECK(Contains(pi, 1, 9, 1));
    BOOST_CHECK(Contains(pi, 3, 13, 1));
}


BOOST_AUTO_TEST_CASE(PairedInfoPairTraverse) {
    MockGraph graph;
    MockIndex pi(graph);
    //Check that an empty index has an empty iterator range
    BOOST_CHECK(pair_begin(pi) == pair_end(pi));
    EdgeDataSet empty;
    BOOST_CHECK_EQUAL(GetEdgePairInfo(pi), empty);
    RawPoint p0 = {0, 0}, p1 = {10, 1}, p2 = {20, 1}, pj1 = {12, 1}, pj2 = {27, 1};
    pi.Add(1, 3, p1);
    pi.Add(1, 9, p2);
    pi.Add(1, 1, p0);
    pi.Add(2, 4, p1);
    BOOST_CHECK(omnigraph::de::pair_begin(pi) != omnigraph::de::pair_end(pi));
    //PrintPI(pi);
    EdgeDataSet test1 = {{1, 1, p0}, {2, 2, p0},
                         {1, 3, p1}, {4, 2, pj1},
                         {1, 9, p2}, {8, 2, pj2},
                         {2, 4, p1}, {3, 1, pj1}};
    BOOST_CHECK_EQUAL(GetEdgePairInfo(pi), test1);
}

using TestIndex = UnclusteredPairedInfoIndexT<Graph>;

BOOST_AUTO_TEST_CASE(PairedInfoRandomSymmetry) {
    Graph graph(55);
    debruijn_graph::RandomGraph<Graph>(graph, /*max_size*/100).Generate(/*iterations*/1000);

    TestIndex pi(graph);
    debruijn_graph::RandomPairedIndex<TestIndex>(pi, 100).Generate(20);

    for (auto it = pair_begin(pi); it != pair_end(pi); ++it) {
        auto info = *it;
        auto conj_info = pi.Get(graph.conjugate(it.second()), graph.conjugate(it.first()));
        BOOST_CHECK_EQUAL(info.size(), conj_info.size());
        auto offset = DEDistance(graph.length(it.first())) - DEDistance(graph.length(it.second()));
        for (auto i = info.begin(), ci = conj_info.begin(); i != info.end(); ++i, ++ci) {
            auto conj_point = *ci;
            conj_point.d += offset;
            BOOST_CHECK_EQUAL(*i, conj_point);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace de

} // namespace omnigraph
