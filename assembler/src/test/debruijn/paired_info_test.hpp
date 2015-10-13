//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <boost/test/unit_test.hpp>
#include "omni/omni_utils.hpp"
#include "de/paired_info_helpers.hpp"

namespace debruijn_graph {

using namespace omnigraph::de;

class MockGraph {

public:
    typedef int EdgeId;

    MockGraph() {
        std::pair<EdgeId, EdgeId> data[] = {{1,  2},
                                            {3,  4},
                                            {8,  9},
                                            {13, 14}};
        for (const auto &pair : data) {
            conjs_[pair.first] = pair.second;
            conjs_[pair.second] = pair.first;
            lengths_[pair.first] = 1;
            lengths_[pair.second] = 1;
        }
    }

    EdgeId conjugate(EdgeId id) const {
        return conjs_.find(id)->second;
    }

    int length(EdgeId id) const {
        return lengths_.find(id)->second;
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
using ClMockIndex = PairedInfoIndexT<MockGraph>;
using EdgeSet = std::set<MockIndex::EdgeId>;

template<typename Index>
EdgeSet GetNeighbours(const Index &pi, MockGraph::EdgeId e) {
    EdgeSet result;
    for (auto i : pi.Get(e))
        result.insert(i.first);
    return result;
}

EdgeSet GetRawNeighbours(const MockIndex &pi, MockGraph::EdgeId e) {
    EdgeSet result;
    for (auto i : pi.RawGet(e))
        result.insert(i.first);
    return result;
}

using EdgeDataSet = std::set<omnigraph::de::PairInfo<MockIndex::EdgeId>>;

EdgeDataSet GetEdgePairInfo(const MockIndex &pi) {
    EdgeDataSet result;
    for (auto i = omnigraph::de::pair_begin(pi); i != omnigraph::de::pair_end(pi); ++i)
        for (auto j : *i)
            result.emplace(i.first(), i.second(), j);
    return result;
}

bool Contains(const MockIndex &pi, MockGraph::EdgeId e1, MockGraph::EdgeId e2, float distance) {
    for (auto p : pi.Get(e1, e2))
        if (math::eq(p.d, distance))
            return true;
    return false;
}

/*std::ostream& operator<<(std::ostream& str, const EdgeSet &set) {
    str << "{";
    for (const auto i : set)
        str << i << ",";
    str << "}";
    return str;
}*/

void PrintPi(const MockIndex &pi) {
    std::cout << "--Size: " << pi.size() << "--\n";
    for (auto i = pi.data_begin(); i != pi.data_end(); ++i ) {
        std::cout << i->first << ": " << i->second.size() << ": ";
        for (auto j : i->second) {
            std::cout << "[" << j.first << ": " << j.second.size() << ": ";
            for (auto k : j.second)
                std::cout << k << ", ";
            std::cout << "]";
        }
        std::cout << std::endl;
    }
    std::cout << "------------\n";
}

BOOST_AUTO_TEST_SUITE(pair_info_tests)

BOOST_AUTO_TEST_CASE(PairedInfoConstruct) {
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 8, RawPoint(1, 3));
    //Check for normal info
    BOOST_CHECK(Contains(pi, 1, 8, 1));
    //Check for reversed info
    BOOST_CHECK(Contains(pi, 8, 1, -1));
    //Check for conjugate info
    BOOST_CHECK(Contains(pi, 9, 2, 1));
    pi.Add(1, 3, {2, 2});
    pi.Add(1, 3, {3, 1});
    BOOST_CHECK(Contains(pi, 1, 8, 1));
    BOOST_CHECK(Contains(pi, 1, 3, 2));
    RawHistogram test1;
    test1.insert({2, 1});
    test1.insert({3, 2});
    BOOST_CHECK_EQUAL(pi.Get(1, 3).Unwrap(), test1);
    RawHistogram test2;
    test2.insert({-2, 1});
    test2.insert({-3, 2});
    BOOST_CHECK_EQUAL(pi.Get(3, 1).Unwrap(), test2);
}

BOOST_AUTO_TEST_CASE(PairedInfoRemove) {
    MockGraph graph;
    ClMockIndex pi(graph); //Deleting currently works only with the clustered index
    pi.Add(1, 2, {1, 1, 0});
    pi.Add(1, 3, {2, 2, 0});
    pi.Add(1, 3, {3, 1, 0});
    //Check for single pair remove
    pi.Remove(1, 3, {2, 1, 0});
    HistogramWithWeight test1;
    test1.insert({3, 2, 0});
    BOOST_CHECK_EQUAL(pi.Get(1, 3).Unwrap(), test1);
    pi.Remove(1, 3, {3, 1, 0});
    BOOST_CHECK(!pi.contains(1, 3));
    //Check for full remove
    pi.Add(1, 3, {2, 2, 0});
    pi.Add(1, 3, {3, 1, 0});
    pi.Remove(1, 3);
    BOOST_CHECK(!pi.contains(1, 3));
    //Check for neighbours remove
    pi.Add(1, 3, {3, 1, 0});
    pi.Add(2, 13, {7, 1, 0});
    pi.Add(8, 14, {3, 1, 0});
    pi.Add(13, 4, {5, 1, 0});
    BOOST_CHECK(pi.contains(1, 14));
    pi.Remove(1);
    EdgeSet test3 = {3, 8};
    BOOST_CHECK_EQUAL(GetNeighbours(pi, 14), test3);
}

BOOST_AUTO_TEST_CASE(PairedInfoUnexistingEdges) {
    //Check that accessing missing edges doesn't broke index
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 2, {1, 1});
    pi.Add(1, 3, {2, 2});
    EdgeSet empty;
    BOOST_CHECK(pi.Get(14).empty());
    BOOST_CHECK_EQUAL(GetNeighbours(pi, 5), empty);
    BOOST_CHECK_EQUAL(GetNeighbours(pi, 8), empty);
    RawHistogram empty_hist;
    BOOST_CHECK_EQUAL(pi.Get(1, 5).Unwrap(), empty_hist);
    BOOST_CHECK_EQUAL(pi.Get(1, 8).Unwrap(), empty_hist);
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

    EdgeSet test1 = {3, 4, 8, 9, 14};
    BOOST_CHECK_EQUAL(GetNeighbours(pi, 1), test1);

    EdgeSet test2 = {1, 3, 8};
    BOOST_CHECK_EQUAL(GetNeighbours(pi, 14), test2);
}

/*BOOST_AUTO_TEST_CASE(PairedInfoRawNeighbours) {
    MockGraph graph;
    MockIndex pi(graph);
    pi.Add(1, 3, RawPoint(1, 1));
    pi.Add(1, 9, RawPoint(2, 1));
    pi.Add(8, 14, RawPoint(3, 1));
    pi.Add(3, 2, RawPoint(4, 1));
    pi.Add(13, 4, RawPoint(5, 1));
    pi.Add(9, 2, RawPoint(6, 1));
    pi.Add(2, 13, RawPoint(7, 1));

    PrintPi(pi);

    EdgeSet test1 = {3, 4, 8, 9};
    BOOST_CHECK_EQUAL(GetRawNeighbours(pi, 1), test1);
    EdgeSet test2 = {};
    BOOST_CHECK_EQUAL(GetRawNeighbours(pi, 2), test2);
    EdgeSet test3 = {1};
    BOOST_CHECK_EQUAL(GetRawNeighbours(pi, 14), test3);
}*/

BOOST_AUTO_TEST_CASE(PairedInfoPairTraverse) {
    MockGraph graph;
    MockIndex pi(graph);
    //Check that an empty index has empty iterator range
    BOOST_CHECK(omnigraph::de::pair_begin(pi) == omnigraph::de::pair_end(pi));
    RawPoint p1 = {1, 1}, p2 = {2, 1};
    pi.Add(1, 3, p1);
    pi.Add(1, 9, p2);

    EdgeDataSet test1 = {{1, 3, p1}, {3, 1, -p1}, {1, 9, p2}, {9, 1, -p2}, {2, 4, -p1}, {4, 2, p1}, {8, 2, p2}, {2, 8, -p2}};
    BOOST_CHECK_EQUAL(GetEdgePairInfo(pi), test1);
}

BOOST_AUTO_TEST_SUITE_END()

}