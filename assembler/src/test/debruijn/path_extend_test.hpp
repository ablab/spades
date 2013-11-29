/*
 * path_extend.hpp
 *
 *  Created on: Nov 29, 2013
 *      Author: andrey
 */


#pragma once

#include <boost/test/unit_test.hpp>

#include "test_utils.hpp"
#include "omni/omni_utils.hpp"
#include "path_extend/path_visualizer.hpp"

namespace path_extend {

BOOST_FIXTURE_TEST_SUITE(path_extend_basic, TmpFolderFixture)

BOOST_AUTO_TEST_CASE( BidirectionalPathConstructor ) {
    Graph g(13);
    IdTrackHandler<Graph> int_ids(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g, int_ids);
    EdgeId start = *g.ConstEdgeBegin();

    BidirectionalPath path1(g);
    BOOST_CHECK_EQUAL(path1.Size(), 0);
    BOOST_CHECK_EQUAL(path1.Length(), 0);
    BOOST_CHECK_EQUAL(path1.Empty(), true);

    BidirectionalPath path2(g, start);
    BOOST_CHECK_EQUAL(path2.Size(), 1);
    BOOST_CHECK_EQUAL(path2.Length(), 426);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 426);
    BOOST_CHECK_EQUAL(path2.GapAt(0), 0);
    BOOST_CHECK_EQUAL(path2[0], start);
    BOOST_CHECK_EQUAL(path2.At(0), start);

    EdgeId e1 = g.conjugate(start);
    EdgeId e2 = *(g.OutgoingEdges(g.EdgeEnd(e1)).begin());
    vector<EdgeId> v;
    v.push_back(e1);
    v.push_back(e2);
    BidirectionalPath path3(g, v);
    BOOST_CHECK_EQUAL(path3.Size(), 2);
    BOOST_CHECK_EQUAL(path3.Length(), 427);
    BOOST_CHECK_EQUAL(path3.LengthAt(0), 427);
    BOOST_CHECK_EQUAL(path3.LengthAt(1), 1);
    BOOST_CHECK_EQUAL(path3.GapAt(1), 0);
    BOOST_CHECK_EQUAL(path3[0], e1);
    BOOST_CHECK_EQUAL(path3.At(1), e2);

    BidirectionalPath path4(path3);
    BOOST_CHECK_EQUAL(path4.Size(), 2);
    BOOST_CHECK_EQUAL(path4.Length(), 427);
    BOOST_CHECK_EQUAL(path4.LengthAt(0), 427);
    BOOST_CHECK_EQUAL(path4.LengthAt(1), 1);
    BOOST_CHECK_EQUAL(path4.GapAt(1), 0);
    BOOST_CHECK_EQUAL(path4[0], e1);
    BOOST_CHECK_EQUAL(path4.At(1), e2);

    path2.Clear();
    BOOST_CHECK_EQUAL(path2.Empty(), true);

    path4.Clear();
    BOOST_CHECK_EQUAL(path4.Empty(), true);
}

BOOST_AUTO_TEST_CASE( BidirectionalPathAddRemove ) {
    Graph g(13);
    IdTrackHandler<Graph> int_ids(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g, int_ids);
    EdgeId start = *g.ConstEdgeBegin();

    BidirectionalPath path1(g);
    BOOST_CHECK_EQUAL(path1.Size(), 0);
    BOOST_CHECK_EQUAL(path1.Length(), 0);

    // 98 26 145 70
    EdgeId e1 = g.conjugate(start);
    EdgeId e2 = *(g.OutgoingEdges(g.EdgeEnd(e1)).begin());
    EdgeId e3 = *(g.OutgoingEdges(g.EdgeEnd(e2)).begin());
    EdgeId e4 = *(g.OutgoingEdges(g.EdgeEnd(e3)).begin());

    path1.PushBack(e1);
    BOOST_CHECK_EQUAL(path1.Size(), 1);
    BOOST_CHECK_EQUAL(path1.Length(), 426);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 426);

    path1.PushBack(e2);
    BOOST_CHECK_EQUAL(path1.Size(), 2);
    BOOST_CHECK_EQUAL(path1.Length(), 427);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 427);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 1);
    BOOST_CHECK_EQUAL(path1.GapAt(1), 0);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    path1.PushBack(e3);
    BOOST_CHECK_EQUAL(path1.Size(), 3);
    BOOST_CHECK_EQUAL(path1.Length(), 1006);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 1006);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 580);
    BOOST_CHECK_EQUAL(path1.LengthAt(2), 579);
    BOOST_CHECK_EQUAL(path1.GapAt(0), 0);
    BOOST_CHECK_EQUAL(path1.GapAt(2), 0);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1[2], e3);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    path1.PopBack();
    BOOST_CHECK_EQUAL(path1.Size(), 2);
    BOOST_CHECK_EQUAL(path1.Length(), 427);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 427);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 1);
    BOOST_CHECK_EQUAL(path1.GapAt(1), 0);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    BidirectionalPath path2(g);
    path2.PushBack(e3);
    BOOST_CHECK_EQUAL(path2.Size(), 1);
    BOOST_CHECK_EQUAL(path2.Length(), 579);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 579);

    path2.PushBack(e4);
    BOOST_CHECK_EQUAL(path2.Empty(), false);
    BOOST_CHECK_EQUAL(path2.Size(), 2);
    BOOST_CHECK_EQUAL(path2.Length(), 635);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 635);
    BOOST_CHECK_EQUAL(path2.LengthAt(1), 56);
    BOOST_CHECK_EQUAL(path2[0], e3);
    BOOST_CHECK_EQUAL(path2[1], e4);

    path1.PushBack(path2);
    BOOST_CHECK_EQUAL(path1.Size(), 4);
    BOOST_CHECK_EQUAL(path1.Length(), 1062);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 1062);
    BOOST_CHECK_EQUAL(path1.LengthAt(2), 635);
    BOOST_CHECK_EQUAL(path1.LengthAt(3), 56);
    BOOST_CHECK_EQUAL(path1.GapAt(0), 0);
    BOOST_CHECK_EQUAL(path1.GapAt(2), 0);
    BOOST_CHECK_EQUAL(path1.GapAt(3), 0);
    BOOST_CHECK_EQUAL(path1[1], e2);
    BOOST_CHECK_EQUAL(path1[2], e3);
    BOOST_CHECK_EQUAL(path1[3], e4);

    path1.PopBack(1);
    BOOST_CHECK_EQUAL(path1.Size(), 3);
    BOOST_CHECK_EQUAL(path1.Length(), 1006);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 1006);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 580);
    BOOST_CHECK_EQUAL(path1.LengthAt(2), 579);
    BOOST_CHECK_EQUAL(path1.GapAt(0), 0);
    BOOST_CHECK_EQUAL(path1.GapAt(2), 0);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1[2], e3);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    path1.PopBack(3);
    BOOST_CHECK_EQUAL(path1.Length(), 0);
    BOOST_CHECK_EQUAL(path1.Empty(), true);

    path1.PushBack(e1);
    BOOST_CHECK_EQUAL(path1.Size(), 1);
    BOOST_CHECK_EQUAL(path1.Length(), 426);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 426);

    path1.PushBack(e2);
    BOOST_CHECK_EQUAL(path1.Size(), 2);
    BOOST_CHECK_EQUAL(path1.Length(), 427);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 427);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 1);
    BOOST_CHECK_EQUAL(path1.GapAt(1), 0);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    path1.PopBack(10);
    BOOST_CHECK_EQUAL(path1.Length(), 0);
    BOOST_CHECK_EQUAL(path1.Empty(), true);

    path1.PopBack(10);
    BOOST_CHECK_EQUAL(path1.Length(), 0);
    BOOST_CHECK_EQUAL(path1.Empty(), true);

    path2.PopBack(1);
    BOOST_CHECK_EQUAL(path2.Size(), 1);
    BOOST_CHECK_EQUAL(path2.Length(), 579);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 579);

    path2.PopBack(2);
    BOOST_CHECK_EQUAL(path2.Length(), 0);
    BOOST_CHECK_EQUAL(path2.Empty(), true);

    path2.PopBack();
    BOOST_CHECK_EQUAL(path2.Length(), 0);
    BOOST_CHECK_EQUAL(path2.Empty(), true);

    path2.PushBack(e3, 0);
    BOOST_CHECK_EQUAL(path2.Empty(), false);
    BOOST_CHECK_EQUAL(path2.Size(), 1);
    BOOST_CHECK_EQUAL(path2.Length(), 579);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 579);

    path2.PushBack(e4, 0);
    BOOST_CHECK_EQUAL(path2.Empty(), false);
    BOOST_CHECK_EQUAL(path2.Size(), 2);
    BOOST_CHECK_EQUAL(path2.Length(), 635);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 635);
    BOOST_CHECK_EQUAL(path2.LengthAt(1), 56);
    BOOST_CHECK_EQUAL(path2[0], e3);
    BOOST_CHECK_EQUAL(path2[1], e4);

    path2.Clear();
    BOOST_CHECK_EQUAL(path2.Empty(), true);
}

BOOST_AUTO_TEST_CASE( BidirectionalPathAddRemoveGaps ) {
    Graph g(13);
    IdTrackHandler<Graph> int_ids(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g, int_ids);
    EdgeId start = *g.ConstEdgeBegin();

    BidirectionalPath path1(g);

    // 98 26 145 70
    EdgeId e1 = g.conjugate(start);
    EdgeId e2 = *(g.OutgoingEdges(g.EdgeEnd(e1)).begin());
    EdgeId e3 = *(g.OutgoingEdges(g.EdgeEnd(e2)).begin());
    EdgeId e4 = *(g.OutgoingEdges(g.EdgeEnd(e3)).begin());

    path1.PushBack(e1);
    BOOST_CHECK_EQUAL(path1.Size(), 1);
    BOOST_CHECK_EQUAL(path1.Length(), 426);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 426);

    path1.PushBack(e2, 10);
    BOOST_CHECK_EQUAL(path1.Size(), 2);
    BOOST_CHECK_EQUAL(path1.Length(), 437);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 437);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 1);
    BOOST_CHECK_EQUAL(path1.GapAt(1), 10);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    path1.PushBack(e3, 10);
    BOOST_CHECK_EQUAL(path1.Size(), 3);
    BOOST_CHECK_EQUAL(path1.Length(), 1026);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 1026);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 590);
    BOOST_CHECK_EQUAL(path1.LengthAt(2), 579);
    BOOST_CHECK_EQUAL(path1.GapAt(0), 0);
    BOOST_CHECK_EQUAL(path1.GapAt(2), 10);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1[2], e3);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    path1.PopBack();
    BOOST_CHECK_EQUAL(path1.Size(), 2);
    BOOST_CHECK_EQUAL(path1.Length(), 437);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 437);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 1);
    BOOST_CHECK_EQUAL(path1.GapAt(1), 10);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    BidirectionalPath path2(g);
    //FRONT GAP SHOULD BE IGNORED
    path2.PushBack(e3, 10);
    BOOST_CHECK_EQUAL(path2.Size(), 1);
    BOOST_CHECK_EQUAL(path2.Length(), 579);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 579);

    path2.PushBack(e4, 100);
    BOOST_CHECK_EQUAL(path2.Empty(), false);
    BOOST_CHECK_EQUAL(path2.Size(), 2);
    BOOST_CHECK_EQUAL(path2.Length(), 735);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 735);
    BOOST_CHECK_EQUAL(path2.LengthAt(1), 56);
    BOOST_CHECK_EQUAL(path2.GapAt(0), 0);
    BOOST_CHECK_EQUAL(path2.GapAt(1), 100);
    BOOST_CHECK_EQUAL(path2[0], e3);
    BOOST_CHECK_EQUAL(path2[1], e4);

    path1.PushBack(path2);
    BOOST_CHECK_EQUAL(path1.Size(), 4);
    BOOST_CHECK_EQUAL(path1.Length(), 1172);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 1172);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 736);
    BOOST_CHECK_EQUAL(path1.LengthAt(2), 735);
    BOOST_CHECK_EQUAL(path1.LengthAt(3), 56);
    BOOST_CHECK_EQUAL(path1.GapAt(0), 0);
    BOOST_CHECK_EQUAL(path1.GapAt(1), 10);
    BOOST_CHECK_EQUAL(path1.GapAt(2), 0);
    BOOST_CHECK_EQUAL(path1.GapAt(3), 100);
    BOOST_CHECK_EQUAL(path1[1], e2);
    BOOST_CHECK_EQUAL(path1[2], e3);
    BOOST_CHECK_EQUAL(path1[3], e4);

    path1.PopBack(1);
    BOOST_CHECK_EQUAL(path1.Size(), 3);
    BOOST_CHECK_EQUAL(path1.Length(), 1016);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 1016);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 580);
    BOOST_CHECK_EQUAL(path1.LengthAt(2), 579);
    BOOST_CHECK_EQUAL(path1.GapAt(0), 0);
    BOOST_CHECK_EQUAL(path1.GapAt(1), 10);
    BOOST_CHECK_EQUAL(path1.GapAt(2), 0);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1[2], e3);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    path1.PopBack(3);
    BOOST_CHECK_EQUAL(path1.Length(), 0);
    BOOST_CHECK_EQUAL(path1.Empty(), true);

    path1.PushBack(e1, 1000);
    BOOST_CHECK_EQUAL(path1.Size(), 1);
    BOOST_CHECK_EQUAL(path1.Length(), 426);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 426);

    path1.PushBack(e2, -10);
    BOOST_CHECK_EQUAL(path1.Size(), 2);
    BOOST_CHECK_EQUAL(path1.Length(), 417);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 417);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 1);
    BOOST_CHECK_EQUAL(path1.GapAt(0), 0);
    BOOST_CHECK_EQUAL(path1.GapAt(1), -10);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    path1.PopBack(10);
    BOOST_CHECK_EQUAL(path1.Length(), 0);
    BOOST_CHECK_EQUAL(path1.Empty(), true);

    path1.PopBack(10);
    BOOST_CHECK_EQUAL(path1.Length(), 0);
    BOOST_CHECK_EQUAL(path1.Empty(), true);

    path2.PopBack(1);
    BOOST_CHECK_EQUAL(path2.Size(), 1);
    BOOST_CHECK_EQUAL(path2.Length(), 579);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 579);

    path2.PopBack(2);
    BOOST_CHECK_EQUAL(path2.Length(), 0);
    BOOST_CHECK_EQUAL(path2.Empty(), true);

    path2.PopBack();
    BOOST_CHECK_EQUAL(path2.Length(), 0);
    BOOST_CHECK_EQUAL(path2.Empty(), true);

    path2.PushBack(e3);
    BOOST_CHECK_EQUAL(path2.Empty(), false);
    BOOST_CHECK_EQUAL(path2.Size(), 1);
    BOOST_CHECK_EQUAL(path2.Length(), 579);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 579);

    path2.PushBack(e4, -1);
    BOOST_CHECK_EQUAL(path2.Empty(), false);
    BOOST_CHECK_EQUAL(path2.Size(), 2);
    BOOST_CHECK_EQUAL(path2.Length(), 634);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 634);
    BOOST_CHECK_EQUAL(path2.LengthAt(1), 56);
    BOOST_CHECK_EQUAL(path2.GapAt(0), 0);
    BOOST_CHECK_EQUAL(path2.GapAt(1), -1);
    BOOST_CHECK_EQUAL(path2[0], e3);
    BOOST_CHECK_EQUAL(path2[1], e4);

    path1.Clear();
    BOOST_CHECK_EQUAL(path1.Empty(), true);
    path1.PushBack(e1, 0);
    path1.PushBack(e2, 9);
    path1.PushBack(path2);
    BOOST_CHECK_EQUAL(path1.Size(), 4);
    BOOST_CHECK_EQUAL(path1.Length(), 1070);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 1070);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 635);
    BOOST_CHECK_EQUAL(path1.LengthAt(2), 634);
    BOOST_CHECK_EQUAL(path1.LengthAt(3), 56);
    BOOST_CHECK_EQUAL(path1.GapAt(0), 0);
    BOOST_CHECK_EQUAL(path1.GapAt(1), 9);
    BOOST_CHECK_EQUAL(path1.GapAt(2), 0);
    BOOST_CHECK_EQUAL(path1.GapAt(3), -1);
    BOOST_CHECK_EQUAL(path1[1], e2);
    BOOST_CHECK_EQUAL(path1[2], e3);
    BOOST_CHECK_EQUAL(path1[3], e4);
}

BOOST_AUTO_TEST_CASE( BidirectionalPathEquals ) {
    Graph g(13);
    IdTrackHandler<Graph> int_ids(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g, int_ids);
    EdgeId start = *g.ConstEdgeBegin();

    BidirectionalPath path1(g);
    BidirectionalPath path2(g);

    // 98 26 145 70 3 139 139
    EdgeId e1 = g.conjugate(start);
    EdgeId e2 = *(g.OutgoingEdges(g.EdgeEnd(e1)).begin());
    EdgeId e3 = *(g.OutgoingEdges(g.EdgeEnd(e2)).begin());
    EdgeId e4 = *(g.OutgoingEdges(g.EdgeEnd(e3)).begin());
    EdgeId e5 = *(g.OutgoingEdges(g.EdgeEnd(e4)).begin());
    auto it5 = g.OutgoingEdges(g.EdgeEnd(e5)).begin();
    EdgeId e7 = *it5;
    INFO(g.int_id(e7));
    ++it5;
    EdgeId e6 = *it5;
    INFO(g.int_id(e6));

    path1.PushBack(e1);
    path1.PushBack(e2);
    path1.PushBack(e3);
    path1.PushBack(e4);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e7);
}

BOOST_AUTO_TEST_CASE( BidirectionalPathSubpaths ) {
    Graph g(13);
    IdTrackHandler<Graph> int_ids(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g, int_ids);
    EdgeId start = *g.ConstEdgeBegin();

    BidirectionalPath path1(g);
    BidirectionalPath path2(g);

    // 98 26 145 70 3 139 139
    EdgeId e1 = g.conjugate(start);
    EdgeId e2 = *(g.OutgoingEdges(g.EdgeEnd(e1)).begin());
    EdgeId e3 = *(g.OutgoingEdges(g.EdgeEnd(e2)).begin());
    EdgeId e4 = *(g.OutgoingEdges(g.EdgeEnd(e3)).begin());
    EdgeId e5 = *(g.OutgoingEdges(g.EdgeEnd(e4)).begin());
    auto it5 = g.OutgoingEdges(g.EdgeEnd(e5)).begin();
    EdgeId e7 = *it5;
    INFO(g.int_id(e7));
    ++it5;
    EdgeId e6 = *it5;
    INFO(g.int_id(e6));

    path1.PushBack(e1);
    path1.PushBack(e2);
    path1.PushBack(e3);
    path1.PushBack(e4);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e7);
}

BOOST_AUTO_TEST_CASE( BidirectionalPathAddRemoveConjugate ) {
    Graph g(13);
    IdTrackHandler<Graph> int_ids(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g, int_ids);
    EdgeId start = *g.ConstEdgeBegin();

    BidirectionalPath p(g);
    // 98 26 145 70
    EdgeId e1 = g.conjugate(start);
    EdgeId e2 = *(g.OutgoingEdges(g.EdgeEnd(e1)).begin());
    EdgeId e3 = *(g.OutgoingEdges(g.EdgeEnd(e2)).begin());
    EdgeId e4 = *(g.OutgoingEdges(g.EdgeEnd(e3)).begin());

    BidirectionalPath cp(g);
    cp.Subscribe(&p);
    p.Subscribe(&cp);

    BOOST_CHECK_EQUAL(p.Conjugate().Empty(), true);
    BOOST_CHECK_EQUAL(cp.Conjugate().Empty(), true);

    p.PushBack(e1, 10);

}


BOOST_AUTO_TEST_CASE( BidirectionalPathSearch ) {
    Graph g(13);
    IdTrackHandler<Graph> int_ids(g);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g, int_ids);
    EdgeId start = *g.ConstEdgeBegin();

    BidirectionalPath path1(g);

    // 98 26 145 70
    EdgeId e1 = g.conjugate(start);
    EdgeId e2 = *(g.OutgoingEdges(g.EdgeEnd(e1)).begin());
    EdgeId e3 = *(g.OutgoingEdges(g.EdgeEnd(e2)).begin());
    EdgeId e4 = *(g.OutgoingEdges(g.EdgeEnd(e3)).begin());
}


BOOST_AUTO_TEST_SUITE_END()

}

