//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * path_extend.hpp
 *
 *  Created on: Nov 29, 2013
 *      Author: andrey
 */


#pragma once

#include <boost/test/unit_test.hpp>

#include "test_utils.hpp"
#include "modules/path_extend/path_visualizer.hpp"
#include "modules/path_extend/pe_utils.hpp"
namespace path_extend {

BOOST_FIXTURE_TEST_SUITE(path_extend_basic, fs::TmpFolderFixture)

BOOST_AUTO_TEST_CASE( BidirectionalPathConstructor ) {
    Graph g(13);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g);
    EdgeId start = *g.ConstEdgeBegin();

    BidirectionalPath path1(g);
    BOOST_CHECK_EQUAL(path1.Size(), 0);
    BOOST_CHECK_EQUAL(path1.Length(), 0);
    BOOST_CHECK_EQUAL(path1.Empty(), true);

    BidirectionalPath path2(g, start);
    BOOST_CHECK_EQUAL(path2.Size(), 1);
    BOOST_CHECK_EQUAL(path2.Length(), 426);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 426);
    BOOST_CHECK_EQUAL(path2.GapAt(0).gap, 0);
    BOOST_CHECK_EQUAL(path2[0], start);
    BOOST_CHECK_EQUAL(path2.At(0), start);

    EdgeId e1 = g.conjugate(start);
    EdgeId e2 = *(g.OutgoingEdges(g.EdgeEnd(e1)).begin());
    std::vector<EdgeId> v;
    v.push_back(e1);
    v.push_back(e2);
    BidirectionalPath path3(g, v);
    BOOST_CHECK_EQUAL(path3.Size(), 2);
    BOOST_CHECK_EQUAL(path3.Length(), 427);
    BOOST_CHECK_EQUAL(path3.LengthAt(0), 427);
    BOOST_CHECK_EQUAL(path3.LengthAt(1), 1);
    BOOST_CHECK_EQUAL(path3.GapAt(1).gap, 0);
    BOOST_CHECK_EQUAL(path3[0], e1);
    BOOST_CHECK_EQUAL(path3.At(1), e2);

    BidirectionalPath path4(path3);
    BOOST_CHECK_EQUAL(path4.Size(), 2);
    BOOST_CHECK_EQUAL(path4.Length(), 427);
    BOOST_CHECK_EQUAL(path4.LengthAt(0), 427);
    BOOST_CHECK_EQUAL(path4.LengthAt(1), 1);
    BOOST_CHECK_EQUAL(path4.GapAt(1).gap, 0);
    BOOST_CHECK_EQUAL(path4[0], e1);
    BOOST_CHECK_EQUAL(path4.At(1), e2);

    path2.Clear();
    BOOST_CHECK_EQUAL(path2.Empty(), true);

    path4.Clear();
    BOOST_CHECK_EQUAL(path4.Empty(), true);
}

BOOST_AUTO_TEST_CASE( BidirectionalPathAddRemove ) {
    Graph g(13);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g);
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
    BOOST_CHECK_EQUAL(path1.GapAt(1).gap, 0);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    path1.PushBack(e3);
    BOOST_CHECK_EQUAL(path1.Size(), 3);
    BOOST_CHECK_EQUAL(path1.Length(), 1006);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 1006);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 580);
    BOOST_CHECK_EQUAL(path1.LengthAt(2), 579);
    BOOST_CHECK_EQUAL(path1.GapAt(0).gap, 0);
    BOOST_CHECK_EQUAL(path1.GapAt(2).gap, 0);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1[2], e3);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    path1.PopBack();
    BOOST_CHECK_EQUAL(path1.Size(), 2);
    BOOST_CHECK_EQUAL(path1.Length(), 427);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 427);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 1);
    BOOST_CHECK_EQUAL(path1.GapAt(1).gap, 0);
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
    BOOST_CHECK_EQUAL(path1.GapAt(0).gap, 0);
    BOOST_CHECK_EQUAL(path1.GapAt(2).gap, 0);
    BOOST_CHECK_EQUAL(path1.GapAt(3).gap, 0);
    BOOST_CHECK_EQUAL(path1[1], e2);
    BOOST_CHECK_EQUAL(path1[2], e3);
    BOOST_CHECK_EQUAL(path1[3], e4);

    path1.PopBack(1);
    BOOST_CHECK_EQUAL(path1.Size(), 3);
    BOOST_CHECK_EQUAL(path1.Length(), 1006);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 1006);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 580);
    BOOST_CHECK_EQUAL(path1.LengthAt(2), 579);
    BOOST_CHECK_EQUAL(path1.GapAt(0).gap, 0);
    BOOST_CHECK_EQUAL(path1.GapAt(2).gap, 0);
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
    BOOST_CHECK_EQUAL(path1.GapAt(1).gap, 0);
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

    path2.PushBack(e3, Gap(0));
    BOOST_CHECK_EQUAL(path2.Empty(), false);
    BOOST_CHECK_EQUAL(path2.Size(), 1);
    BOOST_CHECK_EQUAL(path2.Length(), 579);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 579);

    path2.PushBack(e4, Gap(0));
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
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g);
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

    path1.PushBack(e2, Gap(10));
    BOOST_CHECK_EQUAL(path1.Size(), 2);
    BOOST_CHECK_EQUAL(path1.Length(), 437);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 437);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 1);
    BOOST_CHECK_EQUAL(path1.GapAt(1).gap, 10);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    path1.PushBack(e3, Gap(10));
    BOOST_CHECK_EQUAL(path1.Size(), 3);
    BOOST_CHECK_EQUAL(path1.Length(), 1026);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 1026);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 590);
    BOOST_CHECK_EQUAL(path1.LengthAt(2), 579);
    BOOST_CHECK_EQUAL(path1.GapAt(0).gap, 0);
    BOOST_CHECK_EQUAL(path1.GapAt(2).gap, 10);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1[2], e3);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    path1.PopBack();
    BOOST_CHECK_EQUAL(path1.Size(), 2);
    BOOST_CHECK_EQUAL(path1.Length(), 437);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 437);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 1);
    BOOST_CHECK_EQUAL(path1.GapAt(1).gap, 10);
    BOOST_CHECK_EQUAL(path1[0], e1);
    BOOST_CHECK_EQUAL(path1.At(1), e2);

    BidirectionalPath path2(g);
    path2.PushBack(e3);
    BOOST_CHECK_EQUAL(path2.Size(), 1);
    BOOST_CHECK_EQUAL(path2.Length(), 579);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 579);

    path2.PushBack(e4, Gap(100));
    BOOST_CHECK_EQUAL(path2.Empty(), false);
    BOOST_CHECK_EQUAL(path2.Size(), 2);
    BOOST_CHECK_EQUAL(path2.Length(), 735);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 735);
    BOOST_CHECK_EQUAL(path2.LengthAt(1), 56);
    BOOST_CHECK_EQUAL(path2.GapAt(0).gap, 0);
    BOOST_CHECK_EQUAL(path2.GapAt(1).gap, 100);
    BOOST_CHECK_EQUAL(path2[0], e3);
    BOOST_CHECK_EQUAL(path2[1], e4);

    path1.PushBack(path2, Gap(10));
    BOOST_CHECK_EQUAL(path1.Size(), 4);
    BOOST_CHECK_EQUAL(path1.Length(), 1182);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 1182);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 746);
    BOOST_CHECK_EQUAL(path1.LengthAt(2), 735);
    BOOST_CHECK_EQUAL(path1.LengthAt(3), 56);
    BOOST_CHECK_EQUAL(path1.GapAt(0).gap, 0);
    BOOST_CHECK_EQUAL(path1.GapAt(1).gap, 10);
    BOOST_CHECK_EQUAL(path1.GapAt(2).gap, 10);
    BOOST_CHECK_EQUAL(path1.GapAt(3).gap, 100);
    BOOST_CHECK_EQUAL(path1[1], e2);
    BOOST_CHECK_EQUAL(path1[2], e3);
    BOOST_CHECK_EQUAL(path1[3], e4);

    path1.PopBack(1);
    BOOST_CHECK_EQUAL(path1.Size(), 3);
    BOOST_CHECK_EQUAL(path1.Length(), 1026);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 1026);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 590);
    BOOST_CHECK_EQUAL(path1.LengthAt(2), 579);
    BOOST_CHECK_EQUAL(path1.GapAt(0).gap, 0);
    BOOST_CHECK_EQUAL(path1.GapAt(1).gap, 10);
    BOOST_CHECK_EQUAL(path1.GapAt(2).gap, 10);
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

    path1.PushBack(e2, Gap(-10));
    BOOST_CHECK_EQUAL(path1.Size(), 2);
    BOOST_CHECK_EQUAL(path1.Length(), 417);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 417);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 1);
    BOOST_CHECK_EQUAL(path1.GapAt(0).gap, 0);
    BOOST_CHECK_EQUAL(path1.GapAt(1).gap, -10);
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

    path2.PushBack(e4, Gap(-1));
    BOOST_CHECK_EQUAL(path2.Empty(), false);
    BOOST_CHECK_EQUAL(path2.Size(), 2);
    BOOST_CHECK_EQUAL(path2.Length(), 634);
    BOOST_CHECK_EQUAL(path2.LengthAt(0), 634);
    BOOST_CHECK_EQUAL(path2.LengthAt(1), 56);
    BOOST_CHECK_EQUAL(path2.GapAt(0).gap, 0);
    BOOST_CHECK_EQUAL(path2.GapAt(1).gap, -1);
    BOOST_CHECK_EQUAL(path2[0], e3);
    BOOST_CHECK_EQUAL(path2[1], e4);

    path1.Clear();
    BOOST_CHECK_EQUAL(path1.Empty(), true);
    path1.PushBack(e1, Gap(0));
    path1.PushBack(e2, Gap(9));
    path1.PushBack(path2);
    BOOST_CHECK_EQUAL(path1.Size(), 4);
    BOOST_CHECK_EQUAL(path1.Length(), 1070);
    BOOST_CHECK_EQUAL(path1.LengthAt(0), 1070);
    BOOST_CHECK_EQUAL(path1.LengthAt(1), 635);
    BOOST_CHECK_EQUAL(path1.LengthAt(2), 634);
    BOOST_CHECK_EQUAL(path1.LengthAt(3), 56);
    BOOST_CHECK_EQUAL(path1.GapAt(0).gap, 0);
    BOOST_CHECK_EQUAL(path1.GapAt(1).gap, 9);
    BOOST_CHECK_EQUAL(path1.GapAt(2).gap, 0);
    BOOST_CHECK_EQUAL(path1.GapAt(3).gap, -1);
    BOOST_CHECK_EQUAL(path1[1], e2);
    BOOST_CHECK_EQUAL(path1[2], e3);
    BOOST_CHECK_EQUAL(path1[3], e4);
}

BOOST_AUTO_TEST_CASE( BidirectionalPathEquals ) {
    Graph g(13);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g);
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
    ++it5;
    EdgeId e6 = *it5;

    path1.PushBack(e1);
    path1.PushBack(e2);
    path1.PushBack(e3);
    path1.PushBack(e4);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e6);

    path2.PushBack(e5);
    path2.PushBack(e6);
    path2.PushBack(e5);
    BOOST_CHECK_EQUAL(path1 == path2, false);

    path2.Clear();
    path2.PushBack(e1);
    path2.PushBack(e2);
    path2.PushBack(e3);
    path2.PushBack(e4);
    path2.PushBack(e5);
    path2.PushBack(e6);
    path2.PushBack(e5);
    BOOST_CHECK_EQUAL(path1 != path2, true);

    path2.Clear();
    path2.PushBack(e1);
    path2.PushBack(e2);
    path2.PushBack(e3);
    path2.PushBack(e4);
    path2.PushBack(e5);
    path2.PushBack(e6);
    path2.PushBack(e5);
    path2.PushBack(e6);
    BOOST_CHECK_EQUAL(path1.Equal(path2), true);
}

BOOST_AUTO_TEST_CASE( BidirectionalPathSubpaths ) {
    Graph g(13);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g);
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
    ++it5;
    EdgeId e6 = *it5;

    path1.PushBack(e1);
    path1.PushBack(e2);
    path1.PushBack(e3);
    path1.PushBack(e4);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5); //6
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e7);
    BOOST_CHECK_EQUAL(path1.SubPath(0) == path1, true);

    path2.PushBack(e5);
    path2.PushBack(e6);
    path2.PushBack(e5);
    BOOST_CHECK_EQUAL(path1.SubPath(4, 7) == path2, true);
    BOOST_CHECK_EQUAL(path1.SubPath(4, 6) == path2, false);
    BOOST_CHECK_EQUAL(path1.SubPath(4, 7) == path1.SubPath(6, 9), true);
    BOOST_CHECK_EQUAL(path1.SubPath(4, 7) == path1.SubPath(8, 11), true);
    BOOST_CHECK_EQUAL(path1.SubPath(4, 7) != path1.SubPath(8, 12), true);

    path2.Clear();
    path2.PushBack(e5);
    path2.PushBack(e6);
    path2.PushBack(e5);
    path2.PushBack(e7);
    BOOST_CHECK_EQUAL(path1.SubPath(8) == path2, true);
    BOOST_CHECK_EQUAL(path1.SubPath(8, 12) == path2, true);

    path1.Clear();
    path1.PushBack(e1);
    path1.PushBack(e2, Gap(10));
    path1.PushBack(e3, Gap(10));
    BidirectionalPath path3(path1.SubPath(1));
    BOOST_CHECK_EQUAL(path3.Size(), 2);
    BOOST_CHECK_EQUAL(path3.Length(), 590);
    BOOST_CHECK_EQUAL(path3.LengthAt(0), 590);
    BOOST_CHECK_EQUAL(path3.LengthAt(1), 579);
    BOOST_CHECK_EQUAL(path3.GapAt(0).gap, 0);
    BOOST_CHECK_EQUAL(path3.GapAt(1).gap, 10);
    BOOST_CHECK_EQUAL(path3[0], e2);
    BOOST_CHECK_EQUAL(path3[1], e3);
}

BOOST_AUTO_TEST_CASE( BidirectionalPathAddRemoveConjugate ) {
    Graph g(13);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g);
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

    p.PushBack(e1);
    BOOST_CHECK_EQUAL(cp.Conjugate() == p, true);

    p.PushBack(e2, Gap(100));
    BOOST_CHECK_EQUAL(cp.Conjugate() == p, true);

    p.PushBack(e3, Gap(10));
    BOOST_CHECK_EQUAL(cp.Conjugate() == p, true);

    p.PushBack(e4, Gap(10));
    BOOST_CHECK_EQUAL(cp.Conjugate() == p, true);
    BOOST_CHECK_EQUAL(cp.Front(), g.conjugate(e4));
    BOOST_CHECK_EQUAL(cp[1], g.conjugate(e3));
    BOOST_CHECK_EQUAL(cp[2], g.conjugate(e2));
    BOOST_CHECK_EQUAL(cp.Back(), g.conjugate(e1));
    BOOST_CHECK_EQUAL(cp.GapAt(0).gap, 0);
    BOOST_CHECK_EQUAL(cp.GapAt(1).gap, 10);
    BOOST_CHECK_EQUAL(cp.GapAt(2).gap, 10);
    BOOST_CHECK_EQUAL(cp.GapAt(3).gap, 100);
    BOOST_CHECK_EQUAL(cp.LengthAt(0), 1182);
    BOOST_CHECK_EQUAL(cp.LengthAt(1), 1116);
    BOOST_CHECK_EQUAL(cp.LengthAt(2), 527);
    BOOST_CHECK_EQUAL(cp.LengthAt(3), 426);

    p.PopBack(2);
    BOOST_CHECK_EQUAL(cp.Conjugate() == p, true);
    BOOST_CHECK_EQUAL(cp.Back(), g.conjugate(e1));
    BOOST_CHECK_EQUAL(cp.Front(), g.conjugate(e2));

    p.PopBack();
    BOOST_CHECK_EQUAL(cp.Conjugate() == p, true);
    BOOST_CHECK_EQUAL(cp.Back(), g.conjugate(e1));
    BOOST_CHECK_EQUAL(cp.Front(), g.conjugate(e1));

    p.PushBack(e2, Gap(100));
    BOOST_CHECK_EQUAL(cp.Conjugate() == p, true);
    BOOST_CHECK_EQUAL(cp.Back(), g.conjugate(e1));
    BOOST_CHECK_EQUAL(cp.Front(), g.conjugate(e2));

    p.PushBack(e3, Gap(10));
    BOOST_CHECK_EQUAL(cp.Conjugate() == p, true);
    BOOST_CHECK_EQUAL(cp.Back(), g.conjugate(e1));
    BOOST_CHECK_EQUAL(cp.Front(), g.conjugate(e3));

    p.PushBack(e4, Gap(10));
    BOOST_CHECK_EQUAL(cp.Conjugate() == p, true);
    BOOST_CHECK_EQUAL(cp.Front(), g.conjugate(e4));
    BOOST_CHECK_EQUAL(cp[1], g.conjugate(e3));
    BOOST_CHECK_EQUAL(cp[2], g.conjugate(e2));
    BOOST_CHECK_EQUAL(cp.Back(), g.conjugate(e1));
    BOOST_CHECK_EQUAL(cp.GapAt(0).gap, 0);
    BOOST_CHECK_EQUAL(cp.GapAt(1).gap, 10);
    BOOST_CHECK_EQUAL(cp.GapAt(2).gap, 10);
    BOOST_CHECK_EQUAL(cp.GapAt(3).gap, 100);
    BOOST_CHECK_EQUAL(cp.LengthAt(0), 1182);
    BOOST_CHECK_EQUAL(cp.LengthAt(1), 1116);
    BOOST_CHECK_EQUAL(cp.LengthAt(2), 527);
    BOOST_CHECK_EQUAL(cp.LengthAt(3), 426);
}


BOOST_AUTO_TEST_CASE( BidirectionalPathSearch ) {
    Graph g(13);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g);
    EdgeId start = *g.ConstEdgeBegin();

    BidirectionalPath path1(g);

    // 98 26 145 70 3 139 139
    EdgeId e1 = g.conjugate(start);
    EdgeId e2 = *(g.OutgoingEdges(g.EdgeEnd(e1)).begin());
    EdgeId e3 = *(g.OutgoingEdges(g.EdgeEnd(e2)).begin());
    EdgeId e4 = *(g.OutgoingEdges(g.EdgeEnd(e3)).begin());
    EdgeId e5 = *(g.OutgoingEdges(g.EdgeEnd(e4)).begin());
    auto it5 = g.OutgoingEdges(g.EdgeEnd(e5)).begin();
    EdgeId e7 = *it5;
    ++it5;
    EdgeId e6 = *it5;

    path1.PushBack(e1);
    path1.PushBack(e2);
    path1.PushBack(e3);
    path1.PushBack(e4);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5); //6
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e7);

    BOOST_CHECK_EQUAL(path1.FindFirst(e1), 0);
    BOOST_CHECK_EQUAL(path1.FindFirst(e5), 4);
    BOOST_CHECK_EQUAL(path1.FindFirst(e6), 5);
    BOOST_CHECK_EQUAL(path1.FindLast(e1), 0);
    BOOST_CHECK_EQUAL(path1.FindLast(e5), 10);
    BOOST_CHECK_EQUAL(path1.FindLast(e6), 9);

    auto v = path1.FindAll(e2);
    BOOST_CHECK_EQUAL(v.size(), 1);
    BOOST_CHECK_EQUAL(v[0], 1);
    v = path1.FindAll(e6);
    BOOST_CHECK_EQUAL(v.size(), 3);
    BOOST_CHECK_EQUAL(v[1], 7);
    v = path1.FindAll(e5);
    BOOST_CHECK_EQUAL(v.size(), 4);
    BOOST_CHECK_EQUAL(v[1], 6);
    BOOST_CHECK_EQUAL(v[2], 8);
    v = path1.FindAll(start);
    BOOST_CHECK_EQUAL(v.size(), 0);
    BOOST_CHECK_EQUAL(path1.FindFirst(start), -1);
    BOOST_CHECK_EQUAL(path1.FindLast(start), -1);

    BOOST_CHECK_EQUAL(path1.CompareFrom(0, path1), true);
    BOOST_CHECK_EQUAL(path1.CompareFrom(1, path1), false);
    BOOST_CHECK_EQUAL(path1.CompareFrom(20, path1), false);
    BOOST_CHECK_EQUAL(path1.FindFirst(path1), 0);
    BOOST_CHECK_EQUAL(path1.FindLast(path1), 0);

    BidirectionalPath path2 = path1.SubPath(1, 4);
    BidirectionalPath path3 = path1.SubPath(4, 7);;
    BidirectionalPath path4 = path1.SubPath(9);
    BidirectionalPath path5(g);
    path5.PushBack(e4);
    path5.PushBack(e5);
    path5.PushBack(e6);
    path5.PushBack(e5);
    path5.PushBack(e7);

    BOOST_CHECK_EQUAL(path1.CompareFrom(0, path2), false);
    BOOST_CHECK_EQUAL(path1.CompareFrom(1, path2), true);
    BOOST_CHECK_EQUAL(path1.CompareFrom(2, path2), false);
    BOOST_CHECK_EQUAL(path1.FindFirst(path2), 1);
    BOOST_CHECK_EQUAL(path1.FindLast(path2), 1);
    //BOOST_CHECK_EQUAL(path1.Contains(path2), true);

    BOOST_CHECK_EQUAL(path1.CompareFrom(4, path3), true);
    BOOST_CHECK_EQUAL(path1.CompareFrom(5, path3), false);
    BOOST_CHECK_EQUAL(path1.CompareFrom(6, path3), true);
    BOOST_CHECK_EQUAL(path1.FindFirst(path3), 4);
    BOOST_CHECK_EQUAL(path1.FindLast(path3), 8);
    //BOOST_CHECK_EQUAL(path1.Contains(path3), true);

    BOOST_CHECK_EQUAL(path1.CompareFrom(9, path4), true);
    BOOST_CHECK_EQUAL(path1.CompareFrom(11, path4), false);
    BOOST_CHECK_EQUAL(path1.FindFirst(path4), 9);
    BOOST_CHECK_EQUAL(path1.FindLast(path4), 9);
    //BOOST_CHECK_EQUAL(path1.Contains(path4), true);

    BOOST_CHECK_EQUAL(path1.CompareFrom(0, path5), false);
    BOOST_CHECK_EQUAL(path1.CompareFrom(3, path5), false);
    BOOST_CHECK_EQUAL(path1.CompareFrom(5, path5), false);
    BOOST_CHECK_EQUAL(path1.FindFirst(path5), -1);
    BOOST_CHECK_EQUAL(path1.FindLast(path5), -1);
    //BOOST_CHECK_EQUAL(path1.Contains(path5), false);
}


BOOST_AUTO_TEST_CASE( BidirectionalPathLoopDetector ) {
    Graph g(13);
    graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g);
    EdgeId start = *g.ConstEdgeBegin();

    BidirectionalPath path1(g);
    GraphCoverageMap cover_map(g);
    LoopDetector loop_detect(&path1, cover_map);
    path1.Subscribe(&cover_map);

    // 98 26 145 70 3 139 139
    EdgeId e1 = g.conjugate(start);
    EdgeId e2 = *(g.OutgoingEdges(g.EdgeEnd(e1)).begin());
    EdgeId e3 = *(g.OutgoingEdges(g.EdgeEnd(e2)).begin());
    EdgeId e4 = *(g.OutgoingEdges(g.EdgeEnd(e3)).begin());
    EdgeId e5 = *(g.OutgoingEdges(g.EdgeEnd(e4)).begin());
    auto it5 = g.OutgoingEdges(g.EdgeEnd(e5)).begin();
    EdgeId e7 = *it5;
    ++it5;
    EdgeId e6 = *it5;

    path1.PushBack(e1);
    path1.PushBack(e2);
    path1.PushBack(e3);
    path1.PushBack(e4);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5); //6
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e6);
    size_t skip_identical = 0;
    BOOST_CHECK_EQUAL(loop_detect.IsCycled(4, skip_identical), false);
    BOOST_CHECK_EQUAL(loop_detect.IsCycled(3, skip_identical), true);
    BOOST_CHECK_EQUAL(loop_detect.IsCycled(2, skip_identical), true);
    loop_detect.RemoveLoop(skip_identical, false);
    BOOST_CHECK_EQUAL(path1.Size(), 6);
    BOOST_CHECK_EQUAL(path1.Back(), e6);

    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);

    skip_identical = 0;
    BOOST_CHECK_EQUAL(loop_detect.IsCycled(4, skip_identical), false);
    skip_identical = 0;
    BOOST_CHECK_EQUAL(loop_detect.IsCycled(3, skip_identical), true);
    skip_identical = 0;
    BOOST_CHECK_EQUAL(loop_detect.IsCycled(2, skip_identical), true);

    loop_detect.RemoveLoop(skip_identical, true);
    BOOST_CHECK_EQUAL(path1.Size(), 5);
    BOOST_CHECK_EQUAL(path1.Back(), e5);

    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e7);

    skip_identical = 0;
    BOOST_CHECK_EQUAL(loop_detect.IsCycled(4, skip_identical), false);
    skip_identical = 0;
    BOOST_CHECK_EQUAL(loop_detect.IsCycled(3, skip_identical), false);
    skip_identical = 0;
    BOOST_CHECK_EQUAL(loop_detect.IsCycled(2, skip_identical), false);

    BOOST_CHECK_EQUAL(loop_detect.EdgesToRemove(skip_identical, false), 0);
    loop_detect.RemoveLoop(skip_identical, false);
    BOOST_CHECK_EQUAL(path1.Size(), 12);
    BOOST_CHECK_EQUAL(path1.Back(), e7);
}


BOOST_AUTO_TEST_SUITE_END()

}

