//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************


#include "modules/path_extend/path_visualizer.hpp"
#include "modules/path_extend/pe_utils.hpp"

#include "graphio.hpp"

#include <gtest/gtest.h>

using namespace path_extend;
using namespace debruijn_graph;

TEST( PathExtend, BidirectionalPathConstructor ) {
    Graph g(13);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g));
    EdgeId start = *g.ConstEdgeBegin();

    BidirectionalPath path1(g);
    EXPECT_EQ(path1.Size(), 0);
    EXPECT_EQ(path1.Length(), 0);
    EXPECT_EQ(path1.Empty(), true);

    BidirectionalPath path2(g, start);
    EXPECT_EQ(path2.Size(), 1);
    EXPECT_EQ(path2.Length(), 426);
    EXPECT_EQ(path2.LengthAt(0), 426);
    EXPECT_EQ(path2.GapAt(0).gap, 0);
    EXPECT_EQ(path2[0], start);
    EXPECT_EQ(path2.At(0), start);

    EdgeId e1 = g.conjugate(start);
    EdgeId e2 = *(g.OutgoingEdges(g.EdgeEnd(e1)).begin());
    std::vector<EdgeId> v;
    v.push_back(e1);
    v.push_back(e2);
    BidirectionalPath path3(g, v);
    EXPECT_EQ(path3.Size(), 2);
    EXPECT_EQ(path3.Length(), 427);
    EXPECT_EQ(path3.LengthAt(0), 427);
    EXPECT_EQ(path3.LengthAt(1), 1);
    EXPECT_EQ(path3.GapAt(1).gap, 0);
    EXPECT_EQ(path3[0], e1);
    EXPECT_EQ(path3.At(1), e2);

    BidirectionalPath path4(path3);
    EXPECT_EQ(path4.Size(), 2);
    EXPECT_EQ(path4.Length(), 427);
    EXPECT_EQ(path4.LengthAt(0), 427);
    EXPECT_EQ(path4.LengthAt(1), 1);
    EXPECT_EQ(path4.GapAt(1).gap, 0);
    EXPECT_EQ(path4[0], e1);
    EXPECT_EQ(path4.At(1), e2);

    path2.Clear();
    EXPECT_EQ(path2.Empty(), true);

    path4.Clear();
    EXPECT_EQ(path4.Empty(), true);
}

TEST( PathExtend, BidirectionalPathAddRemove ) {
    Graph g(13);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g));
    EdgeId start = *g.ConstEdgeBegin();

    BidirectionalPath path1(g);
    EXPECT_EQ(path1.Size(), 0);
    EXPECT_EQ(path1.Length(), 0);

    // 98 26 145 70
    EdgeId e1 = g.conjugate(start);
    EdgeId e2 = *(g.OutgoingEdges(g.EdgeEnd(e1)).begin());
    EdgeId e3 = *(g.OutgoingEdges(g.EdgeEnd(e2)).begin());
    EdgeId e4 = *(g.OutgoingEdges(g.EdgeEnd(e3)).begin());

    path1.PushBack(e1);
    EXPECT_EQ(path1.Size(), 1);
    EXPECT_EQ(path1.Length(), 426);
    EXPECT_EQ(path1.LengthAt(0), 426);

    path1.PushBack(e2);
    EXPECT_EQ(path1.Size(), 2);
    EXPECT_EQ(path1.Length(), 427);
    EXPECT_EQ(path1.LengthAt(0), 427);
    EXPECT_EQ(path1.LengthAt(1), 1);
    EXPECT_EQ(path1.GapAt(1).gap, 0);
    EXPECT_EQ(path1[0], e1);
    EXPECT_EQ(path1.At(1), e2);

    path1.PushBack(e3);
    EXPECT_EQ(path1.Size(), 3);
    EXPECT_EQ(path1.Length(), 1006);
    EXPECT_EQ(path1.LengthAt(0), 1006);
    EXPECT_EQ(path1.LengthAt(1), 580);
    EXPECT_EQ(path1.LengthAt(2), 579);
    EXPECT_EQ(path1.GapAt(0).gap, 0);
    EXPECT_EQ(path1.GapAt(2).gap, 0);
    EXPECT_EQ(path1[0], e1);
    EXPECT_EQ(path1[2], e3);
    EXPECT_EQ(path1.At(1), e2);

    path1.PopBack();
    EXPECT_EQ(path1.Size(), 2);
    EXPECT_EQ(path1.Length(), 427);
    EXPECT_EQ(path1.LengthAt(0), 427);
    EXPECT_EQ(path1.LengthAt(1), 1);
    EXPECT_EQ(path1.GapAt(1).gap, 0);
    EXPECT_EQ(path1[0], e1);
    EXPECT_EQ(path1.At(1), e2);

    BidirectionalPath path2(g);
    path2.PushBack(e3);
    EXPECT_EQ(path2.Size(), 1);
    EXPECT_EQ(path2.Length(), 579);
    EXPECT_EQ(path2.LengthAt(0), 579);

    path2.PushBack(e4);
    EXPECT_EQ(path2.Empty(), false);
    EXPECT_EQ(path2.Size(), 2);
    EXPECT_EQ(path2.Length(), 635);
    EXPECT_EQ(path2.LengthAt(0), 635);
    EXPECT_EQ(path2.LengthAt(1), 56);
    EXPECT_EQ(path2[0], e3);
    EXPECT_EQ(path2[1], e4);

    path1.PushBack(path2);
    EXPECT_EQ(path1.Size(), 4);
    EXPECT_EQ(path1.Length(), 1062);
    EXPECT_EQ(path1.LengthAt(0), 1062);
    EXPECT_EQ(path1.LengthAt(2), 635);
    EXPECT_EQ(path1.LengthAt(3), 56);
    EXPECT_EQ(path1.GapAt(0).gap, 0);
    EXPECT_EQ(path1.GapAt(2).gap, 0);
    EXPECT_EQ(path1.GapAt(3).gap, 0);
    EXPECT_EQ(path1[1], e2);
    EXPECT_EQ(path1[2], e3);
    EXPECT_EQ(path1[3], e4);

    path1.PopBack(1);
    EXPECT_EQ(path1.Size(), 3);
    EXPECT_EQ(path1.Length(), 1006);
    EXPECT_EQ(path1.LengthAt(0), 1006);
    EXPECT_EQ(path1.LengthAt(1), 580);
    EXPECT_EQ(path1.LengthAt(2), 579);
    EXPECT_EQ(path1.GapAt(0).gap, 0);
    EXPECT_EQ(path1.GapAt(2).gap, 0);
    EXPECT_EQ(path1[0], e1);
    EXPECT_EQ(path1[2], e3);
    EXPECT_EQ(path1.At(1), e2);

    path1.PopBack(3);
    EXPECT_EQ(path1.Length(), 0);
    EXPECT_EQ(path1.Empty(), true);

    path1.PushBack(e1);
    EXPECT_EQ(path1.Size(), 1);
    EXPECT_EQ(path1.Length(), 426);
    EXPECT_EQ(path1.LengthAt(0), 426);

    path1.PushBack(e2);
    EXPECT_EQ(path1.Size(), 2);
    EXPECT_EQ(path1.Length(), 427);
    EXPECT_EQ(path1.LengthAt(0), 427);
    EXPECT_EQ(path1.LengthAt(1), 1);
    EXPECT_EQ(path1.GapAt(1).gap, 0);
    EXPECT_EQ(path1[0], e1);
    EXPECT_EQ(path1.At(1), e2);

    path1.PopBack(10);
    EXPECT_EQ(path1.Length(), 0);
    EXPECT_EQ(path1.Empty(), true);

    path1.PopBack(10);
    EXPECT_EQ(path1.Length(), 0);
    EXPECT_EQ(path1.Empty(), true);

    path2.PopBack(1);
    EXPECT_EQ(path2.Size(), 1);
    EXPECT_EQ(path2.Length(), 579);
    EXPECT_EQ(path2.LengthAt(0), 579);

    path2.PopBack(2);
    EXPECT_EQ(path2.Length(), 0);
    EXPECT_EQ(path2.Empty(), true);

    path2.PopBack();
    EXPECT_EQ(path2.Length(), 0);
    EXPECT_EQ(path2.Empty(), true);

    path2.PushBack(e3, Gap(0));
    EXPECT_EQ(path2.Empty(), false);
    EXPECT_EQ(path2.Size(), 1);
    EXPECT_EQ(path2.Length(), 579);
    EXPECT_EQ(path2.LengthAt(0), 579);

    path2.PushBack(e4, Gap(0));
    EXPECT_EQ(path2.Empty(), false);
    EXPECT_EQ(path2.Size(), 2);
    EXPECT_EQ(path2.Length(), 635);
    EXPECT_EQ(path2.LengthAt(0), 635);
    EXPECT_EQ(path2.LengthAt(1), 56);
    EXPECT_EQ(path2[0], e3);
    EXPECT_EQ(path2[1], e4);

    path2.Clear();
    EXPECT_EQ(path2.Empty(), true);
}

TEST( PathExtend, BidirectionalPathAddRemoveGaps ) {
    Graph g(13);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g));
    EdgeId start = *g.ConstEdgeBegin();

    BidirectionalPath path1(g);

    // 98 26 145 70
    EdgeId e1 = g.conjugate(start);
    EdgeId e2 = *(g.OutgoingEdges(g.EdgeEnd(e1)).begin());
    EdgeId e3 = *(g.OutgoingEdges(g.EdgeEnd(e2)).begin());
    EdgeId e4 = *(g.OutgoingEdges(g.EdgeEnd(e3)).begin());

    path1.PushBack(e1);
    EXPECT_EQ(path1.Size(), 1);
    EXPECT_EQ(path1.Length(), 426);
    EXPECT_EQ(path1.LengthAt(0), 426);

    path1.PushBack(e2, Gap(10));
    EXPECT_EQ(path1.Size(), 2);
    EXPECT_EQ(path1.Length(), 437);
    EXPECT_EQ(path1.LengthAt(0), 437);
    EXPECT_EQ(path1.LengthAt(1), 1);
    EXPECT_EQ(path1.GapAt(1).gap, 10);
    EXPECT_EQ(path1[0], e1);
    EXPECT_EQ(path1.At(1), e2);

    path1.PushBack(e3, Gap(10));
    EXPECT_EQ(path1.Size(), 3);
    EXPECT_EQ(path1.Length(), 1026);
    EXPECT_EQ(path1.LengthAt(0), 1026);
    EXPECT_EQ(path1.LengthAt(1), 590);
    EXPECT_EQ(path1.LengthAt(2), 579);
    EXPECT_EQ(path1.GapAt(0).gap, 0);
    EXPECT_EQ(path1.GapAt(2).gap, 10);
    EXPECT_EQ(path1[0], e1);
    EXPECT_EQ(path1[2], e3);
    EXPECT_EQ(path1.At(1), e2);

    path1.PopBack();
    EXPECT_EQ(path1.Size(), 2);
    EXPECT_EQ(path1.Length(), 437);
    EXPECT_EQ(path1.LengthAt(0), 437);
    EXPECT_EQ(path1.LengthAt(1), 1);
    EXPECT_EQ(path1.GapAt(1).gap, 10);
    EXPECT_EQ(path1[0], e1);
    EXPECT_EQ(path1.At(1), e2);

    BidirectionalPath path2(g);
    path2.PushBack(e3);
    EXPECT_EQ(path2.Size(), 1);
    EXPECT_EQ(path2.Length(), 579);
    EXPECT_EQ(path2.LengthAt(0), 579);

    path2.PushBack(e4, Gap(100));
    EXPECT_EQ(path2.Empty(), false);
    EXPECT_EQ(path2.Size(), 2);
    EXPECT_EQ(path2.Length(), 735);
    EXPECT_EQ(path2.LengthAt(0), 735);
    EXPECT_EQ(path2.LengthAt(1), 56);
    EXPECT_EQ(path2.GapAt(0).gap, 0);
    EXPECT_EQ(path2.GapAt(1).gap, 100);
    EXPECT_EQ(path2[0], e3);
    EXPECT_EQ(path2[1], e4);

    path1.PushBack(path2, Gap(10));
    EXPECT_EQ(path1.Size(), 4);
    EXPECT_EQ(path1.Length(), 1182);
    EXPECT_EQ(path1.LengthAt(0), 1182);
    EXPECT_EQ(path1.LengthAt(1), 746);
    EXPECT_EQ(path1.LengthAt(2), 735);
    EXPECT_EQ(path1.LengthAt(3), 56);
    EXPECT_EQ(path1.GapAt(0).gap, 0);
    EXPECT_EQ(path1.GapAt(1).gap, 10);
    EXPECT_EQ(path1.GapAt(2).gap, 10);
    EXPECT_EQ(path1.GapAt(3).gap, 100);
    EXPECT_EQ(path1[1], e2);
    EXPECT_EQ(path1[2], e3);
    EXPECT_EQ(path1[3], e4);

    path1.PopBack(1);
    EXPECT_EQ(path1.Size(), 3);
    EXPECT_EQ(path1.Length(), 1026);
    EXPECT_EQ(path1.LengthAt(0), 1026);
    EXPECT_EQ(path1.LengthAt(1), 590);
    EXPECT_EQ(path1.LengthAt(2), 579);
    EXPECT_EQ(path1.GapAt(0).gap, 0);
    EXPECT_EQ(path1.GapAt(1).gap, 10);
    EXPECT_EQ(path1.GapAt(2).gap, 10);
    EXPECT_EQ(path1[0], e1);
    EXPECT_EQ(path1[2], e3);
    EXPECT_EQ(path1.At(1), e2);

    path1.PopBack(3);
    EXPECT_EQ(path1.Length(), 0);
    EXPECT_EQ(path1.Empty(), true);

    path1.PushBack(e1);
    EXPECT_EQ(path1.Size(), 1);
    EXPECT_EQ(path1.Length(), 426);
    EXPECT_EQ(path1.LengthAt(0), 426);

    path1.PushBack(e2, Gap(-10));
    EXPECT_EQ(path1.Size(), 2);
    EXPECT_EQ(path1.Length(), 417);
    EXPECT_EQ(path1.LengthAt(0), 417);
    EXPECT_EQ(path1.LengthAt(1), 1);
    EXPECT_EQ(path1.GapAt(0).gap, 0);
    EXPECT_EQ(path1.GapAt(1).gap, -10);
    EXPECT_EQ(path1[0], e1);
    EXPECT_EQ(path1.At(1), e2);

    path1.PopBack(10);
    EXPECT_EQ(path1.Length(), 0);
    EXPECT_EQ(path1.Empty(), true);

    path1.PopBack(10);
    EXPECT_EQ(path1.Length(), 0);
    EXPECT_EQ(path1.Empty(), true);

    path2.PopBack(1);
    EXPECT_EQ(path2.Size(), 1);
    EXPECT_EQ(path2.Length(), 579);
    EXPECT_EQ(path2.LengthAt(0), 579);

    path2.PopBack(2);
    EXPECT_EQ(path2.Length(), 0);
    EXPECT_EQ(path2.Empty(), true);

    path2.PopBack();
    EXPECT_EQ(path2.Length(), 0);
    EXPECT_EQ(path2.Empty(), true);

    path2.PushBack(e3);
    EXPECT_EQ(path2.Empty(), false);
    EXPECT_EQ(path2.Size(), 1);
    EXPECT_EQ(path2.Length(), 579);
    EXPECT_EQ(path2.LengthAt(0), 579);

    path2.PushBack(e4, Gap(-1));
    EXPECT_EQ(path2.Empty(), false);
    EXPECT_EQ(path2.Size(), 2);
    EXPECT_EQ(path2.Length(), 634);
    EXPECT_EQ(path2.LengthAt(0), 634);
    EXPECT_EQ(path2.LengthAt(1), 56);
    EXPECT_EQ(path2.GapAt(0).gap, 0);
    EXPECT_EQ(path2.GapAt(1).gap, -1);
    EXPECT_EQ(path2[0], e3);
    EXPECT_EQ(path2[1], e4);

    path1.Clear();
    EXPECT_EQ(path1.Empty(), true);
    path1.PushBack(e1, Gap(0));
    path1.PushBack(e2, Gap(9));
    path1.PushBack(path2);
    EXPECT_EQ(path1.Size(), 4);
    EXPECT_EQ(path1.Length(), 1070);
    EXPECT_EQ(path1.LengthAt(0), 1070);
    EXPECT_EQ(path1.LengthAt(1), 635);
    EXPECT_EQ(path1.LengthAt(2), 634);
    EXPECT_EQ(path1.LengthAt(3), 56);
    EXPECT_EQ(path1.GapAt(0).gap, 0);
    EXPECT_EQ(path1.GapAt(1).gap, 9);
    EXPECT_EQ(path1.GapAt(2).gap, 0);
    EXPECT_EQ(path1.GapAt(3).gap, -1);
    EXPECT_EQ(path1[1], e2);
    EXPECT_EQ(path1[2], e3);
    EXPECT_EQ(path1[3], e4);
}

TEST( PathExtend, BidirectionalPathEquals ) {
    Graph g(13);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g));
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
    EXPECT_NE(path1, path2);

    path2.Clear();
    path2.PushBack(e1);
    path2.PushBack(e2);
    path2.PushBack(e3);
    path2.PushBack(e4);
    path2.PushBack(e5);
    path2.PushBack(e6);
    path2.PushBack(e5);
    EXPECT_NE(path1, path2);

    path2.Clear();
    path2.PushBack(e1);
    path2.PushBack(e2);
    path2.PushBack(e3);
    path2.PushBack(e4);
    path2.PushBack(e5);
    path2.PushBack(e6);
    path2.PushBack(e5);
    path2.PushBack(e6);
    EXPECT_TRUE(path1 == path2);
}

TEST( PathExtend, BidirectionalPathSubpaths ) {
    Graph g(13);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g));
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
    EXPECT_EQ(path1.SubPath(0), path1);

    path2.PushBack(e5);
    path2.PushBack(e6);
    path2.PushBack(e5);
    EXPECT_EQ(path1.SubPath(4, 7), path2);
    EXPECT_NE(path1.SubPath(4, 6), path2);
    EXPECT_EQ(path1.SubPath(4, 7), path1.SubPath(6, 9));
    EXPECT_EQ(path1.SubPath(4, 7), path1.SubPath(8, 11));
    EXPECT_NE(path1.SubPath(4, 7), path1.SubPath(8, 12));

    path2.Clear();
    path2.PushBack(e5);
    path2.PushBack(e6);
    path2.PushBack(e5);
    path2.PushBack(e7);
    EXPECT_EQ(path1.SubPath(8), path2);
    EXPECT_EQ(path1.SubPath(8, 12), path2);

    path1.Clear();
    path1.PushBack(e1);
    path1.PushBack(e2, Gap(10));
    path1.PushBack(e3, Gap(10));
    BidirectionalPath path3(path1.SubPath(1));
    EXPECT_EQ(path3.Size(), 2);
    EXPECT_EQ(path3.Length(), 590);
    EXPECT_EQ(path3.LengthAt(0), 590);
    EXPECT_EQ(path3.LengthAt(1), 579);
    EXPECT_EQ(path3.GapAt(0).gap, 0);
    EXPECT_EQ(path3.GapAt(1).gap, 10);
    EXPECT_EQ(path3[0], e2);
    EXPECT_EQ(path3[1], e3);
}

TEST( PathExtend, BidirectionalPathAddRemoveConjugate ) {
    Graph g(13);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g));
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

    EXPECT_TRUE(p.Conjugate().Empty());
    EXPECT_TRUE(cp.Conjugate().Empty());

    p.PushBack(e1);
    EXPECT_EQ(cp.Conjugate(), p);

    p.PushBack(e2, Gap(100));
    EXPECT_EQ(cp.Conjugate(), p);

    p.PushBack(e3, Gap(10));
    EXPECT_EQ(cp.Conjugate(), p);

    p.PushBack(e4, Gap(10));
    EXPECT_EQ(cp.Conjugate(), p);
    EXPECT_EQ(cp.Front(), g.conjugate(e4));
    EXPECT_EQ(cp[1], g.conjugate(e3));
    EXPECT_EQ(cp[2], g.conjugate(e2));
    EXPECT_EQ(cp.Back(), g.conjugate(e1));
    EXPECT_EQ(cp.GapAt(0).gap, 0);
    EXPECT_EQ(cp.GapAt(1).gap, 10);
    EXPECT_EQ(cp.GapAt(2).gap, 10);
    EXPECT_EQ(cp.GapAt(3).gap, 100);
    EXPECT_EQ(cp.LengthAt(0), 1182);
    EXPECT_EQ(cp.LengthAt(1), 1116);
    EXPECT_EQ(cp.LengthAt(2), 527);
    EXPECT_EQ(cp.LengthAt(3), 426);

    p.PopBack(2);
    EXPECT_EQ(cp.Conjugate(), p);
    EXPECT_EQ(cp.Back(), g.conjugate(e1));
    EXPECT_EQ(cp.Front(), g.conjugate(e2));

    p.PopBack();
    EXPECT_EQ(cp.Conjugate(), p);
    EXPECT_EQ(cp.Back(), g.conjugate(e1));
    EXPECT_EQ(cp.Front(), g.conjugate(e1));

    p.PushBack(e2, Gap(100));
    EXPECT_EQ(cp.Conjugate(), p);
    EXPECT_EQ(cp.Back(), g.conjugate(e1));
    EXPECT_EQ(cp.Front(), g.conjugate(e2));

    p.PushBack(e3, Gap(10));
    EXPECT_EQ(cp.Conjugate(), p);
    EXPECT_EQ(cp.Back(), g.conjugate(e1));
    EXPECT_EQ(cp.Front(), g.conjugate(e3));

    p.PushBack(e4, Gap(10));
    EXPECT_EQ(cp.Conjugate(), p);
    EXPECT_EQ(cp.Front(), g.conjugate(e4));
    EXPECT_EQ(cp[1], g.conjugate(e3));
    EXPECT_EQ(cp[2], g.conjugate(e2));
    EXPECT_EQ(cp.Back(), g.conjugate(e1));
    EXPECT_EQ(cp.GapAt(0).gap, 0);
    EXPECT_EQ(cp.GapAt(1).gap, 10);
    EXPECT_EQ(cp.GapAt(2).gap, 10);
    EXPECT_EQ(cp.GapAt(3).gap, 100);
    EXPECT_EQ(cp.LengthAt(0), 1182);
    EXPECT_EQ(cp.LengthAt(1), 1116);
    EXPECT_EQ(cp.LengthAt(2), 527);
    EXPECT_EQ(cp.LengthAt(3), 426);
}


TEST( PathExtend, BidirectionalPathSearch ) {
    Graph g(13);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g));
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

    EXPECT_EQ(path1.FindFirst(e1), 0);
    EXPECT_EQ(path1.FindFirst(e5), 4);
    EXPECT_EQ(path1.FindFirst(e6), 5);
    EXPECT_EQ(path1.FindLast(e1), 0);
    EXPECT_EQ(path1.FindLast(e5), 10);
    EXPECT_EQ(path1.FindLast(e6), 9);

    auto v = path1.FindAll(e2);
    EXPECT_EQ(v.size(), 1);
    EXPECT_EQ(v[0], 1);
    v = path1.FindAll(e6);
    EXPECT_EQ(v.size(), 3);
    EXPECT_EQ(v[1], 7);
    v = path1.FindAll(e5);
    EXPECT_EQ(v.size(), 4);
    EXPECT_EQ(v[1], 6);
    EXPECT_EQ(v[2], 8);
    v = path1.FindAll(start);
    EXPECT_EQ(v.size(), 0);
    EXPECT_EQ(path1.FindFirst(start), -1);
    EXPECT_EQ(path1.FindLast(start), -1);

    EXPECT_TRUE(path1.CompareFrom(0, path1));
    EXPECT_FALSE(path1.CompareFrom(1, path1));
    EXPECT_FALSE(path1.CompareFrom(20, path1));
    EXPECT_EQ(path1.FindFirst(path1), 0);
    EXPECT_EQ(path1.FindLast(path1), 0);

    BidirectionalPath path2 = path1.SubPath(1, 4);
    BidirectionalPath path3 = path1.SubPath(4, 7);;
    BidirectionalPath path4 = path1.SubPath(9);
    BidirectionalPath path5(g);
    path5.PushBack(e4);
    path5.PushBack(e5);
    path5.PushBack(e6);
    path5.PushBack(e5);
    path5.PushBack(e7);

    EXPECT_FALSE(path1.CompareFrom(0, path2));
    EXPECT_TRUE(path1.CompareFrom(1, path2));
    EXPECT_FALSE(path1.CompareFrom(2, path2));
    EXPECT_EQ(path1.FindFirst(path2), 1);
    EXPECT_EQ(path1.FindLast(path2), 1);
    //EXPECT_EQ(path1.Contains(path2), true);

    EXPECT_TRUE(path1.CompareFrom(4, path3));
    EXPECT_FALSE(path1.CompareFrom(5, path3));
    EXPECT_TRUE(path1.CompareFrom(6, path3));
    EXPECT_EQ(path1.FindFirst(path3), 4);
    EXPECT_EQ(path1.FindLast(path3), 8);
    //EXPECT_EQ(path1.Contains(path3), true);

    EXPECT_TRUE(path1.CompareFrom(9, path4));
    EXPECT_FALSE(path1.CompareFrom(11, path4));
    EXPECT_EQ(path1.FindFirst(path4), 9);
    EXPECT_EQ(path1.FindLast(path4), 9);
    //EXPECT_EQ(path1.Contains(path4), true);

    EXPECT_FALSE(path1.CompareFrom(0, path5));
    EXPECT_FALSE(path1.CompareFrom(3, path5));
    EXPECT_FALSE(path1.CompareFrom(5, path5));
    EXPECT_EQ(path1.FindFirst(path5), -1);
    EXPECT_EQ(path1.FindLast(path5), -1);
    //EXPECT_EQ(path1.Contains(path5), false);
}


TEST( PathExtend, BidirectionalPathLoopDetector ) {
    Graph g(13);
    ASSERT_TRUE(graphio::ScanBasicGraph("./src/test/debruijn/graph_fragments/path_extend/distance_estimation", g));
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
    EXPECT_FALSE(loop_detect.IsCycled(4, skip_identical));
    EXPECT_TRUE(loop_detect.IsCycled(3, skip_identical));
    EXPECT_TRUE(loop_detect.IsCycled(2, skip_identical));
    loop_detect.RemoveLoop(skip_identical, false);
    EXPECT_EQ(path1.Size(), 6);
    EXPECT_EQ(path1.Back(), e6);

    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);

    skip_identical = 0;
    EXPECT_FALSE(loop_detect.IsCycled(4, skip_identical));
    skip_identical = 0;
    EXPECT_TRUE(loop_detect.IsCycled(3, skip_identical));
    skip_identical = 0;
    EXPECT_TRUE(loop_detect.IsCycled(2, skip_identical));

    loop_detect.RemoveLoop(skip_identical, true);
    EXPECT_EQ(path1.Size(), 5);
    EXPECT_EQ(path1.Back(), e5);

    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e6);
    path1.PushBack(e5);
    path1.PushBack(e7);

    skip_identical = 0;
    EXPECT_FALSE(loop_detect.IsCycled(4, skip_identical));
    skip_identical = 0;
    EXPECT_FALSE(loop_detect.IsCycled(3, skip_identical));
    skip_identical = 0;
    EXPECT_FALSE(loop_detect.IsCycled(2, skip_identical));

    EXPECT_EQ(loop_detect.EdgesToRemove(skip_identical, false), 0);
    loop_detect.RemoveLoop(skip_identical, false);
    EXPECT_EQ(path1.Size(), 12);
    EXPECT_EQ(path1.Back(), e7);
}
