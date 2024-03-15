//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "modules/path_extend/overlap_analysis.hpp"
#include "modules/path_extend/overlap_remover.hpp"
#include "sequence/sequence.hpp"

#include <gtest/gtest.h>

using namespace debruijn_graph;

TEST( OverlapAnalysis, Trivial ) {
    SWOverlapAnalyzer analyzer(-1ul);

    Sequence s1("ACGTACGT");
    Sequence s2("ACGTACGT");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    EXPECT_NE(overlap.match_cnt, 0) << overlap;
    EXPECT_EQ(overlap.match_cnt, 8) << overlap;
}

TEST( OverlapAnalysis, Mismatch ) {
    SWOverlapAnalyzer analyzer(-1ul);

    Sequence s1("GTACTTACG");
    Sequence s2("GTACGTACG");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    EXPECT_EQ(overlap.match_cnt, 8) << overlap;
}

TEST( OverlapAnalysis, Mismatch2 ) {
    SWOverlapAnalyzer analyzer(-1ul);
    //two nucleotide match is not enough to compensate mismatch
    Sequence s1("GTACTTA");
    Sequence s2("GTACGTA");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    EXPECT_EQ(overlap.match_cnt, 4) << overlap;
}

TEST( OverlapAnalysis, Deletion ) {
    SWOverlapAnalyzer analyzer(-1ul);

    Sequence s1("ACGTACTACGT");
    Sequence s2("ACGTACGTACGT");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    EXPECT_EQ(overlap.match_cnt, 11) << overlap;
}

TEST( OverlapAnalysis, Insertion ) {
    SWOverlapAnalyzer analyzer(-1ul);

    Sequence s1("ACGTACGTACGT");
    Sequence s2("ACGTACTACGT");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    EXPECT_EQ(overlap.match_cnt, 11) << overlap;
}

TEST( OverlapAnalysis, QueryLonger ) {
    SWOverlapAnalyzer analyzer(-1ul);

    Sequence s1("AACGTACGTT");
    Sequence s2("ACGTACGT");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    EXPECT_NE(overlap.match_cnt, 0) << overlap;
    EXPECT_EQ(overlap.match_cnt, 8) << overlap;
}

TEST( OverlapAnalysis, RefLonger ) {
    SWOverlapAnalyzer analyzer(-1ul);

    Sequence s1("ACGTACGT");
    Sequence s2("AACGTACGTT");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    EXPECT_NE(overlap.match_cnt, 0) << overlap;
    EXPECT_EQ(overlap.match_cnt, 8) << overlap;
}

TEST( OverlapAnalysis, TrivialBoundedFlank ) {
    SWOverlapAnalyzer analyzer(4);

    Sequence s1("CCCCACGT");
    Sequence s2("ACGTAAAA");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    EXPECT_NE(overlap.match_cnt, 0) << overlap;
    EXPECT_EQ(overlap.match_cnt, 4) << overlap;
}

TEST( OverlapAnalysis, PartialPerfectAlignment ) {
    SWOverlapAnalyzer analyzer(6);

    Sequence s1("ACGTACGT");
    Sequence s2("ACGTACGT");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    EXPECT_NE(overlap.match_cnt, 0) << overlap;
    EXPECT_EQ(overlap.match_cnt, 4) << overlap;
}

TEST( OverlapAnalysis, Unalignable ) {
    SWOverlapAnalyzer analyzer(6);

    Sequence s1("CCCCCCCC");
    Sequence s2("AAAAAAAA");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    EXPECT_EQ(overlap, OverlapInfo()) << overlap;
}

TEST( OverlapAnalysis, Large ) {
    SWOverlapAnalyzer analyzer(25);

    Sequence s1("CCCCCCCCACGTTCGTACATACGTGGGGGGG");
    Sequence s2("AAAAAAAAACGTACGTACGTACGTCCCCCCC");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    EXPECT_EQ(overlap.match_cnt, 14) << overlap;
    EXPECT_EQ(overlap.r1, Range(8, 24)) << overlap;
    EXPECT_EQ(overlap.r2, Range(8, 24)) << overlap;
}

inline path_extend::Gap MimicLAGapAnalyzer(size_t k, Sequence &s1, Sequence &s2) {
    const int INVALID_GAP = -1000000;
    constexpr static double IDENTITY_RATIO = 0.9;

    SWOverlapAnalyzer overlap_analyzer_(10000);
    auto overlap_info = overlap_analyzer_.AnalyzeOverlap(s1, s2);
    size_t min_la_length_ = 4;
    if (overlap_info.size() < min_la_length_) {
        DEBUG("Low alignment size");
        return path_extend::Gap(INVALID_GAP);
    }
    if (overlap_info.identity() < IDENTITY_RATIO) {
        DEBUG("Low identity score");
        return path_extend::Gap(INVALID_GAP);
    }
    // std::cout << overlap_info;

    return path_extend::Gap(
            (int) (k - overlap_info.r1.size() - s1.size() + overlap_info.r1.end_pos - overlap_info.r2.start_pos),
            {(uint32_t) (s1.size() - overlap_info.r1.end_pos), (uint32_t) overlap_info.r2.start_pos});
}

//TODO what does it test?! Why is it not a test of SWOverlapAnalyzer?
TEST( OverlapAnalysis, SimpleGap ) {
    Sequence s1("AAAAAAAACGCGCTTTCGCTTTAA");
    Sequence s2("GGGGCGCGCTTTCGCTAAAAAAAAAA");
    size_t k = 5;
    path_extend::Gap g = MimicLAGapAnalyzer(k, s1, s2);
    EXPECT_EQ(14, 14);
    EXPECT_EQ(g.trash.current, 4);
    EXPECT_EQ(g.trash.previous, 4);
    EXPECT_EQ(g.gap, -15);
}

TEST( OverlapAnalysis, UndefinedBehavior ) {
    SWOverlapAnalyzer analyzer(166);
    Sequence s1("ACGCAAGTAAGTGACGAAGGACATCCTCCCGCCCTCCCTTCCTCCCTGTCTTCATTCGCCTCCCTTCCCCGGTCTTCGCATTTCTGCAAGCGCTTTACCGAGCGGTCAGCGTGCGATAGACTCGCGCCGATCGCTTCTGCTGGCCTCACGCAGGCTGGGGGTGTTT");
    Sequence s2("CTCCCTCCTCCCGTCCTCCCCCTCCCTGTCTTCATTCGCCTCCCTTCCCCGGTCTTCGCATTTCTGCAAGCGCTTTACCGAGCGGTCAGCGTGCGATAGACTCGCGCCGATCGCTTCTGCTGGCCTCACGCAGGCTGGGGGTGTTTGTGTTTCGTCTGGACGCCGA");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    EXPECT_EQ(overlap.match_cnt, 140) << overlap;
    EXPECT_EQ(overlap.r1, Range(23, 166)) << overlap;
    EXPECT_EQ(overlap.r2, Range(5, 146)) << overlap;
}
