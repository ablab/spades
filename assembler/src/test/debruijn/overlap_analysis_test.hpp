#pragma once

#include <boost/test/unit_test.hpp>
#include "overlap_analysis.hpp"
//#include "repeat_resolving_routine.hpp"

namespace debruijn_graph {

BOOST_AUTO_TEST_SUITE(overlap_analysis_tests)

BOOST_AUTO_TEST_CASE( TrivialTest ) {
    SWOverlapAnalyzer analyzer(-1ul);

    Sequence s1("ACGTACGT");
    Sequence s2("ACGTACGT");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    std::cout << overlap << std::endl;
    BOOST_CHECK_NE(overlap.match_cnt, 0);
    BOOST_CHECK_EQUAL(overlap.match_cnt, 8);
}

BOOST_AUTO_TEST_CASE( MismatchTest ) {
    SWOverlapAnalyzer analyzer(-1ul);

    Sequence s1("GTACTTACG");
    Sequence s2("GTACGTACG");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    std::cout << overlap << std::endl;
    BOOST_CHECK_EQUAL(overlap.match_cnt, 8);
}

BOOST_AUTO_TEST_CASE( MismatchTest2 ) {
    SWOverlapAnalyzer analyzer(-1ul);
    //two nucleotide match is not enough to compensate mismatch
    Sequence s1("GTACTTA");
    Sequence s2("GTACGTA");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    std::cout << overlap << std::endl;
    BOOST_CHECK_EQUAL(overlap.match_cnt, 4);
}

BOOST_AUTO_TEST_CASE( DelTest ) {
    SWOverlapAnalyzer analyzer(-1ul);

    Sequence s1("ACGTACTACGT");
    Sequence s2("ACGTACGTACGT");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    std::cout << overlap << std::endl;
    BOOST_CHECK_EQUAL(overlap.match_cnt, 11);
}

BOOST_AUTO_TEST_CASE( InsTest ) {
    SWOverlapAnalyzer analyzer(-1ul);

    Sequence s1("ACGTACGTACGT");
    Sequence s2("ACGTACTACGT");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    std::cout << overlap << std::endl;
    BOOST_CHECK_EQUAL(overlap.match_cnt, 11);
}

BOOST_AUTO_TEST_CASE( QueryLongerTest ) {
    SWOverlapAnalyzer analyzer(-1ul);

    Sequence s1("AACGTACGTT");
    Sequence s2("ACGTACGT");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    std::cout << overlap << std::endl;
    BOOST_CHECK_NE(overlap.match_cnt, 0);
    BOOST_CHECK_EQUAL(overlap.match_cnt, 8);
}

BOOST_AUTO_TEST_CASE( RefLongerTest ) {
    SWOverlapAnalyzer analyzer(-1ul);

    Sequence s1("ACGTACGT");
    Sequence s2("AACGTACGTT");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    std::cout << overlap << std::endl;
    BOOST_CHECK_NE(overlap.match_cnt, 0);
    BOOST_CHECK_EQUAL(overlap.match_cnt, 8);
}

BOOST_AUTO_TEST_CASE( TrivialBoundedFlankTest ) {
    SWOverlapAnalyzer analyzer(4);

    Sequence s1("CCCCACGT");
    Sequence s2("ACGTAAAA");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    std::cout << overlap << std::endl;
    BOOST_CHECK_NE(overlap.match_cnt, 0);
    BOOST_CHECK_EQUAL(overlap.match_cnt, 4);
}

BOOST_AUTO_TEST_CASE( PartialPerfectAlignmentTest ) {
    SWOverlapAnalyzer analyzer(6);

    Sequence s1("ACGTACGT");
    Sequence s2("ACGTACGT");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    std::cout << overlap << std::endl;
    BOOST_CHECK_NE(overlap.match_cnt, 0);
    BOOST_CHECK_EQUAL(overlap.match_cnt, 4);
}

BOOST_AUTO_TEST_CASE( UnalignableTest ) {
    SWOverlapAnalyzer analyzer(6);

    Sequence s1("CCCCCCCC");
    Sequence s2("AAAAAAAA");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    std::cout << overlap << std::endl;
    BOOST_CHECK_EQUAL(overlap, OverlapInfo());
}

BOOST_AUTO_TEST_CASE( LargeTest ) {
    SWOverlapAnalyzer analyzer(25);

    Sequence s1("CCCCCCCCACGTTCGTACATACGTGGGGGGG");
    Sequence s2("AAAAAAAAACGTACGTACGTACGTCCCCCCC");
    OverlapInfo overlap = analyzer.AnalyzeOverlap(s1, s2);
    std::cout << overlap << std::endl;
    BOOST_CHECK_EQUAL(overlap.match_cnt, 14);
    BOOST_CHECK_EQUAL(overlap.r1, Range(8, 24));
    BOOST_CHECK_EQUAL(overlap.r1, Range(8, 24));
}

BOOST_AUTO_TEST_SUITE_END()

}
