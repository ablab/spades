//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include <boost/test/unit_test.hpp>
#include "sequence/seq.hpp"
#include "sequence/sequence.hpp"
#include "sequence/nucl.hpp"
#include <string>

typedef unsigned long long ull;

BOOST_AUTO_TEST_CASE( TestSeqSelector ) {
    BOOST_CHECK_EQUAL('G', nucl(Seq<10>("ACGTACGTAC")[2]));
    BOOST_CHECK_EQUAL('G', nucl(Seq<60,ull>("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[2]));
    BOOST_CHECK_EQUAL('G', nucl(Seq<60,ull>("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[16]));
    BOOST_CHECK_EQUAL('T', nucl(Seq<60,ull>("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[17]));
    BOOST_CHECK_EQUAL('A', nucl(Seq<60,ull>("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[18]));
    BOOST_CHECK_EQUAL('C', nucl(Seq<60,ull>("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[19]));
}

BOOST_AUTO_TEST_CASE( TestSeqShiftLeft ) {
    Seq<10> s("ACGTACGTAC");
    BOOST_CHECK_EQUAL("CGTACGTACA", (s << dignucl('A')).str());
    BOOST_CHECK_EQUAL("CGTACGTACC", (s << dignucl('C')).str());
    BOOST_CHECK_EQUAL("CGTACGTACG", (s << dignucl('G')).str());
    BOOST_CHECK_EQUAL("CGTACGTACT", (s << dignucl('T')).str());
    Seq<60,ull> s2("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    BOOST_CHECK_EQUAL("CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACA", (s2 << dignucl('A')).str());
    BOOST_CHECK_EQUAL("CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACC", (s2 << dignucl('C')).str());
    BOOST_CHECK_EQUAL("CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACG", (s2 << dignucl('G')).str());
    BOOST_CHECK_EQUAL("CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACT", (s2 << dignucl('T')).str());
}

BOOST_AUTO_TEST_CASE( TestSeqShiftRight ) {
    Seq<10> s("ACGTACGTAC");
    BOOST_CHECK_EQUAL("AACGTACGTA", (s >> dignucl('A')).str());
    BOOST_CHECK_EQUAL("CACGTACGTA", (s >> dignucl('C')).str());
    BOOST_CHECK_EQUAL("GACGTACGTA", (s >> dignucl('G')).str());
    BOOST_CHECK_EQUAL("TACGTACGTA", (s >> dignucl('T')).str());
    Seq<60,ull> s2("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    BOOST_CHECK_EQUAL("AACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA", (s2 >> dignucl('A')).str());
    BOOST_CHECK_EQUAL("CACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA", (s2 >> dignucl('C')).str());
    BOOST_CHECK_EQUAL("GACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA", (s2 >> dignucl('G')).str());
    BOOST_CHECK_EQUAL("TACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA", (s2 >> dignucl('T')).str());
    BOOST_CHECK(Seq<5>("ACACA") == Seq<5>("CACAC")>>dignucl('A'));
}

BOOST_AUTO_TEST_CASE( TestSeqStr ) {
    Seq<10> s("ACGTACGTAC");
    BOOST_CHECK_EQUAL("ACGTACGTAC", s.str());
    Seq<60,ull> s2("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    BOOST_CHECK_EQUAL("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", s2.str());
}

BOOST_AUTO_TEST_CASE( TestSeqHeadAndTail ) {
    Seq<10> s("ACGTACGTAC");
    BOOST_CHECK_EQUAL("CGTACGTAC", Seq<9>(s, 1).str()); // tail
    BOOST_CHECK_EQUAL("ACGTACGTA", Seq<9>(s).str()); // head
}

BOOST_AUTO_TEST_CASE( TestSeqFromBiggerSeq ) {
    Seq<10> s("ACGTACGTAC");
    BOOST_CHECK_EQUAL("ACGTA", Seq<5>(s).str());
}

BOOST_AUTO_TEST_CASE( TestSeqFromType ) {
    Sequence s("ACGTACGTAC");
    BOOST_CHECK_EQUAL("ACGTA", Seq<5>(s).str());
    BOOST_CHECK_EQUAL("GTACG", Seq<5>(s, 2).str());
}

BOOST_AUTO_TEST_CASE( TestSeqPushBack ) {
    {
        Seq<4> s("ACGT");
        BOOST_CHECK_EQUAL("ACGTC", s.pushBack('C').str());
    }
    {
        Seq<4, unsigned char> s("ACGT");
        BOOST_CHECK_EQUAL("ACGTC", s.pushBack('C').str());
    }
    {
        Seq<3> s("ACG");
        BOOST_CHECK_EQUAL("ACGC", s.pushBack('C').str());
    }
}

BOOST_AUTO_TEST_CASE( TestSeqNull ) {
        Seq<0> s("");
    BOOST_CHECK_EQUAL("", s.str());
}

BOOST_AUTO_TEST_CASE( TestSeqEndNull ) {
        Seq<0> s("");
    BOOST_CHECK_EQUAL("", s.end<0>().str());
}

BOOST_AUTO_TEST_CASE( TestSeqStartNull ) {
        Seq<0> s("");
    BOOST_CHECK_EQUAL("", s.start<0>().str());
}

BOOST_AUTO_TEST_CASE( TestSeqAddSymbolForNullValue ) {
    Seq<1> s1("G");
    Seq<1> s2 = (s1 << 'A');
    Seq<1> s3("A");
    BOOST_CHECK_EQUAL(s3.str(), s2.str());
}

BOOST_AUTO_TEST_CASE( TestSeqEnd ) {
    Seq<5> s1("ACGTA");
    BOOST_CHECK_EQUAL("CGTA", s1.end<4>().str());
}

BOOST_AUTO_TEST_CASE( TestSeqStart ) {
    Seq<5> s1("ACGTA");
    BOOST_CHECK_EQUAL("ACGT", s1.start<4>().str());
}

BOOST_AUTO_TEST_CASE( TestSeqComplex ) {
    Sequence s1("ACAAA");
    Sequence s2("CAAAC");
    BOOST_CHECK_EQUAL((!(Seq<4>(!s1))).str(), Seq<4>(s2).str());
    BOOST_CHECK_EQUAL(!(Seq<4>(!s1)), Seq<4>(s2));
}

BOOST_AUTO_TEST_CASE( TestSeqFromCharArray ) {
    std::string s = "ACGTACGTAC";
    BOOST_CHECK_EQUAL("ACGTACGTAC", Seq<10>(s.c_str()).str());
}

BOOST_AUTO_TEST_CASE( TestSeqReverseComplement ) {
    Seq<10> s("ACGTACGTAC");
    BOOST_CHECK_EQUAL("GTACGTACGT", (!s).str());
    Seq<60,ull> s2("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    BOOST_CHECK_EQUAL("GTACGTACGTGTACGTACGTGTACGTACGTGTACGTACGTGTACGTACGTGTACGTACGT", (!s2).str());
    Seq<9> s3("ACGTACGTA");
    BOOST_CHECK_EQUAL("TACGTACGT", (!s3).str());
}

BOOST_AUTO_TEST_CASE( Test16 ) {
    Seq<16> s("AAAAAAAAAAAAAAAA");
    BOOST_CHECK_EQUAL(s << 'C', Seq<16>("AAAAAAAAAAAAAAAC"));
}

BOOST_AUTO_TEST_CASE( Test16_2 ) {
    Seq<16> s("TTTTTTTTTTTTTTTT");
    BOOST_CHECK_EQUAL(Seq<16>("TTTTTTTTTTTTTTTA"), s << 'A');
}

BOOST_AUTO_TEST_CASE( TestFirstLast ) {
    Seq<7> s1("ACGTACT");
    BOOST_CHECK_EQUAL(0, s1.first());
    BOOST_CHECK_EQUAL(3, s1.last());
    Seq<7> s2("TTTTTTT");
    BOOST_CHECK_EQUAL(3, s2.first());
    BOOST_CHECK_EQUAL(3, s2.last());
}
