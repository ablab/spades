//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "sequence/seq.hpp"
#include "sequence/sequence.hpp"
#include "sequence/nucl.hpp"

#include <string>
#include <gtest/gtest.h>

typedef unsigned long long ull;

TEST( Seq, Selector) {
    ASSERT_EQ('G', nucl(Seq<10>("ACGTACGTAC")[2]));
    ASSERT_EQ('G', nucl(Seq<60,ull>("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[2]));
    ASSERT_EQ('G', nucl(Seq<60,ull>("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[16]));
    ASSERT_EQ('T', nucl(Seq<60,ull>("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[17]));
    ASSERT_EQ('A', nucl(Seq<60,ull>("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[18]));
    ASSERT_EQ('C', nucl(Seq<60,ull>("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")[19]));
}

TEST( Seq, ShiftLeft ) {
    Seq<10> s("ACGTACGTAC");
    ASSERT_EQ("CGTACGTACA", (s << dignucl('A')).str());
    ASSERT_EQ("CGTACGTACC", (s << dignucl('C')).str());
    ASSERT_EQ("CGTACGTACG", (s << dignucl('G')).str());
    ASSERT_EQ("CGTACGTACT", (s << dignucl('T')).str());
    Seq<60,ull> s2("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    ASSERT_EQ("CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACA", (s2 << dignucl('A')).str());
    ASSERT_EQ("CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACC", (s2 << dignucl('C')).str());
    ASSERT_EQ("CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACG", (s2 << dignucl('G')).str());
    ASSERT_EQ("CGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACT", (s2 << dignucl('T')).str());
}

TEST( Seq, ShiftRight ) {
    Seq<10> s("ACGTACGTAC");
    ASSERT_EQ("AACGTACGTA", (s >> dignucl('A')).str());
    ASSERT_EQ("CACGTACGTA", (s >> dignucl('C')).str());
    ASSERT_EQ("GACGTACGTA", (s >> dignucl('G')).str());
    ASSERT_EQ("TACGTACGTA", (s >> dignucl('T')).str());
    Seq<60,ull> s2("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    ASSERT_EQ("AACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA", (s2 >> dignucl('A')).str());
    ASSERT_EQ("CACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA", (s2 >> dignucl('C')).str());
    ASSERT_EQ("GACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA", (s2 >> dignucl('G')).str());
    ASSERT_EQ("TACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA", (s2 >> dignucl('T')).str());
    ASSERT_TRUE(Seq<5>("ACACA") == Seq<5>("CACAC")>>dignucl('A'));
}

TEST( Seq, Str ) {
    Seq<10> s("ACGTACGTAC");
    ASSERT_EQ("ACGTACGTAC", s.str());
    Seq<60,ull> s2("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    ASSERT_EQ("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC", s2.str());
}

TEST( Seq, HeadAndTail ) {
    Seq<10> s("ACGTACGTAC");
    ASSERT_EQ("CGTACGTAC", Seq<9>(s, 1).str()); // tail
    ASSERT_EQ("ACGTACGTA", Seq<9>(s).str()); // head
}

TEST( Seq, FromBiggerSeq ) {
    Seq<10> s("ACGTACGTAC");
    ASSERT_EQ("ACGTA", Seq<5>(s).str());
}

TEST( Seq, FromType ) {
    Sequence s("ACGTACGTAC");
    ASSERT_EQ("ACGTA", Seq<5>(s).str());
    ASSERT_EQ("GTACG", Seq<5>(s, 2).str());
}

TEST( Seq, PushBack ) {
    {
        Seq<4> s("ACGT");
        ASSERT_EQ("ACGTC", s.pushBack('C').str());
    }
    {
        Seq<4, unsigned char> s("ACGT");
        ASSERT_EQ("ACGTC", s.pushBack('C').str());
    }
    {
        Seq<3> s("ACG");
        ASSERT_EQ("ACGC", s.pushBack('C').str());
    }
}

TEST( Seq, Null ) {
        Seq<0> s("");
    ASSERT_EQ("", s.str());
}

TEST( Seq, EndNull ) {
        Seq<0> s("");
    ASSERT_EQ("", s.end<0>().str());
}

TEST( Seq, StartNull ) {
        Seq<0> s("");
    ASSERT_EQ("", s.start<0>().str());
}

TEST( Seq, AddSymbolForNullValue ) {
    Seq<1> s1("G");
    Seq<1> s2 = (s1 << 'A');
    Seq<1> s3("A");
    ASSERT_EQ(s3.str(), s2.str());
}

TEST( Seq, End ) {
    Seq<5> s1("ACGTA");
    ASSERT_EQ("CGTA", s1.end<4>().str());
}

TEST( Seq, Start ) {
    Seq<5> s1("ACGTA");
    ASSERT_EQ("ACGT", s1.start<4>().str());
}

TEST( Seq, Complex ) {
    Sequence s1("ACAAA");
    Sequence s2("CAAAC");
    ASSERT_EQ((!(Seq<4>(!s1))).str(), Seq<4>(s2).str());
    ASSERT_EQ(!(Seq<4>(!s1)), Seq<4>(s2));
}

TEST( Seq, FromCharArray ) {
    std::string s = "ACGTACGTAC";
    ASSERT_EQ("ACGTACGTAC", Seq<10>(s.c_str()).str());
}

TEST( Seq, ReverseComplement ) {
    Seq<10> s("ACGTACGTAC");
    ASSERT_EQ("GTACGTACGT", (!s).str());
    Seq<60,ull> s2("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC");
    ASSERT_EQ("GTACGTACGTGTACGTACGTGTACGTACGTGTACGTACGTGTACGTACGTGTACGTACGT", (!s2).str());
    Seq<9> s3("ACGTACGTA");
    ASSERT_EQ("TACGTACGT", (!s3).str());
}

TEST( Seq, Test16 ) {
    Seq<16> s("AAAAAAAAAAAAAAAA");
    ASSERT_EQ(s << 'C', Seq<16>("AAAAAAAAAAAAAAAC"));
}

TEST( Seq, Test16_2 ) {
    Seq<16> s("TTTTTTTTTTTTTTTT");
    ASSERT_EQ(Seq<16>("TTTTTTTTTTTTTTTA"), s << 'A');
}

TEST( Seq, TestFirstLast ) {
    Seq<7> s1("ACGTACT");
    ASSERT_EQ(0, s1.first());
    ASSERT_EQ(3, s1.last());
    Seq<7> s2("TTTTTTT");
    ASSERT_EQ(3, s2.first());
    ASSERT_EQ(3, s2.last());
}
