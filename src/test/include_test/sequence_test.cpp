//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "sequence/sequence.hpp"
#include "sequence/nucl.hpp"
#include <string>
#include <gtest/gtest.h>

TEST( Sequence, Selector ) {
    Sequence s("TTATTAGGGAT");
    EXPECT_EQ('G', nucl(Sequence("ACGTACGTAC")[2]));
    EXPECT_EQ('A', nucl(Sequence("A")[0]));
}

TEST( Sequence, ZeroLength ) {
    Sequence s("");
    EXPECT_EQ(0, s.size());
}

TEST( Sequence, NullValue ) {
        Sequence s("");
    EXPECT_EQ("", (!s).str());
}

TEST( Sequence, Sum ) {
    EXPECT_EQ("ACG", (Sequence("A") + Sequence("CG")).str());
    EXPECT_EQ("ACGTTGCA", (Sequence("ACGT") + Sequence("TGCA")).str());
    EXPECT_EQ("ACGTACGTTGCATGCA", (Sequence("ACGTACGT") + Sequence("TGCATGCA")).str());
}

TEST( Sequence, Str ) {
    EXPECT_EQ("ACGTACGTAC", Sequence("ACGTACGTAC").str());
    EXPECT_EQ("ACG", Sequence("ACG").str());
}

TEST( Sequence, ReverseComplement ) {
    Sequence s = Sequence("AACCGGTTAA");
    EXPECT_EQ("TTAACCGGTT", (!s).str());
    Sequence s2 = Sequence("ACG");
    EXPECT_EQ("CGT", (!s2).str());
}
