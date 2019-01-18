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
    ASSERT_EQ('G', nucl(Sequence("ACGTACGTAC")[2]));
    ASSERT_EQ('A', nucl(Sequence("A")[0]));
}

TEST( Sequence, ZeroLength ) {
    Sequence s("");
    ASSERT_EQ(0, s.size());
}

TEST( Sequence, NullValue ) {
        Sequence s("");
    ASSERT_EQ("", (!s).str());
}

TEST( Sequence, Sum ) {
    ASSERT_EQ("ACG", (Sequence("A") + Sequence("CG")).str());
    ASSERT_EQ("ACGTTGCA", (Sequence("ACGT") + Sequence("TGCA")).str());
    ASSERT_EQ("ACGTACGTTGCATGCA", (Sequence("ACGTACGT") + Sequence("TGCATGCA")).str());
}

TEST( Sequence, Str ) {
    ASSERT_EQ("ACGTACGTAC", Sequence("ACGTACGTAC").str());
    ASSERT_EQ("ACG", Sequence("ACG").str());
}

TEST( Sequence, ReverseComplement ) {
    Sequence s = Sequence("AACCGGTTAA");
    ASSERT_EQ("TTAACCGGTT", (!s).str());
    Sequence s2 = Sequence("ACG");
    ASSERT_EQ("CGT", (!s2).str());
}
