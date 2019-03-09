//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "sequence/rtseq.hpp"
#include "sequence/sequence.hpp"
#include "adt/cyclichash.hpp"
#include <gtest/gtest.h>

TEST( CyclicHash, TestBasic ) {
    unsigned k = 4;
    typedef uint64_t digest;
    rolling_hash::SymmetricCyclicHash<rolling_hash::DNASeqHash> hasher(k);
    auto h1 = hasher(std::string("ACCG"));
    auto h2 = hasher.hash(std::string("ACCG"));
    auto h3 = hasher.hash(std::string("CGGT"));
    auto h4 = hasher.hash(std::string("GGTA"));

    EXPECT_EQ((digest) h1, (digest) h1);
    EXPECT_EQ((digest) h1, h1.value());
    EXPECT_EQ((digest) h1, (digest) h2);
    EXPECT_EQ((digest) h1, (digest) h3);

    ASSERT_NE((digest) h3, (digest) h4);

    EXPECT_EQ((digest) hasher.hash_update(h3, 'C', 'A'), (digest) h4);
}

TEST( CyclicHash, TestBasicRtSeq ) {
    unsigned k = 4;
    typedef uint64_t digest;
    rolling_hash::SymmetricCyclicHash<rolling_hash::DNASeqHash> hasher(k);
    auto h1 = hasher(RtSeq(4, "ACCG"));
    auto h2 = hasher.hash(RtSeq(4, "ACCG"));
    auto h3 = hasher.hash(RtSeq(4, "CGGT"));
    auto h4 = hasher.hash(RtSeq(4, "GGTA"));

    EXPECT_EQ((digest) h1, (digest) h1);
    EXPECT_EQ((digest) h1, h1.value());
    EXPECT_EQ((digest) h1, (digest) h2);
    EXPECT_EQ((digest) h1, (digest) h3);
    ASSERT_NE((digest) h3, (digest) h4);

    //EXPECT_EQ((digest) hasher.hash_update(h3, 'C', 'A'), (digest) h4);
    EXPECT_EQ((digest) hasher.hash_update(h3, 'C', 'A'), (digest) h4);

}

TEST( CyclicHash, TestRC ) {
    unsigned k = 4;
    typedef uint64_t digest;
    rolling_hash::SymmetricCyclicHash<rolling_hash::DNASeqHash> hasher(k);

    Sequence s("AACCTTGGACGTCGTAACGACT");
    Sequence s2 = !s;
    size_t kmer_cnt = s.size() - k + 1;
    for (size_t i = 0; i < kmer_cnt; ++i) {
        EXPECT_EQ((digest) hasher(RtSeq(k, s, i)), (digest) hasher(RtSeq(k, s2, kmer_cnt - i - 1)));
    }
}

TEST( CyclicHash, TestRoll ) {
    unsigned k = 4;
    typedef uint64_t digest;
    rolling_hash::SymmetricCyclicHash<rolling_hash::NDNASeqHash> hasher(k);

    Sequence s("AACCTTGGACGTCGTAACGACT");
    size_t kmer_cnt = s.size() - k + 1;
    auto hash = hasher(s);
    for (size_t i = 1; i < kmer_cnt; ++i) {
        hash = hasher.hash_update(hash, s[i - 1], s[i - 1 + k]);
        EXPECT_EQ((digest) hasher(RtSeq(k, s, i)), (digest) hash);
    }
}
