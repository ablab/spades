#pragma once
#include <boost/test/unit_test.hpp>
#include "sequence/rtseq.hpp"
#include "sequence/sequence.hpp"
#include "adt/cyclichash.hpp"
#include "utils/verify.hpp"

BOOST_AUTO_TEST_CASE( TestBasic ) {
    unsigned k = 4;
    typedef uint64_t digest;
    rolling_hash::SymmetricCyclicHash<rolling_hash::DNASeqHash> hasher(k);
    auto h1 = hasher(std::string("ACCG"));
    auto h2 = hasher.hash(std::string("ACCG"));
    auto h3 = hasher.hash(std::string("CGGT"));
    auto h4 = hasher.hash(std::string("GGTA"));

    BOOST_CHECK_EQUAL((digest) h1, (digest) h1);
    BOOST_CHECK_EQUAL((digest) h1, h1.value());
    BOOST_CHECK_EQUAL((digest) h1, (digest) h2);
    BOOST_CHECK_EQUAL((digest) h1, (digest) h3);

    BOOST_CHECK((digest) h3 != (digest) h4);

    BOOST_CHECK_EQUAL((digest) hasher.hash_update(h3, 'C', 'A'), (digest) h4);
}

BOOST_AUTO_TEST_CASE( TestBasicRtSeq ) {
    unsigned k = 4;
    typedef uint64_t digest;
    rolling_hash::SymmetricCyclicHash<rolling_hash::DNASeqHash> hasher(k);
    auto h1 = hasher(RtSeq(4, "ACCG"));
    auto h2 = hasher.hash(RtSeq(4, "ACCG"));
    auto h3 = hasher.hash(RtSeq(4, "CGGT"));
    auto h4 = hasher.hash(RtSeq(4, "GGTA"));

    BOOST_CHECK_EQUAL((digest) h1, (digest) h1);
    BOOST_CHECK_EQUAL((digest) h1, h1.value());
    BOOST_CHECK_EQUAL((digest) h1, (digest) h2);
    BOOST_CHECK_EQUAL((digest) h1, (digest) h3);
    BOOST_CHECK((digest) h3 != (digest) h4);

    //BOOST_CHECK_EQUAL((digest) hasher.hash_update(h3, 'C', 'A'), (digest) h4);
    BOOST_CHECK_EQUAL((digest) hasher.hash_update(h3, 'C', 'A'), (digest) h4);

}

BOOST_AUTO_TEST_CASE( TestRC ) {
    unsigned k = 4;
    typedef uint64_t digest;
    rolling_hash::SymmetricCyclicHash<rolling_hash::DNASeqHash> hasher(k);

    Sequence s("AACCTTGGACGTCGTAACGACT");
    Sequence s2 = !s;
    size_t kmer_cnt = s.size() - k + 1;
    for (size_t i = 0; i < kmer_cnt; ++i) {
        VERIFY((digest) hasher(RtSeq(k, s, i)) == (digest) hasher(RtSeq(k, s2, kmer_cnt - i - 1)));
        BOOST_CHECK_EQUAL((digest) hasher(RtSeq(k, s, i)), (digest) hasher(RtSeq(k, s2, kmer_cnt - i - 1)));
    }
}

BOOST_AUTO_TEST_CASE( TestRoll ) {
    unsigned k = 4;
    typedef uint64_t digest;
    rolling_hash::SymmetricCyclicHash<rolling_hash::NDNASeqHash> hasher(k);

    Sequence s("AACCTTGGACGTCGTAACGACT");
    size_t kmer_cnt = s.size() - k + 1;
    auto hash = hasher(s);
    for (size_t i = 1; i < kmer_cnt; ++i) {
        hash = hasher.hash_update(hash, s[i - 1], s[i - 1 + k]);
        VERIFY((digest) hasher(RtSeq(k, s, i)) == (digest) hash);
        BOOST_CHECK_EQUAL((digest) hasher(RtSeq(k, s, i)), (digest) hash);
    }
}
