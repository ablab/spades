//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include <boost/test/unit_test.hpp>
#include "sequence/sequence.hpp"
#include "sequence/nucl.hpp"
#include <string>
#include "utils/perf/memory.hpp"
#include <ctime>

BOOST_AUTO_TEST_CASE( TestSequenceSelector ) {
    Sequence s("TTATTAGGGAT");
    BOOST_CHECK_EQUAL('G', nucl(Sequence("ACGTACGTAC")[2]));
    BOOST_CHECK_EQUAL('A', nucl(Sequence("A")[0]));
}

BOOST_AUTO_TEST_CASE( TestZeroLengthSequence ) {
    Sequence s("");
    BOOST_CHECK_EQUAL(0, s.size());
}

BOOST_AUTO_TEST_CASE( TestSequenceNullValue ) {
        Sequence s("");
    BOOST_CHECK_EQUAL("", (!s).str());
}

BOOST_AUTO_TEST_CASE( TestSequenceSum ) {
    BOOST_CHECK_EQUAL("ACG", (Sequence("A") + Sequence("CG")).str());
    BOOST_CHECK_EQUAL("ACGTTGCA", (Sequence("ACGT") + Sequence("TGCA")).str());
    BOOST_CHECK_EQUAL("ACGTACGTTGCATGCA", (Sequence("ACGTACGT") + Sequence("TGCATGCA")).str());
}

BOOST_AUTO_TEST_CASE( TestSequenceStr ) {
    BOOST_CHECK_EQUAL("ACGTACGTAC", Sequence("ACGTACGTAC").str());
    BOOST_CHECK_EQUAL("ACG", Sequence("ACG").str());
}

BOOST_AUTO_TEST_CASE( TestSequenceReverseComplement ) {
    Sequence s = Sequence("AACCGGTTAA");
    BOOST_CHECK_EQUAL("TTAACCGGTT", (!s).str());
    Sequence s2 = Sequence("ACG");
    BOOST_CHECK_EQUAL("CGT", (!s2).str());
}

//todo strange test
BOOST_AUTO_TEST_CASE( TestSequenceRefCount ) {
    Sequence s("AAAAAAA");
    Sequence s2(s);
    Sequence s3 = !s;
    Sequence s4 = s;
    Sequence ss = s.Subseq(3);
    BOOST_CHECK_MESSAGE(true, s.str() + s2.str() + s3.str() + s4.str() + ss.str());
}

//todo strange test
BOOST_AUTO_TEST_CASE( TestSequenceRefCount2 ) {
    Sequence *s = new Sequence("AAAAAAA");
    Sequence *s2 = new Sequence(*s);
    Sequence *s3 = new Sequence(!(*s));
    Sequence *ss = new Sequence(s->Subseq(3));
    BOOST_CHECK_MESSAGE(true, s->str() << s2->str() << s3->str() << ss->str());
    delete s;
    delete s2;
    delete s3;
    delete ss;
}

//todo is it suitable here???
//BOOST_AUTO_TEST_CASE( TestSequenceMemory ) {
//    time_t now = time(NULL);
//    int N = 1000;
//    int SIZE = 300;
//    vector<Sequence*> vs(N);
//    unsigned long vm1, vm2, vm3;
//      long rss1, rss2, rss3;
//    process_mem_usage(vm1, rss1);
//    for (int i = 0; i < N; ++i) {
//        string s(SIZE,'-');
//        for (int j = 0; j < SIZE; ++j) {
//            s[j] = nucl(rand() % 4);
//        }
//        vs[i] = new Sequence(s);
//    }
//    process_mem_usage(vm2, rss2);
//    cout << "Memory after creation for " <<  N << " Sequences of size " << SIZE << ": VM = "<< (vm2 - vm1) << " KB., RSS = "<< (rss2 - rss1) << " KB." << endl;
//    for (int i = 0; i < N; ++i) {
//        delete vs[i];
//        vs[i] = 0;
//    }
//    process_mem_usage(vm3, rss3);
//    cout << "Memory after deletion for " <<  N << " Sequences of size " << SIZE << ": VM = "<< (vm3 - vm1) << " KB., RSS = "<< (rss3 - rss1) << " KB." << endl;
//    cout << "Time: " <<  time(NULL) - now << endl;
//}
//

