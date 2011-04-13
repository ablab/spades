#include "cute.h"
#include "hash.hpp"
#include "seq.hpp"
#include "sequence.hpp"

void TestIt() {
	HashSym<Sequence> hh;
	Sequence s("AACATGCTGCACTGGGTATGCATGACTGCAATTATACGCGCGCTACGATCATTACGGTATCATGACATTCATCGGATCATCGTACTGCATCGTATAGATCACATATGATCATATACCTTC");
	s = s.Subseq(0, MPSIZE);
	ASSERT_EQUAL(hh(s), hh(s));
	ASSERT_EQUAL(hh(s), hh(!s));
	Sequence t = s.Subseq(K, 2 * K);
	ASSERT_EQUAL(hh(t), hh(t));
	ASSERT_EQUAL(hh(t), hh(!t));
	hash_t ha[MPSIZE - K + 1];
	hh.kmers(s, ha);
	for (int i = 0; i <= MPSIZE - K; i++) {
		ASSERT_EQUAL(ha[i], hh(s.Subseq(i, i + K)));
	}
}

cute::suite HashSuite(){
	cute::suite s;
	s.push_back(CUTE(TestIt));
	return s;
}
