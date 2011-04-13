#include "cute.h"
#include "hash.hpp"
#include "seq.hpp"
#include "sequence.hpp"

void TestIt() {
	HashSym<Sequence> hh;
	Sequence s("AACATGCTGCACTGGGT");
	ASSERT_EQUAL(hh(s), hh(s));
	Sequence t = s.Subseq(1, 6);
	ASSERT_EQUAL(hh(t), hh(t));
	ASSERT_EQUAL(hh(t), hh(!t));
	ASSERT_EQUAL(hh(s), hh(!s));
}

cute::suite HashSuite(){
	cute::suite s;
	s.push_back(CUTE(TestIt));
	return s;
}
