//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "cute.h"
#include "hash.hpp"
#include "seq.hpp"
#include "sequence.hpp"

using namespace hashing;

int const MPSIZE = 100;

void TestIt() {
	HashSym<Sequence> hh;
	Sequence s("AACATGCTGCACTGGGTATGCATGACTGCAATTATACGCGCGCTACGATCATTACGGTATCATGACATTCATCGGATCATCGTACTGCATCGTATAGATCACATATGATCATATACCTTC");
	s = s.Subseq(0, MPSIZE);
	ASSERT_EQUAL(hh(s), hh(s));
	ASSERT_EQUAL(hh(s), hh(!s));
	Sequence t = s.Subseq(K, 2 * K);
	ASSERT_EQUAL(hh(t), hh(t));
	ASSERT_EQUAL(hh(t), hh(!t));
  /* it's impossible to get ha.size() in assert in kmers()... 	
  hash_t ha[MPSIZE - K + 1];
	hh.kmers(s, ha);
	for (int i = 0; i <= MPSIZE - K; i++) {
		ASSERT_EQUAL(ha[i], hh(s.Subseq(i, i + K)));
	}
  */
}

cute::suite HashSuite(){
	cute::suite s;
	s.push_back(CUTE(TestIt));
	return s;
}
