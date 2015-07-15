//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef READGTEST_HPP_
#define READGTEST_HPP_
#include "cute/cute.h"
#include "io/read.hpp"

void TestGetSubsequence() {
  Read r("TestRead1", "ACGTACGT", "BBBBBBBB");
  ASSERT_EQUAL("CG", r.getSubSequence(1, 2).str());
  ASSERT_EQUAL("AC", r.getSubSequence(0, 2).str());
  ASSERT_EQUAL("A", r.getSubSequence(0, 1).str());
  ASSERT_EQUAL("ACGTACGT", r.getSubSequence(0, 8).str());
}

void TestTrimBadQuality() {
  Read r("TestRead1", "ACGTACGT", "\1\2\3\3\2\1\3\1");
  ASSERT_EQUAL(5, r.trimBadQuality());
  ASSERT_EQUAL("GTACG", r.getSequenceString());
  ASSERT_EQUAL("\3\3\2\1\3", r.getQualityString());
  Read r2("TestRead2", "ACGTACGT", "\1\2\2\2\2\1\1\1");
  ASSERT_EQUAL(0, r2.trimBadQuality());
  ASSERT_EQUAL("", r2.getSequenceString());
}

void TestFirstValidKmer() {
  Read r("TestRead1", "ACGTACGT", "\1\2\3\3\2\1\3\1");
  ASSERT_EQUAL(-1, r.firstValidKmer(0, 9));
  ASSERT_EQUAL(-1, r.firstValidKmer(1, 8));
  ASSERT_EQUAL(0, r.firstValidKmer(0, 7));
  ASSERT_EQUAL(1, r.firstValidKmer(1, 7));
  ASSERT_EQUAL(0, r.firstValidKmer(0, 8));
  Read r2("TestRead1", "ACNTACGT", "\1\2\3\3\2\1\3\1");
  ASSERT_EQUAL(0, r2.firstValidKmer(0, 2));
  ASSERT_EQUAL(3, r2.firstValidKmer(0, 3));
  ASSERT_EQUAL(3, r2.firstValidKmer(1, 2));
  ASSERT_EQUAL(3, r2.firstValidKmer(1, 5));
  ASSERT_EQUAL(3, r2.firstValidKmer(0, 5));
  ASSERT_EQUAL(3, r2.firstValidKmer(2, 5));
  ASSERT_EQUAL(3, r2.firstValidKmer(2, 1));
  ASSERT_EQUAL(-1, r2.firstValidKmer(0, 6));
}

cute::suite ReadSuite() {
  cute::suite s;
  s.push_back(CUTE(TestGetSubsequence));
  s.push_back(CUTE(TestTrimBadQuality));
  s.push_back(CUTE(TestFirstValidKmer));
  return s;
}
#endif /* READTEST_HPP_ */
