//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef TEST_SINGLEREADTEST_HPP_
#define TEST_SINGLEREADTEST_HPP_

#include "cute/cute.h"
#include "io/single_read.hpp"

using namespace io;

void TestSingleRead() {
  SingleRead sr("Read1", "ATGCATGC", "\3\3\4\4\3\3\4\4");
  ASSERT_EQUAL(true, sr.IsValid());
  ASSERT_EQUAL(8, sr.size());
  ASSERT_EQUAL("$$%%$$%%", sr.GetPhredQualityString());
  ASSERT_EQUAL(dignucl('A'), sr[0]);
  ASSERT_EQUAL(dignucl('T'), sr[5]);
  ASSERT_EQUAL("GCATGCAT", (!sr).GetSequenceString());
}

cute::suite SingleReadSuite() {
  cute::suite s;
  s.push_back(CUTE(TestSingleRead));
  return s;
}
#endif /* TEST_SINGLEREADTEST_HPP_ */
