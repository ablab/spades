#ifndef TEST_SINGLEREADTEST_HPP_
#define TEST_SINGLEREADTEST_HPP_

#include "cute/cute.h"
#include "io/single_read.hpp"

using namespace io;

void TestSingleRead() {
  SingleRead sr("Read1", "ATGCATGC", "aabbaabb");
  ASSERT_EQUAL(true, sr.IsValid());
  ASSERT_EQUAL(8, sr.size());
  ASSERT_EQUAL("ccddccdd", sr.GetPhredQualityString(2));
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
