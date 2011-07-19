#ifndef TEST_FASTQGZPARSERTEST_HPP_
#define TEST_FASTQGZPARSERTEST_HPP_

#include "cute/cute.h"
#include "common/io/fastqgz_parser.hpp"

void TestFastqgzParserNoFile() {
  FastqgzParser parser("./no-file");
  ASSERT(!parser.is_open());
}

void TestFastqgzParserReading() {
  FastqgzParser parser("./test/data/s_test.fastq.gz", 33);
  ASSERT(parser.is_open());
  SingleRead read;
  parser >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_768/1", read.name());
  ASSERT_EQUAL("ATGCATGCATGC", read.GetSequenceString());
  ASSERT_EQUAL("HGHIHHHGHECH", read.GetPhredQualityString());
  parser >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_1700/1", read.name());
  ASSERT_EQUAL("AAAAAAAAAAAC", read.GetSequenceString());
  ASSERT_EQUAL("GGGGCGGGGEGG", read.GetPhredQualityString());
  parser >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_468/1", read.name());
  ASSERT_EQUAL("TGTGTGTGTGTG", read.GetSequenceString());
  ASSERT_EQUAL("DADDA8<?>@HH", read.GetPhredQualityString());
  ASSERT(parser.eof());
}

void TestFastqgzParserFull() {
  FastqgzParser parser("./test/data/s_test.fastq.gz", 33);
  ASSERT(parser.is_open());
  SingleRead read;
  while (!parser.eof()) {
    parser >> read;
  }
  ASSERT_EQUAL("EAS20_8_6_1_2_468/1", read.name());
  ASSERT_EQUAL("TGTGTGTGTGTG", read.GetSequenceString());
  ASSERT_EQUAL("DADDA8<?>@HH", read.GetPhredQualityString());
}

cute::suite FastqgzParserSuite(){
  cute::suite s;
  s.push_back(CUTE(TestFastqgzParserNoFile));
  s.push_back(CUTE(TestFastqgzParserReading));
  s.push_back(CUTE(TestFastqgzParserFull));
  return s;
}

#endif /* TEST_FASTQGZPARSERTEST_HPP_ */
