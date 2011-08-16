#ifndef TEST_SFFPARSERTEST_HPP
#define TEST_SFFPARSERTEST_HPP

#include "cute/cute.h"
#include "common/io/sff_parser.hpp"

using namespace io;

void TestSffParserReading() {
  SffParser parser("./test/data/test.sff", 33);
  ASSERT(parser.is_open());
  SingleRead read;
  parser >> read;
  ASSERT_EQUAL("Name", read.name());
  ASSERT_EQUAL("GACTCTTTCATTTCCTACTGTAGCTTTTAGTCTCTTCAAATACAAGGCACACAGGGATAGG", 
               read.GetSequenceString());
  ASSERT_EQUAL("", 
               read.GetPhredQualityString());
  parser >> read;
  ASSERT(!parser.eof());
}

void TestSffParserFull() {
  SffParser parser("./test/data/test.sff", 33);
  ASSERT(parser.is_open());
  parser.reset();
  ASSERT(parser.is_open());
  SingleRead read;
  while (!parser.eof()) {
    parser >> read;
  }
  ASSERT_EQUAL("Name", read.name());
  ASSERT_EQUAL("GACTCTTTCATTTCCTACTGTAGCTTTTAGTCTCTTCAAATACAAGGCACACAGGGATAGG", 
               read.GetSequenceString());
  ASSERT_EQUAL("", 
               read.GetPhredQualityString());
  parser.close();
  ASSERT(!parser.is_open());
}

cute::suite SffParserSuite() {
  cute::suite s;
  s.push_back(CUTE(TestSffParserReading));
  s.push_back(CUTE(TestSffParserFull));
  return s;
}

#endif /* TEST_SFFPARSERTEST_HPP */
