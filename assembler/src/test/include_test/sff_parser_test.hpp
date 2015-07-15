//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef TEST_SFFPARSERTEST_HPP
#define TEST_SFFPARSERTEST_HPP

#include "cute/cute.h"
#include "io/sff_parser.hpp"

using namespace io;

void TestSffParserReading() {
  SffParser parser("./src/common_test/data/test.sff");
  ASSERT(parser.is_open());
  SingleRead read;
  parser >> read;
  ASSERT_EQUAL("GSV1ISZ08GSHS3", read.name());
  ASSERT_EQUAL("GACTCTTTCATTTCCTACTGTAGCTTTTAGTCTCTTCAAATACAAGGCACACAGGGATAGG", 
               read.GetSequenceString());
  ASSERT_EQUAL("IIIEB=@@GB555>>IIIIIIIHI>>>>IIIIIIIII666IHIHHHIIIIIH>554:AA66", 
               read.GetPhredQualityString());
  parser >> read;
  ASSERT(!parser.eof());
}

void TestSffParserFull() {
  SffParser parser("./src/common_test/data/test.sff");
  ASSERT(parser.is_open());
  parser.reset();
  ASSERT(parser.is_open());
  SingleRead read;
  while (!parser.eof()) {
    parser >> read;
  }
  ASSERT_EQUAL("GSV1ISZ08GXRMP", read.name());
  ASSERT_EQUAL("GACTCTTTCATTTCCTACTGTAGCTTTTAGTCTCTTCAAATACAAGGCACACAGGGAGAGTG", 
               read.GetSequenceString());
  ASSERT_EQUAL("IIHEEEEHHF:99AAHIHHIHIIIIIIGHHBHHIHHH999DCHHHHHHHHHIIG<<757655", 
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
