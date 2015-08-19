//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef TEST_SAMBAMPARSERTEST_HPP
#define TEST_SAMBAMPARSERTEST_HPP

#include "cute/cute.h"
#include "io/sam_bam_parser.hpp"

using namespace io;

void TestSamBamParserReading() {
  SamBamParser parser("./src/common_test/data/ex1.sam");
  ASSERT(parser.is_open());
  SingleRead read;
  parser >> read;
  ASSERT_EQUAL("B7_591:4:96:693:509", read.name());
  ASSERT_EQUAL("CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCG", 
               read.GetSequenceString());
  ASSERT_EQUAL("<<<<<<<<<<<<<<<;<<<<<<<<<5<<<<<;:<;7", 
               read.GetPhredQualityString());
  parser >> read;
  ASSERT(!parser.eof());
}

void TestSamBamParserFull() {
  SamBamParser parser("./src/common_test/data/ex1.sam");
  ASSERT(parser.is_open());
  parser.reset();
  ASSERT(parser.is_open());
  SingleRead read;
  while (!parser.eof()) {
    parser >> read;
  }
  ASSERT_EQUAL("EAS114_26:7:37:79:581", read.name());
  ASSERT_EQUAL("TTTTTTTTTTTTTTTTTTTTTTTCATGCCAGAAAA", 
               read.GetSequenceString());
  ASSERT_EQUAL("3,,,===6===<===<;=====-============", 
               read.GetPhredQualityString());
  parser.close();
  ASSERT(!parser.is_open());
}

cute::suite SamBamParserSuite() {
  cute::suite s;
  s.push_back(CUTE(TestSamBamParserReading));
  s.push_back(CUTE(TestSamBamParserFull));
  return s;
}

#endif /* TEST_SAMBAMPARSERTEST_HPP */
