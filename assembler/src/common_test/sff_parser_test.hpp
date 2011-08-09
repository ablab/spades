#ifndef TEST_SFFPARSERTEST_HPP
#define TEST_SFFPARSERTEST_HPP

#include "cute/cute.h"
#include "common/io/sff_parser.hpp"

using namespace io;

void TestSffParserNoFile() {
  SffParser parser("./no-file");
  ASSERT(!parser.is_open());
}

void TestSffParserReading() {
  SffParser parser("./test/data/proc.srf", 33);
  ASSERT(parser.is_open());
  SingleRead read;
  parser >> read;
  ASSERT_EQUAL("Read", read.name());
  ASSERT_EQUAL("GTATAAGTCAAAGCACCTTTAGCGTTAAGGTACTGAATCTGTGCCGCGTTTCTTTGTTCCTGAGCATGGCACTA", 
               read.GetSequenceString());
  ASSERT_EQUAL("", 
               read.GetPhredQualityString());
  parser >> read;
  ASSERT(!parser.eof());
}

void TestSffParserFull() {
  // SffParser parser("./test/data/proc.srf", 33);
  // ASSERT(parser.is_open());
  // parser.reset();
  // ASSERT(parser.is_open());
  // SingleRead read;
  // while (!parser.eof()) {
  //   parser >> read;
  // }
  // ASSERT_EQUAL("EAS114_26:7:37:79:581", read.name());
  // ASSERT_EQUAL("TTTTTTTTTTTTTTTTTTTTTTTCATGCCAGAAAA", 
  //              read.GetSequenceString());
  // ASSERT_EQUAL("3,,,===6===<===<;=====-============", 
  //              read.GetPhredQualityString());
  // parser.close();
  // ASSERT(!parser.is_open());
}

cute::suite SffParserSuite() {
  cute::suite s;
  s.push_back(CUTE(TestSffParserNoFile));
  s.push_back(CUTE(TestSffParserReading));
  s.push_back(CUTE(TestSffParserFull));
  return s;
}

#endif /* TEST_SFFPARSERTEST_HPP */
