#ifndef TEST_PARSERTEST_HPP_
#define TEST_PARSERTEST_HPP_

#include "cute/cute.h"
#include "common/io/parser.hpp"
#include "common/io/fasta_fastq_gz_parser.hpp"
#include "common/io/parser.cpp"  // DURTY HACK!!!

using namespace io;

void TestGetExtension() {
  ASSERT_EQUAL("fastq", GetExtension("test.fastq"));
  ASSERT_EQUAL("fastq.gz", GetExtension("test.fastq.gz"));
  ASSERT_EQUAL("", GetExtension("README"));
}

void TestParserSelector() {
  ASSERT_EQUAL(reinterpret_cast<Parser*>(NULL), SelectParser("README"));
}

cute::suite ParserSuite() {
  cute::suite s;
  s.push_back(CUTE(TestGetExtension));
  s.push_back(CUTE(TestParserSelector));
  return s;
}

#endif /* TEST_FASTQGZPARSERTEST_HPP_ */
