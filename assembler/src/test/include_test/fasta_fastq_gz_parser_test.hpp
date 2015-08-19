//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef TEST_FASTAFASTQGZPARSERTEST_HPP_
#define TEST_FASTAFASTQGZPARSERTEST_HPP_

#include "cute/cute.h"
#include "io/fasta_fastq_gz_parser.hpp"

using namespace io;

void TestFastaFastqGzParserNoFile() {
  FastaFastqGzParser parser("./no-file");
  ASSERT(!parser.is_open());
}

void TestFastaFastqGzParserReading() {
  FastaFastqGzParser parser("./src/common_test/data/s_test.fastq.gz");
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

void TestFastaFastqGzParserFull() {
  FastaFastqGzParser parser("./src/common_test/data/s_test.fastq.gz");
  ASSERT(parser.is_open());
  parser.reset();
  ASSERT(parser.is_open());
  SingleRead read;
  while (!parser.eof()) {
    parser >> read;
  }
  ASSERT_EQUAL("EAS20_8_6_1_2_468/1", read.name());
  ASSERT_EQUAL("TGTGTGTGTGTG", read.GetSequenceString());
  ASSERT_EQUAL("DADDA8<?>@HH", read.GetPhredQualityString());
  parser.close();
  ASSERT(!parser.is_open());
}

void TestFastaFastqGzParserFastaReading() {
  FastaFastqGzParser parser("./src/common_test/data/test.fasta", SolexaOffset);
  ASSERT(parser.is_open());
  SingleRead read;
  parser >> read;
  ASSERT_EQUAL("GSV1ISZ08GSHS3", read.name());
  ASSERT_EQUAL("CTTTCATTTCCTACTGTAGCTTTTAGTCTCTTCAAATACAAGGCACACA", 
               read.GetSequenceString());
  ASSERT_EQUAL("#################################################", 
               read.GetPhredQualityString());
  parser >> read;
  parser >> read;
  parser >> read;
  ASSERT_EQUAL("GSV1ISZ08GXRMP", read.name());
  ASSERT_EQUAL("CTTTCATTTCCTACTGTAGCTTTTAGTCTCTTCAAATACAAGGCACACAGGGAG", 
               read.GetSequenceString());
  ASSERT_EQUAL("######################################################", 
               read.GetPhredQualityString());
  Sequence seq;
  seq = read.sequence();
  Quality qual(read.quality());
  ASSERT(parser.eof());
  parser.close();
}

cute::suite FastaFastqGzParserSuite() {
  cute::suite s;
  s.push_back(CUTE(TestFastaFastqGzParserNoFile));
  s.push_back(CUTE(TestFastaFastqGzParserReading));
  s.push_back(CUTE(TestFastaFastqGzParserFull));
  s.push_back(CUTE(TestFastaFastqGzParserFastaReading));
  return s;
}

#endif /* TEST_FASTAFASTQGZPARSERTEST_HPP_ */
