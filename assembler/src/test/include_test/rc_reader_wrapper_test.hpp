//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef TEST_RCREADERWRAPPERTEST_HPP_
#define TEST_RCREADERWRAPPERTEST_HPP_

#include "cute/cute.h"
#include "io/single_read.hpp"
#include "io/reader.hpp"
#include "io/rc_reader_wrapper.hpp"

using namespace io;

void TestRCReaderWrapperNoFile() {
  Reader<SingleRead> internal_reader("./no-file.fa");
  RCReaderWrapper<SingleRead> reader(internal_reader);
  ASSERT(!reader.is_open());
}

void TestRCReaderWrapperReading() {
  Reader<SingleRead> internal_reader("./src/common_test/data/s_test.fastq.gz");
  RCReaderWrapper<SingleRead> reader(internal_reader);
  ASSERT(reader.is_open());
  ASSERT(!reader.eof());
  SingleRead read;
  reader >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_768/1", read.name());
  ASSERT_EQUAL("ATGCATGCATGC", read.GetSequenceString());
  ASSERT_EQUAL("HGHIHHHGHECH", read.GetPhredQualityString());
  reader >> read;
  ASSERT_EQUAL("!EAS20_8_6_1_2_768/1", read.name());
  ASSERT_EQUAL("GCATGCATGCAT", read.GetSequenceString());
  ASSERT_EQUAL("HCEHGHHHIHGH", read.GetPhredQualityString());
  size_t number = 2;
  while (!reader.eof()) {
    reader >> read;
    ++number;
  }
  ASSERT_EQUAL(6, number);
}

cute::suite RCReaderWrapperSuite() {
  cute::suite s;
  s.push_back(CUTE(TestRCReaderWrapperNoFile));
  s.push_back(CUTE(TestRCReaderWrapperReading));
  return s;
}

#endif /* TEST_RCREADERWRAPPERTEST_HPP_ */
