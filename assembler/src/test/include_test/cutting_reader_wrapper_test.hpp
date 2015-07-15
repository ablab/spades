//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef TEST_CUTTINGREADERWRAPPERTEST_HPP_
#define TEST_CUTTINGREADERWRAPPERTEST_HPP_

#include "cute/cute.h"
#include "io/single_read.hpp"
#include "io/reader.hpp"
#include "io/cutting_reader_wrapper.hpp"

using namespace io;

void TestCuttingReaderWrapperNoFile() {
  Reader<SingleRead> internal_reader("./no-file.fa");
  CuttingReaderWrapper<SingleRead> reader(internal_reader);
  ASSERT(!reader.is_open());
}

void TestCuttingReaderWrapperReading() {
  Reader<SingleRead> internal_reader("./src/common_test/data/s_test.fastq.gz");
  CuttingReaderWrapper<SingleRead> reader(internal_reader);
  ASSERT(reader.is_open());
  ASSERT(!reader.eof());
  SingleRead read;
  reader >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_768/1", read.name());
  ASSERT_EQUAL("ATGCATGCATGC", read.GetSequenceString());
  ASSERT_EQUAL("HGHIHHHGHECH", read.GetPhredQualityString());
  reader >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_1700/1", read.name());
  ASSERT_EQUAL("AAAAAAAAAAAC", read.GetSequenceString());
  ASSERT_EQUAL("GGGGCGGGGEGG", read.GetPhredQualityString());
  reader >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_468/1", read.name());
  ASSERT_EQUAL("TGTGTGTGTGTG", read.GetSequenceString());
  ASSERT_EQUAL("DADDA8<?>@HH", read.GetPhredQualityString());
  ASSERT(reader.eof());
  internal_reader.reset();
  CuttingReaderWrapper<SingleRead> reader2(internal_reader, 1);
  reader2 >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_768/1", read.name());
  ASSERT_EQUAL("ATGCATGCATGC", read.GetSequenceString());
  ASSERT_EQUAL("HGHIHHHGHECH", read.GetPhredQualityString());
  ASSERT(reader2.eof());
}

cute::suite CuttingReaderWrapperSuite() {
  cute::suite s;
  s.push_back(CUTE(TestCuttingReaderWrapperNoFile));
  s.push_back(CUTE(TestCuttingReaderWrapperReading));
  return s;
}

#endif /* TEST_CUTTINGREADERWRAPPERTEST_HPP_ */
