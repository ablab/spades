//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef TEST_CONVERTINGREADERWRAPPERTEST_HPP_
#define TEST_CONVERTINGREADERWRAPPERTEST_HPP_

#include <utility>
#include <string>
#include "cute/cute.h"
#include "io/single_read.hpp"
#include "io/paired_read.hpp"
#include "io/reader.hpp"
#include "io/converting_reader_wrapper.hpp"

using namespace io;

void TestConvertingReaderWrapperNoFile() {
  Reader<PairedRead> internal_reader(
      std::pair<std::string, std::string>(
          "./no-file.fa", "./no_file.fa"));
  ConvertingReaderWrapper reader(internal_reader);
  ASSERT(!reader.is_open());
}

void TestConvertingReaderWrapperReading() {
  Reader<PairedRead> internal_reader(
      std::pair<std::string, std::string>(
          "./src/common_test/data/s_test.fastq.gz",
          "./src/common_test/data/s_test_2.fastq.gz"));
  ConvertingReaderWrapper reader(internal_reader);
  ASSERT(reader.is_open());
  ASSERT(!reader.eof());
  SingleRead read;
  reader >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_768/1", read.name());
  ASSERT_EQUAL("ATGCATGCATGC", read.GetSequenceString());
  ASSERT_EQUAL("HGHIHHHGHECH", read.GetPhredQualityString());
  reader >> read;
  ASSERT_EQUAL("!EAS20_8_6_1_2_1700/1", read.name());
  ASSERT_EQUAL("GTTTTTTTTTTT", read.GetSequenceString());
  ASSERT_EQUAL("GGEGGGGCGGGG", read.GetPhredQualityString());
  reader >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_1700/1", read.name());
  ASSERT_EQUAL("AAAAAAAAAAAC", read.GetSequenceString());
  ASSERT_EQUAL("GGGGCGGGGEGG", read.GetPhredQualityString());
  reader >> read;
  ASSERT_EQUAL("!EAS20_8_6_1_2_468/1", read.name());
  ASSERT_EQUAL("CACACACACACA", read.GetSequenceString());
  ASSERT_EQUAL("HH@>?<8ADDAD", read.GetPhredQualityString());
  ASSERT(reader.eof());
}

cute::suite ConvertingReaderWrapperSuite() {
  cute::suite s;
  s.push_back(CUTE(TestConvertingReaderWrapperNoFile));
  s.push_back(CUTE(TestConvertingReaderWrapperReading));
  return s;
}

#endif /* TEST_CONVERTINGREADERWRAPPERTEST_HPP_ */
