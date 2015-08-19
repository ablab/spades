//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef TEST_MULTIFILEREADERTEST_HPP_
#define TEST_MULTIFILEREADERTEST_HPP_

#include <vector>
#include <string>
#include "cute/cute.h"
#include "io/single_read.hpp"
#include "io/multifile_reader.hpp"

using namespace io;

void TestMultifileReaderNoFile() {
  std::vector<SingleRead::FilenameType> filenames;
  filenames.push_back("./no-file.fa");
  filenames.push_back("./src/common_test/data/s_test.fastq.gz");
  MultifileReader<SingleRead> reader(filenames);
  ASSERT(reader.is_open());
  filenames.pop_back();
  filenames.push_back("./no-file.fa");
  MultifileReader<SingleRead> reader2(filenames);
  ASSERT(!reader2.is_open());
}

void TestMultifileReaderReadingFrom1File() {
  std::vector<SingleRead::FilenameType> filenames;
  filenames.push_back("./src/common_test/data/s_test_2.fastq.gz");
  MultifileReader<SingleRead> reader(filenames);
  ASSERT(reader.is_open());
  ASSERT(!reader.eof());
  SingleRead read;
  reader >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_1700/1", read.name());
  ASSERT_EQUAL("AAAAAAAAAAAC", read.GetSequenceString());
  ASSERT_EQUAL("GGGGCGGGGEGG", read.GetPhredQualityString());
  reader >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_468/1", read.name());
  ASSERT_EQUAL("TGTGTGTGTGTG", read.GetSequenceString());
  ASSERT_EQUAL("DADDA8<?>@HH", read.GetPhredQualityString());
  ASSERT(reader.eof());
}

void TestMultifileReaderReadingFrom2Files() {
  std::vector<SingleRead::FilenameType> filenames;
  filenames.push_back("./src/common_test/data/s_test.fastq.gz");
  filenames.push_back("./no-file.fa");
  filenames.push_back("./src/common_test/data/s_test_2.fastq.gz");
  MultifileReader<SingleRead> reader(filenames);
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
  reader >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_1700/1", read.name());
  ASSERT_EQUAL("AAAAAAAAAAAC", read.GetSequenceString());
  ASSERT_EQUAL("GGGGCGGGGEGG", read.GetPhredQualityString());
  reader >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_468/1", read.name());
  ASSERT_EQUAL("TGTGTGTGTGTG", read.GetSequenceString());
  ASSERT_EQUAL("DADDA8<?>@HH", read.GetPhredQualityString());
  ASSERT(reader.eof());
}

cute::suite MultifileReaderSuite() {
  cute::suite s;
  s.push_back(CUTE(TestMultifileReaderNoFile));
  s.push_back(CUTE(TestMultifileReaderReadingFrom1File));
  s.push_back(CUTE(TestMultifileReaderReadingFrom2Files));
  return s;
}

#endif /* TEST_MULTIFILEREADERTEST_HPP_ */
