#ifndef TEST_MULTIFILEREADERWRAPPERTEST_HPP_
#define TEST_MULTIFILEREADERWRAPPERTEST_HPP_

#include <vector>
#include <string>
#include "cute/cute.h"
#include "common/io/single_read.hpp"
#include "common/io/reader.hpp"
#include "common/io/multifile_reader_wrapper.hpp"

void TestMultifileReaderWrapperNoFile() {
  std::vector<SingleRead::FilenameType> filenames;
  filenames.push_back("./no-file");
  filenames.push_back("./test/data/s_test.fastq.gz");
  MultifileReaderWrapper<SingleRead> reader(filenames);
  ASSERT(reader.is_open());
  filenames.pop_back();
  filenames.push_back("./no-file");
  MultifileReaderWrapper<SingleRead> reader2(filenames);
  ASSERT(!reader2.is_open());
}

void TestMultifileReaderWrapperReadingFrom1File() {
  std::vector<SingleRead::FilenameType> filenames;
  filenames.push_back("./test/data/s_test_2.fastq.gz");
  MultifileReaderWrapper<SingleRead> reader(filenames);
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

void TestMultifileReaderWrapperReadingFrom2Files() {
  std::vector<SingleRead::FilenameType> filenames;
  filenames.push_back("./test/data/s_test.fastq.gz");
  filenames.push_back("./no-file");
  filenames.push_back("./test/data/s_test_2.fastq.gz");
  MultifileReaderWrapper<SingleRead> reader(filenames);
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

cute::suite MultifileReaderWrapperSuite(){
  cute::suite s;
  s.push_back(CUTE(TestMultifileReaderWrapperNoFile));
  s.push_back(CUTE(TestMultifileReaderWrapperReadingFrom1File));
  s.push_back(CUTE(TestMultifileReaderWrapperReadingFrom2Files));
  return s;
}

#endif /* TEST_MULTIFILEREADERWRAPPERTEST_HPP_ */
