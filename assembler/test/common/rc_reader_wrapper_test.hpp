#ifndef TEST_RCREADERWRAPPERTEST_HPP_
#define TEST_RCREADERWRAPPERTEST_HPP_

#include "cute/cute.h"
#include "common/io/single_read.hpp"
#include "common/io/reader.hpp"
#include "common/io/rc_reader_wrapper.hpp"

void TestRCReaderWrapperNoFile() {
  Reader<SingleRead> internal_reader("./no-file");
  RCReaderWrapper<SingleRead> reader(&internal_reader);
  ASSERT(!reader.is_open());
}

void TestRCReaderWrapperReading() {
  Reader<SingleRead> internal_reader("./test/data/s_test.fastq.gz");
  RCReaderWrapper<SingleRead> reader(&internal_reader);
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

cute::suite RCReaderWrapperSuite(){
  cute::suite s;
  s.push_back(CUTE(TestRCReaderWrapperNoFile));
  s.push_back(CUTE(TestRCReaderWrapperReading));
  return s;
}

#endif /* TEST_RCREADERWRAPPERTEST_HPP_ */
