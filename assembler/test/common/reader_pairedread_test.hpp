#ifndef TEST_READERPAIREDREAD_HPP_
#define TEST_READERPAIREDREAD_HPP_

#include <utility>
#include <string>
#include "cute/cute.h"
#include "common/io/reader.hpp"

using namespace io;

void TestReaderPairedReadNoFile() {
  Reader<PairedRead> reader(std::pair<std::string, std::string>
                            ("./no-file", "./test/data/s_test.fastq.gz"),
                            100, 33);
  ASSERT(!reader.is_open());
  Reader<PairedRead> reader2(std::pair<std::string, std::string>
                             ("./no-file", "./no-file"), 100, 33);
  ASSERT(!reader.is_open());
}

void TestReaderPairedReadReading() {
  Reader<PairedRead> reader(std::pair<std::string, std::string>
                            ("./test/data/s_test.fastq.gz",
                             "./test/data/s_test_2.fastq.gz"), 100, 33);
  ASSERT(reader.is_open());
  PairedRead read;
  reader >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_768/1", read[0].name());
  ASSERT_EQUAL("ATGCATGCATGC", read[0].GetSequenceString());
  ASSERT_EQUAL("HGHIHHHGHECH", read[0].GetPhredQualityString());
  ASSERT_EQUAL("!EAS20_8_6_1_2_1700/1", read[1].name());
  ASSERT_EQUAL(100, read.distance());
  reader >> read;
  ASSERT_EQUAL("EAS20_8_6_1_2_1700/1", read[0].name());
  ASSERT_EQUAL("AAAAAAAAAAAC", read[0].GetSequenceString());
  ASSERT_EQUAL("GGGGCGGGGEGG", read[0].GetPhredQualityString());
  ASSERT_EQUAL("!EAS20_8_6_1_2_468/1", read[1].name());
  ASSERT_EQUAL(100, read.distance());
  ASSERT(reader.eof());
  reader.reset();
  ASSERT(!reader.eof());
  ASSERT(reader.is_open());
  reader.close();
  ASSERT(reader.eof());
  ASSERT(!reader.is_open());
}

cute::suite ReaderPairedReadSuite() {
  cute::suite s;
  s.push_back(CUTE(TestReaderPairedReadNoFile));
  s.push_back(CUTE(TestReaderPairedReadReading));
  return s;
}

#endif /* TEST_READERPAIREDREAD_HPP_ */
