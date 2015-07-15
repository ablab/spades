//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef TEST_PAIREDREADTEST_HPP_
#define TEST_PAIREDREADTEST_HPP_

#include "cute/cute.h"
#include "io/paired_read.hpp"

using namespace io;

void TestPairedRead() {
  SingleRead sr1("Read1", "ATGCATGC", "aabbaabb");
  SingleRead sr2("Read2", "AATTGGCC", "aabbaabb");
  PairedRead pr1(sr1, sr2, 100);
  ASSERT_EQUAL(true, pr1.IsValid());
  ASSERT_EQUAL(sr1, pr1[0]);
  ASSERT_EQUAL(sr2, pr1[1]);
  PairedRead pr2(!sr2, !sr1, 100);
  ASSERT_EQUAL(pr2, !pr1);
}

cute::suite PairedReadSuite() {
  cute::suite s;
  s.push_back(CUTE(TestPairedRead));
  return s;
}

#endif /* TEST_PAIREDREADTEST_HPP_ */
