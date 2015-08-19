//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef HAMMER_KMERFUNCTIONSTEST_HPP_
#define HAMMER_KMERFUNCTIONSTEST_HPP

#include <cmath>
#include <map>
#include "cute/cute.h"
#include "io/read.hpp"
#include "sequence/seq.hpp"
#include "valid_kmer_generator.hpp"

void TestValidKMerGenerator() {
  Read r("TestRead1", "BACNTACGT", "\1\3\2\3\3\2\1\3\1");
  ValidKMerGenerator<2> gen(r, 2);
  const double eps = 1e-10;
  ASSERT(gen.HasMore());
  ASSERT_EQUAL("AC", gen.kmer().str());
  ASSERT(abs(gen.correct_probability() - 0.18408318790937242) < eps)
  gen.Next();
  ASSERT(gen.HasMore());
  ASSERT_EQUAL("TA", gen.kmer().str());
  ASSERT(abs(gen.correct_probability() - 0.18408318790937242) < eps)
  gen.Next();
  ASSERT(gen.HasMore());
  ASSERT_EQUAL("AC", gen.kmer().str());
  ASSERT(abs(gen.correct_probability() - 0.07590165442279753) < eps)
  gen.Next();
  ASSERT(gen.HasMore());
  ASSERT_EQUAL("CG", gen.kmer().str());
  ASSERT(abs(gen.correct_probability() - 0.10259170220194348) < eps)
  gen.Next();
  ASSERT(!gen.HasMore());
}

cute::suite ValidKMerGeneratorSuite() {
  cute::suite s;
  s.push_back(CUTE(TestValidKMerGenerator));
  return s;
}

#endif  // HAMMER_KMERFUNCTIONSTEST_HPP_
