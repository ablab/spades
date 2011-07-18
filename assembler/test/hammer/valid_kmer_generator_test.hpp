#ifndef HAMMER_KMERFUNCTIONSTEST_HPP_
#define HAMMER_KMERFUNCTIONSTEST_HPP

#include <map>
#include "cute/cute.h"
#include "common/sequence/seq.hpp"
#include "hammer/valid_kmer_generator.hpp"

void TestValidKMerGenerator() {
  Read r("TestRead1", "BACNTACGT", "\1\3\2\3\3\2\1\3\1");
  ValidKMerGenerator<2> gen(r, 2);
  ASSERT(gen.HasMore());
  ASSERT_EQUAL("AC", gen.kmer().str());
  gen.Next();
  ASSERT(gen.HasMore());
  ASSERT_EQUAL("TA", gen.kmer().str());
  gen.Next();
  ASSERT(gen.HasMore());
  ASSERT_EQUAL("AC", gen.kmer().str());
  gen.Next();
  ASSERT(gen.HasMore());
  ASSERT_EQUAL("CG", gen.kmer().str());
  gen.Next();
  ASSERT(!gen.HasMore());
}

cute::suite ValidKMerGeneratorSuite() {
  cute::suite s;
  s.push_back(CUTE(TestValidKMerGenerator));
  return s;
}

#endif  // HAMMER_KMERFUNCTIONSTEST_HPP_
