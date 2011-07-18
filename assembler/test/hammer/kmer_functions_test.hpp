#ifndef HAMMER_KMERFUNCTIONSTEST_HPP_
#define HAMMER_KMERFUNCTIONSTEST_HPP

#include <map>
#include "cute/cute.h"
#include "common/sequence/seq.hpp"
#include "hammer/kmer_functions.hpp"

typedef Seq<2> KMer;
typedef map<KMer, KMerStat, KMer::less2> KMerStatMap;

void TestGetSubsequence() {
  Read r("TestRead1", "ACGTACGT", "BBBBBBBB");
  ASSERT_EQUAL("CG", GetSubSequence(r, 1, 2).str());
  ASSERT_EQUAL("AC", GetSubSequence(r, 0, 2).str());
  ASSERT_EQUAL("A", GetSubSequence(r, 0, 1).str());
  ASSERT_EQUAL("ACGTACGT", GetSubSequence(r, 0, 8).str());
}

void TestTrimBadQuality() {
  Read r("TestRead1", "ACGTACGT", "\1\2\3\3\2\1\3\1");
  ASSERT_EQUAL(5, TrimBadQuality(&r));
  ASSERT_EQUAL("GTACG", r.getSequenceString());
  ASSERT_EQUAL("\3\3\2\1\3", r.getQualityString());
  Read r2("TestRead2", "ACGTACGT", "\1\2\2\2\2\1\1\1");
  ASSERT_EQUAL(0, TrimBadQuality(&r2));
  ASSERT_EQUAL("", r2.getSequenceString());
  ASSERT_EQUAL("", r2.getQualityString());
}

void TestFirstValidKmerPos() {
  Read r("TestRead1", "ACGTACGT", "\1\2\3\3\2\1\3\1");
  ASSERT_EQUAL(r.getSequenceString().size(), FirstValidKmerPos(r, 0, 9));
  ASSERT_EQUAL(r.getSequenceString().size(), FirstValidKmerPos(r, 1, 8));
  ASSERT_EQUAL(0, FirstValidKmerPos(r, 0, 7));
  ASSERT_EQUAL(1, FirstValidKmerPos(r, 1, 7));
  ASSERT_EQUAL(0, FirstValidKmerPos(r, 0, 8));
  Read r2("TestRead1", "ACNTACGT", "\1\2\3\3\2\1\3\1");
  ASSERT_EQUAL(0, FirstValidKmerPos(r2, 0, 2));
  ASSERT_EQUAL(3, FirstValidKmerPos(r2, 0, 3));
  ASSERT_EQUAL(3, FirstValidKmerPos(r2, 1, 2));
  ASSERT_EQUAL(3, FirstValidKmerPos(r2, 1, 5));
  ASSERT_EQUAL(3, FirstValidKmerPos(r2, 0, 5));
  ASSERT_EQUAL(3, FirstValidKmerPos(r2, 2, 5));
  ASSERT_EQUAL(3, FirstValidKmerPos(r2, 2, 1));
  ASSERT_EQUAL(r2.getSequenceString().size(), FirstValidKmerPos(r2, 0, 6));
}

void TestAddKMers() {
  Read r("TestRead1", "ACNTACGT", "\1\2\3\3\2\1\3\1");
  KMerStatMap m, m2;
  AddKMers<2>(r, 666, &m); //666 makes absolutely no sense here
  m2[KMer("AC")].count = 2;
  m2[KMer("TA")].count = 1;
  m2[KMer("CG")].count = 1;
  m2[KMer("GT")].count = 1;
  for (KMerStatMap::iterator it = m.begin(); it != m.end(); ++it) {
    ASSERT_EQUAL(m2[it->first].count, it->second.count);
  }
  for (KMerStatMap::iterator it = m2.begin(); it != m2.end(); ++it) {
    ASSERT_EQUAL(m[it->first].count, it->second.count);
  }
}

cute::suite KMerFunctionsSuite() {
  cute::suite s;
  s.push_back(CUTE(TestGetSubsequence));
  s.push_back(CUTE(TestTrimBadQuality));
  s.push_back(CUTE(TestFirstValidKmerPos));
  s.push_back(CUTE(TestAddKMers));
  return s;
}

#endif  // HAMMER_KMERFUNCTIONSTEST_HPP_
