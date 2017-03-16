//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include <boost/test/unit_test.hpp>
#include "longseq.hpp"
#include "sequence/sequence.hpp"
#include "sequence/nucl.hpp"
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <string>

namespace cap {

typedef unsigned long long ull;

size_t ReadNumTests(const char *filename) {
  ifstream fs(filename);
  size_t res;
  fs >> res;
  return res;
}

void ReadTestData(const std::string &filename, std::string &s1, std::string &s2,
                  vector<pair<int, int> > &inversions) {
  std::ifstream fs(filename);
  size_t L, N;

  fs >> L >> N;
  while (N--) {
    pair<int, int> inversion;
    fs >> inversion.first >> inversion.second;
    inversions.push_back(inversion);
  }

  std::sort(inversions.begin(), inversions.end());

  fs >> s1 >> s2;
}

void ConductInversionTest(size_t testn) {
  typedef std::pair<int, int> pii;

  const std::string &test_filename = "tests/" + std::to_string(testn);
  std::string s1, s2;
  std::vector<pii> inversions;

  ReadTestData(test_filename, s1, s2, inversions);
}

BOOST_AUTO_TEST_CASE( SyntheticInversions ) {
  int T = ReadNumTests("tests/summary");

  for (size_t nt = 1; nt <= T; ++nt) {
    ConductInversionTest(nt);
  }
}

}
