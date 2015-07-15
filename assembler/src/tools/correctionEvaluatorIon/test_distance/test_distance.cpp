//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "hkmer_distance.hpp"
#include "HSeq.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

struct Tester {
  bool verbose;
  unsigned failed;
  unsigned passed;
  unsigned equal_run_nucls;
  Tester() : verbose(false), failed(0), passed(0), equal_run_nucls(0) {}
  ~Tester() {
    std::cout << "total:  " << failed + passed << std::endl;
    std::cout << "passed: " << passed << std::endl;
    std::cout << "failed: " << failed << std::endl;
    std::cout << std::endl << "no run indels: " << equal_run_nucls << std::endl;
  }
  void operator()(const std::string& s1, const std::string& s2, int expected) {
    std::vector<hammer::HomopolymerRun> k1, k2;
    hammer::iontorrent::toHomopolymerRuns(s1, k1);
    hammer::iontorrent::toHomopolymerRuns(s2, k2);

    bool all_nucls_are_equal = k1.size() == k2.size();
    for (size_t i = 0; i < k1.size(); i++) {
      if (k1[i].nucl != k2[i].nucl) {
        all_nucls_are_equal = false;
        break;
      }
    }

    if (all_nucls_are_equal)
      ++equal_run_nucls;

    int dist = hammer::distanceHKMer(k1.begin(), k1.end(), k2.begin(), k2.end());
    if (dist != expected) {
      ++failed;
      if (verbose) {
        std::cerr << "dist(" << s1 << ", " << s2 << ") = " << dist
                  << " (expected " << expected << ")" << std::endl;
      }
    } else {
      ++passed;
    }
  }
} test;

// Each line in input file must be of form "<kmer1> <kmer2> <expected_distance>"
int main(int argc, char** argv) {
  if (argc < 2) { std::cout << "Usage: ./test_distance <tests.dat>" << std::endl; return 0; }
  test.verbose = true;
  std::ifstream in(argv[1]);
  while (!in.eof()) {
    std::string s1, s2;
    int expected;
    in >> s1 >> s2 >> expected;
    if (s1.empty())
      return 0;
    test(s1, s2, expected);
  }
}
