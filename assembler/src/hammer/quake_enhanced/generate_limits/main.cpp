//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <stdint.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

using std::unordered_map;
using std::min;
using std::max;
using std::vector;
using std::string;

namespace {
/**
 * @variable Length of string buffer which will store k-mer.
 */
const uint32_t kMaxK = 100;

struct Options {
  string ifile;
  string ofile;
  long double threshold;
  uint32_t k;
  bool valid;
  Options()
      : ifile(""),
        ofile(""),
        threshold(0),
        valid(true) {}
};

void PrintHelp(char *progname) {
  printf("Usage: %s ifile.hist ofile.limits threshold k\n", progname);
  printf("Where:\n");
  printf("\tkmer.hist\tfile with histogram statistics for trusted k[q]-mer distribution\n");
  printf("\tofile.limits\tq-value limits for every k-value will be outputted here\n");
  printf("\tthreshold\t k-mer will be treated as erroneous if probability of being erroneous is greater then threshold\n");
  printf("\tk\t k-mer size\n");
}

Options ParseOptions(int argc, char *argv[]) {
  Options ret;
  if (argc != 5) {
    ret.valid = false;
  } else {
    ret.ifile = argv[1];
    ret.ofile = argv[2];
    sscanf(argv[3], "%Lf", &ret.threshold);
    ret.k = atoi(argv[4]);
  }
  return ret;
}
}

long double Factorial(int n) {
  if (n == 0) {
    return 1;
  }
  static unordered_map<int, long double> ans;
  if (ans.count(n) == 0) {
    ans[n] = Factorial(n - 1) * n;
  }
  return ans[n];
}

long double CNK(int n, int k) {
  return Factorial(n) / (Factorial(k) * Factorial(n - k));
}

long double Bernoulli(int k, int n, long double p) {
  return pow(p, k) * pow(1 - p, n - k) * CNK(n, k);
}

long double GetLimit(uint32_t x, const vector<uint32_t> &hist,
                     long double threshold, uint64_t total, uint32_t k) {
  long double l = 0.5;
  long double r = 1;
  const long double eps = 1e-10;
  while (r - l > eps) {
    long double m = l + (r - l) / 2;
    long double err_prob = 0;
    long double k_prob = 0;
    for (uint32_t i = x; i < hist.size(); ++i) {
      long double de = k * Bernoulli(x, i, (1 - m) / k) * hist[i];
      long double dk = Bernoulli(x, i, m) * hist[i];
      if (de != de || dk != dk) {
        continue;
      }
      err_prob += de;
      k_prob += dk;
    }
    //    k_prob = hist[x];
    if (k_prob < threshold * err_prob) {
      l = m;
    } else {
      r = m;
    }
  }
  if (l + 2 * eps > 1) {
    return 1;
  }
  return l;
}

int main(int argc, char *argv[]) {
  Options opts = ParseOptions(argc, argv);
  if (!opts.valid) {
    PrintHelp(argv[0]);
    return 1;
  }
  FILE *ifile = fopen(opts.ifile.c_str(), "r");
  FILE *ofile = fopen(opts.ofile.c_str(), "w");
  vector<uint32_t> hist;
  uint64_t total = 0;
  float x;
  uint32_t count;
  while (fscanf(ifile, "%f %u", &x, &count) == 2) {
    uint32_t r = (uint32_t)(x + 0.5);
    if (r >= hist.size()) {
      hist.resize(r * 1.5 + 1);
    }
    hist[r] += count;
    total += count;
  }
  for (uint32_t i = 0; i < 25; ++i) {
    fprintf(ofile, "%d %Lf\n", i,
            GetLimit(i, hist, opts.threshold, total, opts.k));
  }
  fclose(ifile);
  fclose(ofile);
  return 0;
}
