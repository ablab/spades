//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <stdint.h>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <string>
#include <vector>

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
  string kmer_hist;
  string trusted_hist;
  string bad_hist;
  size_t level;
  bool valid;
  Options()
      : kmer_hist(""),
        trusted_hist(""),
        bad_hist(""),
        level(0),
        valid(true) {}
};

void PrintHelp(char *progname) {
  printf("Usage: %s kmer.hist trusted.hist bad.hist level\n", progname);
  printf("Where:\n");
  printf("\tkmer.hist\tfile with histogram statistics for k[q]-mer distribution\n");
  printf("\ttrusted.hist\thistogram for trusted k-mers will be outputted here\n");
  printf("\tbad.hist\thistogram for bad k-mers will be outputted here\n");
  printf("\tlevel\tGauss must be at least level times higher then previous local min\n");
}

Options ParseOptions(int argc, char *argv[]) {
  Options ret;
  if (argc != 5) {
    ret.valid = false;
  } else {
    ret.kmer_hist = argv[1];
    ret.trusted_hist = argv[2];
    ret.bad_hist = argv[3];
    ret.level = atoi(argv[4]);
  }
  return ret;
}
}

int main(int argc, char *argv[]) {
  Options opts = ParseOptions(argc, argv);
  if (!opts.valid) {
    PrintHelp(argv[0]);
    return 1;
  }
  FILE *kmer_hist = fopen(opts.kmer_hist.c_str(), "r");
  FILE *trusted_hist = fopen(opts.trusted_hist.c_str(), "w");
  FILE *bad_hist = fopen(opts.bad_hist.c_str(), "w");
  vector<uint32_t> hist;
  float x;
  int count;
  while (fscanf(kmer_hist, "%f %d", &x, &count) == 2) {
    uint32_t r = (uint32_t)(x + 0.5);
    if (r >= hist.size()) {
      hist.resize(r * 1.5 + 1);
    }
    hist[r] += count;
  }
  int fmin = -1;
  for (uint32_t i = 1; i < hist.size() - 1; ++i) {
    if (hist[i + 1] > hist[i]) {
      fmin = i;
      break;
    }
  }
  if (fmin == -1) {
    printf("Bad histogram");
    return 0;
  }
  int fmax = -1;
  for (uint32_t i = fmin; i < hist.size() - 1; ++i) {
    if (hist[i + 1] < hist[i] && hist[i] > hist[fmin] * opts.level) {
      fmax = i;
      break;
    }
  }
  if (fmax == -1) {
    printf("Bad histogram");
    return 0;
  }
  int lborder = fmax;
  int rborder = fmax;
  // ToDo 0.9 ==>> command line options
  while (hist[lborder] > hist[fmax] * 0.9) {
    --lborder;
  }
  while (hist[rborder] > hist[fmax] * 0.9) {
    ++rborder;
  }

  uint32_t mass = 0;
  uint64_t mass_pos = 0;

  for (int i = lborder; i <= rborder; ++i) {
    mass += hist[i];
    mass_pos += hist[i] * i;
  }
  float average = mass_pos / static_cast<double>(mass);
  printf("Gauss median is at %f\n", average);
  vector<uint32_t> hist_trusted(hist);
  for (uint32_t i = 0; static_cast<int>(average - i) >=0; ++i) {
    int where = static_cast<int>(average - i);
    int from = static_cast<int>(average + i + 0.5);
    if (where == from) {
      continue;
    }
    hist_trusted[where] = min(hist_trusted[where + 1], hist_trusted[from]);
  }
  for (uint32_t i = 0; i < hist_trusted.size(); ++i) {
    fprintf(trusted_hist, "%d %d\n", i, hist_trusted[i]);
    fprintf(bad_hist, "%d %d\n", i, hist[i] - hist_trusted[i]);
  }
  fclose(kmer_hist);
  fclose(trusted_hist);
  fclose(bad_hist);
  return 0;
}
