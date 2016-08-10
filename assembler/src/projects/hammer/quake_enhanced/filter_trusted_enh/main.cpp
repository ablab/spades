//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <stdint.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <unordered_map>
#include "utils/logger/logger.hpp"

using std::string;
using std::unordered_map;

namespace {
/**
 * @variable Length of string buffer which will store k-mer.
 */
const uint32_t kMaxK = 100;
/**
 * @variable Every kStep k-mer will appear in the log.
 */
const int kStep = 1e5;

DECL_LOGGER("filter_trusted_enh")

struct Options {
  string ifile;
  string ofile;
  string badfile;
  string limits;
  bool valid;
  Options()
      : ifile(""),
        ofile(""),
        badfile(""),
        limits(""),
        valid(true) {}
};

void PrintHelp(char *progname) {
  printf("Usage: %s ifile.[q]cst ofile.trust ofile.bad file.limits\n", progname);
  printf("Where:\n");
  printf("\tifile.[q]cst\tfile with k|q-mer statistics\n");
  printf("\tofile.trust\ta filename where filtered data will be outputted\n");
  printf("\tofile.bud\ta filename where filtered garbage will be outputted\n");
  printf("\tfile.limits\tfile with q-value limits for k-mers\n");
}

Options ParseOptions(int argc, char *argv[]) {
  Options ret;
  if (argc != 5) {
    ret.valid = false;
  } else {
    ret.ifile = argv[1];
    ret.ofile = argv[2];
    ret.badfile = argv[3];
    ret.limits = argv[4];
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
  BasicConfigurator::configure();
  INFO(logger, "Starting filter_trusted: evaluating "
               << opts.ifile << ".");
  FILE *ifile = fopen(opts.ifile.c_str(), "r");
  FILE *ofile = fopen(opts.ofile.c_str(), "w");
  FILE *badfile = fopen(opts.badfile.c_str(), "w");
  FILE *limits_file = fopen(opts.limits.c_str(), "r");
  unordered_map<uint32_t, long double> limits;
  uint32_t x;
  long double limit;
  while (fscanf(limits_file, "%u %Lf", &x, &limit) == 2) {
    limits[x] = limit;
  }
  char kmer[kMaxK];
  char format[20];
  float freq = -1;
  int count;
  float q_count;
  snprintf(format, sizeof(format), "%%%ds%%d%%f%%f", kMaxK);
  uint64_t read_number = 0;
  while (fscanf(ifile, format, kmer, &count, &q_count, &freq) != EOF) {
    ++read_number;
    if (read_number % kStep == 0) {
      INFO(logger, "Reading k-mer " << read_number << ".");
    }
    if (q_count / count > limits[count]) {
      fprintf(ofile, "%s %d %f %f\n", kmer, count, q_count, freq);
    } else {
      fprintf(badfile, "%s %d %f %f\n", kmer, count, q_count, freq);
    }
  }
  return 0;
}
