//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <stdint.h>
#include <string.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <unordered_set>

using std::unordered_set;
using std::string;

namespace {
/**
 * @variable Length of string buffer which will store k-mer.
 */
const uint32_t kMaxK = 100;

struct Options {
  string genom_file;
  string trust_file;
  string bad_file;
  bool full;
  float threshold;
  bool valid;
  Options()
      : genom_file(""),
        trust_file(""),
        bad_file(""),
    full(false),
        valid(true) {}  
};

void PrintHelp(char *progname) {
  printf("Usage: %s genom.[q]cst ifile.trust ifile.bad [--full]\n", progname);
  printf("Where:\n");
  printf("\tgenom.[q]cst\tfile with k|q-mer statistics from real genom\n");
  printf("\tifile.trust\ta filename where filtered data is\n");
  printf("\tifile.bud\ta filename where filtered garbage is\n");
  printf("\t--full\tpass this option to output all incorrect k-mers with their names to stdout\n");
}

Options ParseOptions(int argc, char *argv[]) {
  Options ret;
  if (argc < 4 || argc > 5) {
    ret.valid = false;
  } else {
    ret.genom_file = argv[1];
    ret.trust_file = argv[2];
    ret.bad_file = argv[3];
    if (argc == 5 && ( !strcmp(argv[4], "--full") || !strcmp(argv[4], "-f") ) )
      ret.full = true;
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
  FILE *genom_file = fopen(opts.genom_file.c_str(), "r");
  FILE *trust_file = fopen(opts.trust_file.c_str(), "r");
  FILE *bad_file = fopen(opts.bad_file.c_str(), "r");
  char kmer[kMaxK];
  char format[20];
  float freq = -1;
  int count;
  float q_count;
  snprintf(format, sizeof(format), "%%%ds%%d%%f%%f", kMaxK);
  unordered_set<string> real_kmers;
  while (fscanf(genom_file, format, kmer, &count, &q_count, &freq) != EOF) {
    real_kmers.insert(string(kmer));
  }
  int trusted = 0;
  int trusted_fail = 0;
  int bad = 0;
  int bad_fail = 0;
  while (fscanf(trust_file, format, kmer, &count, &q_count, &freq) != EOF) {
    if (real_kmers.count(string(kmer)) > 0) {
      ++trusted;
    } else {
      ++trusted_fail;
      if ( opts.full ) printf("  %s\t%d\t%f\t%f\n", kmer, count, q_count, freq);
    }
  }
  printf("trusted: %d\n", trusted + trusted_fail);
  printf("erroneous: %d\n", trusted_fail);
  while (fscanf(bad_file, format, kmer, &count, &q_count, &freq) != EOF) {
    if (real_kmers.count(string(kmer)) > 0) {
      ++bad_fail;
      if ( opts.full ) printf("  %s\t%d\t%f\t%f\n", kmer, count, q_count, freq);
    } else {
      ++bad;
    }
  }
  printf("bad: %d\n", bad + bad_fail);
  printf("erroneous: %d\n", bad_fail);
  return 0;
}
