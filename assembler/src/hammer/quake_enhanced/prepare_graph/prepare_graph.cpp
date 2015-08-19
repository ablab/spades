//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <stdint.h>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <string>
#include "logging.hpp"

using std::map;
using std::string;

namespace {

DECL_LOGGER("count")

/**
 * @variable Length of string buffer which will store k-mer.
 */
const uint32_t kMaxK = 100;
/**
 * @variable Every kStep k-mer will appear in the log.
 */
const int kStep = 1e5;

typedef map<uint64_t, uint32_t> Map;
struct Options {
  string ifile;
  string ofile;
  bool binary;
  bool qmer;
  int ticks_per_step;
  bool valid;
  Options()
      : ifile(""),
        ofile(""),
        binary(false),
        qmer(false),
        ticks_per_step(1),
        valid(true) {}
};

void PrintHelp() {
  printf("Usage: ./prepare_graph ifile.[q]cst ofile.kmer ticks_per_step [b]\n");
  printf("Where:\n");
  printf("\tifile.[q]cst\tfile woth k|q-mer statistics\n");
  printf("\tofile.kmer\ta filename where data prepared for plotting will be outputted\n");
  printf("\tticks_per_step\tone over horizontal scale\n");
  printf("\tb\t\tif ifile.[q]cst is a binary file instead of text[currently not working]\n");
}

Options ParseOptions(int argc, char *argv[]) {
  Options ret;
  if (argc != 4 && argc != 5) {
    ret.valid = false;
  } else {
    ret.ifile = argv[1];
    ret.ofile = argv[2];
    ret.ticks_per_step = atoi(argv[3]);
    if (ret.ticks_per_step <= 0) {
      ret.valid = false;
    }
    if (argc == 5) {
      if (string(argv[4]) == "b") {
        ret.binary = true;
      } else {
        ret.valid = false;
      }
    }
  }
  return ret;
}

void run(const Options &opts) {
  INFO("Starting prepare_graph: evaluating "
       << opts.ifile << ".");
  FILE *ifile = fopen(opts.ifile.c_str(), "r");
  FILE *ofile = fopen(opts.ofile.c_str(), "w");
  Map freq_to_num;
  char kmer[kMaxK];
  char format[20];
  snprintf(format, sizeof(format), "%%%ds%%d%%f%%f", kMaxK);
  float freq = -1;
  int count;
  float q_count;
  uint64_t read_number = 0;
  while (fscanf(ifile, format, kmer, &count, &q_count, &freq) != EOF) {
    ++read_number;
    if (read_number % kStep == 0) {
      INFO("Reading k-mer " << read_number << ".");
    }
    ++freq_to_num[
        static_cast<uint64_t>(freq * opts.ticks_per_step + 0.5)];
  }
  for (Map::iterator it = freq_to_num.begin();
       it != freq_to_num.end(); ++it) {
    fprintf(ofile, "%f %d\n",
            static_cast<float>(it->first) / opts.ticks_per_step,
            it->second);
  }
  fclose(ofile);
  fclose(ifile);
}
}


int main(int argc, char *argv[]) {
  Options opts = ParseOptions(argc, argv);
  if (!opts.valid) {
    PrintHelp();
    return 1;
  }
  run(opts);
  return 0;
}
