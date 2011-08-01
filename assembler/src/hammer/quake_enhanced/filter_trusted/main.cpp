#include <stdint.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
using log4cxx::LoggerPtr;
using log4cxx::Logger;
using log4cxx::BasicConfigurator;

using std::map;
using std::string;

namespace {
/**
 * @variable Length of string buffer which will store k-mer.
 */
const uint32_t kMaxK = 100;
/**
 * @variable Every kStep k-mer will appear in the log.
 */
const int kStep = 1e5;

LoggerPtr logger(Logger::getLogger("prepare_graph"));

typedef map<uint64_t, uint32_t> Map;
struct Options {
  float threshold;
  string ifile;
  string ofile;
  Options()
      : ifile(""),
        ofile(""),
        threshold(-1);
};

void PrintHelp(char *progname) {
  printf("Usage: %s ifile.[q]cst ofile.trust threshold\n", progname);
  printf("Where:\n");
  printf("\tifile.[q]cst\tfile with k|q-mer statistics\n");
  printf("\tofile.trus\ta filename where filtered data will be outputted\n");
  printf("\tthreshold\tq-mer threshold\n");
}

Options ParseOptions(int argc, char *argv[]) {
  Options ret;
  if (argc != 4) {
    ret.valid = false;
  } else {
    ret.ifile = argv[1];
    ret.ofile = argv[2];
    ret.threshold = atof(argv[3]);
    if (ret.ticks_per_step <= 0) {
      ret.valid = false;
    }
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
  FILE *ifile = fopen(opts.ifile.c_str(), "r");
  FILE *ofile = fopen(opts.ofile.c_str(), "w");
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
      LOG4CXX_INFO(logger, "Reading k-mer " << read_number << ".");
    }
    if (q_count > opts.threshold) {
      fprintf(ofile, "%s\n", kmer);
    }
  }
  return 0;
}
