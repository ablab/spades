#include <stdint.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
using log4cxx::LoggerPtr;
using log4cxx::Logger;
using log4cxx::BasicConfigurator;

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

LoggerPtr logger(Logger::getLogger("filter_trusted"));

struct Options {
  string ifile;
  string ofile;
  string badfile;
  float threshold;
  bool valid;
  Options()
      : ifile(""),
        ofile(""),
        badfile(""),
        threshold(-1),
        valid(true) {}  
};

void PrintHelp(char *progname) {
  printf("Usage: %s ifile.[q]cst ofile.trust ofile.bad threshold\n", progname);
  printf("Where:\n");
  printf("\tifile.[q]cst\tfile with k|q-mer statistics\n");
  printf("\tofile.trust\ta filename where filtered data will be outputted\n");
  printf("\tofile.bud\ta filename where filtered garbage will be outputted\n");
  printf("\tthreshold\tq-mer threshold\n");
}

Options ParseOptions(int argc, char *argv[]) {
  Options ret;
  if (argc != 5) {
    ret.valid = false;
  } else {
    ret.ifile = argv[1];
    ret.ofile = argv[2];
    ret.badfile = argv[3];
    ret.threshold = atof(argv[4]);
    if (ret.threshold <= -1e-5) {
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
  BasicConfigurator::configure();
  LOG4CXX_INFO(logger, "Starting filter_trusted: evaluating "
               << opts.ifile << ".");
  FILE *ifile = fopen(opts.ifile.c_str(), "r");
  FILE *ofile = fopen(opts.ofile.c_str(), "w");
  FILE *badfile = fopen(opts.badfile.c_str(), "w");
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
    } else {
      fprintf(badfile, "%s\n", kmer);
    }
  }
  return 0;
}
