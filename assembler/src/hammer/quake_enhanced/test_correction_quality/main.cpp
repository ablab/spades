#include <stdint.h>
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
/**
 * @variable Every kStep k-mer will appear in the log.
 */
const int kStep = 1e5;

struct Options {
  string genom_file;
  string trust_file;
  string bad_file;
  float threshold;
  bool valid;
  Options()
      : genom_file(""),
        trust_file(""),
        bad_file(""),
        valid(true) {}  
};

void PrintHelp(char *progname) {
  printf("Usage: %s genom.[q]cst ifile.trust ifile.bad\n", progname);
  printf("Where:\n");
  printf("\tgenom.[q]cst\tfile with k|q-mer statistics from real genom\n");
  printf("\tifile.trust\ta filename where filtered data is\n");
  printf("\tifile.bud\ta filename where filtered garbage is\n");
}

Options ParseOptions(int argc, char *argv[]) {
  Options ret;
  if (argc != 4) {
    ret.valid = false;
  } else {
    ret.genom_file = argv[1];
    ret.trust_file = argv[2];
    ret.bad_file = argv[3];
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
  snprintf(format, sizeof(format), "%%%ds", kMaxK);
  while (fscanf(trust_file, format, kmer) != EOF) {
    if (real_kmers.count(string(kmer)) > 0) {
      ++trusted;
    } else {
      ++trusted_fail;
    }
  }
  printf("trusted: %d\n", trusted + trusted_fail);
  printf("erroneous: %d\n", trusted_fail);
  snprintf(format, sizeof(format), "%%%ds", kMaxK);
  while (fscanf(bad_file, format, kmer) != EOF) {
    if (real_kmers.count(string(kmer)) > 0) {
      ++bad_fail;
    } else {
      ++bad;
    }
  }
  printf("bad: %d\n", trusted + trusted_fail);
  printf("erroneous: %d\n", trusted_fail);
  return 0;
}
