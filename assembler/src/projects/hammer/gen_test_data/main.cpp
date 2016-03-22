//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
using namespace std;
struct Options {
  size_t len;
  size_t k;
  size_t min_quality;
  size_t coverage;
  size_t seed;
  bool valid;
};

void PrintHelp(string prog_name) {
  cout << "Usage: " << prog_name << " len k mean_quality coverage seed" << endl;
}

Options ParseOptions(int argc, char **argv) {
  Options ret;
  ret.valid = false;
  if (argc != 5) {
    return ret;
  }
  ret.len = atoi(argv[0]);
  ret.k = atoi(argv[1]);
  if (ret.k > ret.len) {
    return ret;
  }
  ret.min_quality = atof(argv[2]);
  ret.coverage = atoi(argv[3]);
  ret.seed = atoi(argv[4]);
  ret.valid = true;
  return ret;
}

int main(int argc, char **argv) {
  Options opt = ParseOptions(argc - 1, argv + 1);
  if (!opt.valid) {
    PrintHelp(argv[0]);
    return 0; 
  }
  srand(opt.seed);
  string acgt = "ACGT";
  string ref_genome(opt.len, ' ');
  for (size_t i = 0; i < opt.len; ++i) {
    ref_genome[i] = acgt[rand() % acgt.size()];
  }
  ofstream ref_file("reference.fasta");
  ofstream read_file("reads.fastq");
  ref_file << "> Random genome generated with len: " << opt.len 
           << " seed: " << opt.seed << endl;
  ref_file << ref_genome << endl;
  read_file << "; Reads for random genome generated with len: " << opt.len 
            << " seed: " << opt.seed << endl;
  size_t n = opt.len * opt.coverage / opt.k;
  for (size_t i = 0 ; i < n; ++i) {
    int s = rand() % (opt.len - opt.k);
    read_file << "@start:" << s + 1 << " length=" << opt.k << endl;
    string read(opt.k, ' ');
    string qual(opt.k, 'B');
    for (size_t i = 0; i < opt.k; ++i) {
      qual[i] = 33 + opt.min_quality + rand() % (41 - opt.min_quality);
      double erp = pow(10.0, - qual[i] - 33 / 10.0);
      if (rand() < erp * RAND_MAX) {
        read[i] = acgt[rand() % acgt.size()];
      } else {
        read[i] = ref_genome[i + s];
      }
    }
    read_file << read << endl;
    read_file << "+" << endl;
    read_file << qual << endl;
  }
  return 0;
}
