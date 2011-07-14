/*
 * preproc.cpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */
#include "hammer_config.hpp" 
#include <omp.h>
#include <string>
#include <cstdlib>
#include <vector>
#include "hammer/defs.hpp"
#include "hammer/kmer_functions.hpp"
#include "common/read/read.hpp"
#include "common/read/ireadstream.hpp"

using std::string;
using std::vector;

namespace {

const int kStep = 1e5;

struct Options {
  int qvoffset;
  string ifile;
  string ofile;
  int nthreads;
  int read_batch_size;
  int file_number;
  bool valid;
  Options() : nthreads(1), read_batch_size(1e6), file_number(100), valid(true) {}
};


void PrintHelp() {
  printf("Usage: ./preproc qvoffset k ifile.fastq ofile.kmer [nthreads]\n");
  printf("Where:\n");
  printf("\tqvoffset\tan offset of fastq quality data\n");
  printf("\tifile.fastq\tan input file with reads in fastq format\n");
  printf("\tofile.kmer\ta filename where k-mer statistics will be outputed\n");
  printf("\tnthreads\ta number of threads (one by default)\n");
}
  
Options ParseOptions(int argc, char * argv[]) {
  Options ret;
  if (argc != 4 && argc != 5) {
    ret.valid =  false;
  } else {
    ret.qvoffset = atoi(argv[1]);  
    ret.valid &= (ret.qvoffset >= 0 && ret.qvoffset <= 255);
    ret.ifile = argv[2];
    ret.ofile = argv[3];
    if (argc == 5) {
      ret.nthreads = atoi(argv[4]);
    }
    ret.valid &= ret.nthreads > 0;
  }
  return ret;
}

}

int main(int argc, char * argv[]) {

  Options opts = ParseOptions(argc, argv);
  if (!opts.valid) {
    PrintHelp();
    return 1;
  }
  printf("Starting preproc: evaluating %s in %d threads.\n", opts.ifile.c_str(), opts.nthreads);

  ireadstream ifs(opts.ifile.c_str(), opts.qvoffset);
  vector<KMerStatMap> vv(opts.nthreads);     
  vector<FILE*> files(opts.file_number);
  for (int i = 0; i < opts.file_number; ++i) {
    files[i] = fopen(((string)itoa(i)) +  ".kmer.part", "w");
  }
  size_t read_number = 0;
  while (!ifs.eof()) {
    // reading a batch of reads
    ++read_number;
    if (read_number % kStep == 0) {
      printf("Reading read %d.\n", (unsigned int)read_number);
    }
    Read r;      
    ifs >> r; 
    vector<KMer> kmers = GetKMers(r);
    fprintf(files[r], "%s", r.str().c_str());
  }    
  printf("Reads wroten to separate files.\n");

  for (int i = 0; i < opts.file_number; ++i) {
    fclose(files[i]);
  }

  return 0;
}


