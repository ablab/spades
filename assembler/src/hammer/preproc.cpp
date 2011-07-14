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

struct Options {
  int qvoffset;
  int k;
  string ifile;
  string ofile;
  int nthreads;
  int read_batch_size;
  bool good;
  Options() : k(K), nthreads(1), read_batch_size(1e6), good(true) {}
};


void PrintHelp() {
  printf("Usage: ./preproc qvoffset ifile.fastq ofile.kmer [nthreads]\n");
  printf("Where:\n");
  printf("\tqvoffset\tan offset of fastq quality data\n");
  printf("\tifile.fastq\tan input file with reads in fastq format\n");
  printf("\tofile.kmer\ta filename where k-mer statistics will be outputed\n");
  printf("\tnthreads\ta number of threads (one by default)\n");
}
  
Options ParseOptions(int argc, char * argv[]) {
  Options ret;
  if (argc != 5 && argc != 4) {
    ret.good =  false;
  } else {
    ret.qvoffset = atoi(argv[1]);  
    ret.good &= (ret.qvoffset >= 0 && ret.qvoffset <= 255);
    //    ret.k = atoi(argv[2]);
    //    ret.good &= (ret.k > 0);
    ret.ifile = argv[2];
    ret.ofile = argv[3];
    if (argc == 5) {
      ret.nthreads = atoi(argv[4]);
    }
    ret.good &= ret.nthreads > 0;
  }
  return ret;
}

}

int main(int argc, char * argv[]) {
  Options opts = ParseOptions(argc, argv);
  if (!opts.good) {
    PrintHelp();
    return 1;
  }
  printf("Starting preproc: evaluating %s in %d threads.\n", opts.ifile.c_str(), opts.nthreads);

  ireadstream ifs(opts.ifile.c_str(), opts.qvoffset);
  size_t batch_number = 0;
  vector<KMerStatMap> vv(opts.nthreads);     
  while (!ifs.eof()) {
    // reading a batch of reads
    ++batch_number;
    printf("Reading batch %d.\n", (unsigned int)batch_number);
    vector<Read> rv;
    for (int read_number = 0; read_number < opts.read_batch_size; ++read_number) {
      Read r;      
      ifs >> r; 
      // trim the reads for bad quality and process only the ones with at least K "reasonable" elements
      if (TrimBadQuality(r) >= K) {
	rv.push_back(r);
      }
      if (ifs.eof()) {
	break;
      }
    }    

    printf("Batch %u read.\n", (unsigned int)batch_number);
    // ToDo: add multithreading (map and reduce)
    for(int i = 0; i < (int)rv.size(); ++i) {
      AddKMers(rv[i], vv[0], opts.k);
      AddKMers(!(rv[i]), vv[0], opts.k);
    }
    printf("Batch %u added.\n", (unsigned int)batch_number);
  }
  ifs.close();
  printf("All k-mers added to maps.\n");
  for (int i = 0; i < (int)vv.size(); ++i) {
    printf("size(%d) = %u\n", i, (unsigned int)vv[i].size());
  }
  
  FILE* f = fopen(opts.ofile.c_str(), "w");
  for (KMerStatMap::iterator it = vv[0].begin(); it != vv[0].end(); ++it) {
    fprintf(f, "%s %u\n", it->first.str().data(), (unsigned int) it->second.count);
  }
  fclose(f);
  return 0;
}


