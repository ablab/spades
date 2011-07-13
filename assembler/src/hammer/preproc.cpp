/*
 * preproc.cpp
 *
 *  Created on: 11.05.2011
 *      Author: snikolenko
 */
#include "hammer_config.hpp" 
#include <omp.h>
#include <cmath>
#include <string>
#include <cstdlib>
#include <vector>
#include <utility>
#include "hammer/defs.hpp"
#include "common/read/read.hpp"
#include "common/read/ireadstream.hpp"
#include "common/sequence/seq.hpp"

using std::cout;
using std::string;
using std::vector;
using std::endl;
using std::pair;

int qvoffset;

/**
 * add k-mers from read to map
 */
void AddKMers(const Read & r, KMerStatMap & v) {
  KMerStatMap::iterator it;
  string s = r.getSequenceString();
  int i = 0;
  while (true) {
    i = r.firstValidKmer(i, K);
    if (i == -1) break;
    KMer kmer = KMer(r.getSubSequence(i, K));
    while (true) {
      it = v.find(kmer);
      if (it != v.end()) {
	it->second.count++;
      } else {
	pair<KMer, KMerStat> p;
	p.first = kmer;
	p.second.count = 1; 
	v.insert(p);
      }
      if (i + K < (int)r.size() && is_nucl(s[i + K])) {
	kmer = kmer << r[i + K];
	++i;
      } else {
	i += K;
	break;
      }
    }
  }
}

const int kReadBatchSize = 1e6;

void PrintHelp() {
  printf("Usage: ./preproc qvoffset ifile.fastq ofile.kmer [nthreads]\n");
  printf("Where:\n");
  printf("\tqvoffset\tan offset of fastq quality data\n");
  printf("\tifile.fastq\tan input file with reads in fastq format\n");
  printf("\tofile.kmer\ta filename where kmer statistics will be outputed\n");
  printf("\tnthreads\ta number of threads (one by default)\n");
}

int main(int argc, char * argv[]) {
  // Read command line arguments. Exit and print help if anything goes wrong.
  if (argc != 5 && argc != 4) {
    PrintHelp();
    return 1;
  }
  qvoffset = atoi(argv[1]);  
  if (qvoffset < 0 || qvoffset > 255) {
    PrintHelp();
    return 1;
  }
  string reads_filename = argv[2];
  string kmer_filename = argv[3];
  int nthreads = 1;
  if (argc == 5) {
    nthreads = atoi(argv[4]);
  }
  if (nthreads <= 0) {
    PrintHelp();
    return 1;
  }
  //

  printf("Starting preproc: evaluating %s in %d threads.\n", reads_filename.c_str(), nthreads);

  ireadstream ifs(reads_filename.data(), qvoffset);
  size_t batch_number = 0;
  vector<KMerStatMap> vv(nthreads);     
  while (!ifs.eof()) {
    // reading a batch of reads
    ++batch_number;
    printf("Reading batch %d.\n", (unsigned int)batch_number);
    vector<Read> rv;
    for (int read_number = 0; read_number < kReadBatchSize; ++read_number) {
      Read r;      
      ifs >> r; 
      // trim the reads for bad quality and process only the ones with at least K "reasonable" elements
      if (r.trimBadQuality() >= K) {
	rv.push_back(r);
      }
      if (ifs.eof()) {
	break;
      }
    }    

    printf("Batch %u read.\n", (unsigned int)batch_number);
    // ToDo: add multithreading (map and reduce)
    for(int i = 0; i < (int)rv.size(); ++i) {
      AddKMers(rv[i], vv[0]);
      AddKMers(!(rv[i]), vv[0]);
    }
    printf("Batch %u added.\n", (unsigned int)batch_number);
  }
  ifs.close();
  printf("All k-mers added to maps.\n");
  for (int i = 0; i < (int)vv.size(); ++i) {
    printf("size(%d) = %u\n", i, (unsigned int)vv[i].size());
  }
  
  FILE* f = fopen(kmer_filename.data(), "w");
  for (KMerStatMap::iterator it = vv[0].begin(); it != vv[0].end(); ++it) {
    fprintf(f, "%s %u\n", it->first.str().data(), (unsigned int) it->second.count);
  }
  fclose(f);
  return 0;
}


