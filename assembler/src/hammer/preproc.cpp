/**
 * @file    preproc.cpp
 * @author  Alex Davydow
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * @section DESCRIPTION
 *
 * For each k-mer this programm calculates number of occuring in
 * the reads provided. Reads file is supposed to be in fastq 
 * format.
 */
#include <omp.h>
#include <cstdlib>
#include <string>
#include <set>
#include <utility>
#include <vector>
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#include "common/read/ireadstream.hpp"
#include "common/read/read.hpp"
#include "hammer/defs.hpp"
#include "hammer/hammer_config.hpp"
#include "hammer/kmer_functions.hpp"
#include "hammer/kmer_part_joiner.hpp"

using std::make_pair;
using std::pair;
using std::string;
using std::set;
using std::vector;
using log4cxx::LoggerPtr;
using log4cxx::Logger;
using log4cxx::BasicConfigurator;

namespace {

LoggerPtr logger(Logger::getLogger("preproc"));
const int kStep = 1e5;

struct Options {
  uint32_t qvoffset;
  string ifile;
  string ofile;
  uint32_t nthreads;
  uint32_t read_batch_size;
  uint32_t file_number;
  bool valid;
  Options()
      : qvoffset(0),
        ifile(""),
        ofile(""),
        nthreads(1),
        read_batch_size(1e6),
        file_number(2),
        valid(true) {}
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
  }
  return ret;
}

void SplitToFiles(const string &ifile, size_t qvoffset, size_t file_number) {
  ireadstream ifs(ifile.c_str(), qvoffset);
  vector<FILE*> files(file_number);
  for (uint32_t i = 0; i < file_number; ++i) {
    char filename[50];
    snprintf(filename, sizeof(filename), "%u.kmer.part", i);
    files[i] = fopen(filename, "w");
  }
  size_t read_number = 0;
  while (!ifs.eof()) {
    // reading a batch of reads
    ++read_number;
    if (read_number % kStep == 0) {
      LOG4CXX_INFO(logger, "Reading read " << read_number << ".");
    }
    Read r;
    ifs >> r;
    if (TrimBadQuality(&r) >= K) {
      vector<KMer> kmers = GetKMers(r);
      KMer::hash hash_function;
      for (size_t i = 0; i < kmers.size(); ++i) {
        int file_id = hash_function(kmers[i]) % file_number;
        fprintf(files[file_id], "%s\n", kmers[i].str().c_str());
      }
    }
  }
  for (size_t i = 0; i < file_number; ++i) {
    fclose(files[i]);
  }
  ifs.close();
  LOG4CXX_INFO(logger, "Reads wroten to separate files.");
}

void EvalFile(FILE *ifile, FILE *ofile) {
  char buffer[K + 1];
  KMerStatMap stat_map;
  while (fscanf(ifile, "%s", buffer) != EOF) {
#pragma message("Warning about uninitialized _M_instance looks like a fake")
    // Next line produces a misterious warning saying that _M_instance is
    // undeifileed
    // Looks like it is in some way connected to the line
    //       int file_id = hash_function(kmers[i]) % file_number;
    // line in SplitToFiles
    KMer kmer(buffer);
    ++stat_map[kmer].count;
  }
  for (KMerStatMap::iterator it = stat_map.begin();
       it != stat_map.end();
       ++it) {
    fprintf(ofile,
            "%s %u\n", it->first.str().c_str(),
            (unsigned int)it->second.count);
  }
}

void MergeAndSort(const vector<FILE*> &ifiles, FILE *ofile) {
  KMerPartJoiner joiner(ifiles);
  while (!joiner.IsEmpty()) {
    pair<string, int> kmer_stat = joiner.Next();
    fprintf(ofile, "%s %d\n", kmer_stat.first.c_str(), kmer_stat.second);
  }
}
}

int main(int argc, char * argv[]) {
  Options opts = ParseOptions(argc, argv);
  if (!opts.valid) {
    PrintHelp();
    return 1;
  }
  BasicConfigurator::configure();
  LOG4CXX_INFO(logger, "Starting preproc: evaluating " << opts.ifile <<
               " in " << opts.nthreads << " threads.");
  SplitToFiles(opts.ifile, opts.qvoffset, opts.file_number);

  for (uint32_t i = 0; i < opts.file_number; ++i) {
    char ifile_name[50];
    char ofile_name[50];
    snprintf(ifile_name, sizeof(ifile_name), "%u.kmer.part", i);
    snprintf(ofile_name, sizeof(ofile_name), "%u.result.part", i);
    FILE *ifile = fopen(ifile_name, "r");
    FILE *ofile = fopen(ofile_name, "w");
    LOG4CXX_INFO(logger, "Processing " << ifile_name << ".");
    EvalFile(ifile, ofile);
    LOG4CXX_INFO(logger, "Processed " << ifile_name << ". " <<
                 "You can find the result in " << ofile_name <<
                 ".");
    fclose(ifile);
    fclose(ofile);
  }
  vector<FILE*> ifiles;
  for (uint32_t i = 0; i < opts.file_number; ++i) {
    char ifile_name[50];
    snprintf(ifile_name, sizeof(ifile_name), "%u.result.part", i);
    FILE *ifile = fopen(ifile_name, "r");
    ifiles.push_back(ifile);
  }
  FILE *ofile = fopen(opts.ofile.c_str(), "w");

  LOG4CXX_INFO(logger, "Starting message.");
  MergeAndSort(ifiles, ofile);
  for (uint32_t i = 0; i < opts.file_number; ++i) {
    char ifile_name[50];
    snprintf(ifile_name, sizeof(ifile_name), "%u.result.part", i);
    FILE *ifile = fopen(ifile_name, "r");
    fclose(ifile);
  }
  fclose(ofile);
  LOG4CXX_INFO(logger,
               "Preprocessing done. You can find results in " <<
               opts.ofile << ".");
  return 0;
}
