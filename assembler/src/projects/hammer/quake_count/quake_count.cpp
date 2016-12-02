//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

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
 * For each k-mer this program calculates number of occurring in
 * the reads provided. Reads file is supposed to be in fastq
 * format.
 */

#include "standard.hpp"

#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include <set>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <iomanip>
#include "io/ireadstream.hpp"
#include "io/read.hpp"
#include "sequence/seq.hpp"
#include "kmer_freq_info.hpp"
#include "valid_kmer_generator.hpp"
#define SUPPRESS_UNUSED(X) ((void) (X))

using std::string;
using std::set;
using std::vector;
using std::unordered_map;
using std::map;
using std::ofstream;
using std::ifstream;
using std::pair;

namespace {

const uint32_t kK = 55;
typedef Seq<kK> KMer;
typedef unordered_map<KMer, KMerFreqInfo, KMer::hash> UnorderedMap;

void print_time() {
    time_t rawtime;
    tm * ptm;
    time ( &rawtime );
    ptm = gmtime( &rawtime );
    std::cout << std::setfill('0') << "[ " << std::setw(2) << ptm->tm_hour << ":" << std::setw(2) << ptm->tm_min
            << ":" << std::setw(2) << ptm->tm_sec << " ] ";
}

#define LOG(a) print_time(); std::cout << a << std::endl

/**
 * @variable Every kStep k-mer will appear in the log.
 */
const int kStep = 1e5;

struct Options {
  /**
   * @variable An offset for quality in a fastq file.
   */
  uint32_t qvoffset;
  string ifile;
  string ofile;
  uint32_t error_threshold;
  /**
   * @variable How many files will be used when splitting k-mers.
   */
  uint32_t file_number;
  bool q_mers;
  bool valid;
  Options()
      : qvoffset(0),
        ifile(""),
        ofile(""),
        error_threshold(0),
        file_number(3),
        q_mers(false),
        valid(true) {}
};

void PrintHelp(char *program_name) {
  printf("Usage: %s qvoffset ifile.fastq ofile.[q]cst file_number error_threshold [q]\n",
         program_name);
  printf("Where:\n");
  printf("\tqvoffset\tan offset of fastq quality data\n");
  printf("\tifile.fastq\tan input file with reads in fastq format\n");
  printf("\tofile.[q]cst\ta filename where k-mer statistics will be outputted\n");
  printf("\terror_threshold\tnucliotides with quality lower then threshold will be cut from the ends of reads\n");
  printf("\tfile_number\thow many files will be used when splitting k-mers\n");
  printf("\tq\t\tif you want to count q-mers instead of k-mers.\n");
}

Options ParseOptions(int argc, char *argv[]) {
  Options ret;
  if (argc != 6 && argc != 7) {
    ret.valid =  false;
  } else {
    ret.qvoffset = atoi(argv[1]);
    ret.valid &= (ret.qvoffset >= 0 && ret.qvoffset <= 255);
    ret.ifile = argv[2];
    ret.ofile = argv[3];
    ret.error_threshold = atoi(argv[4]);
    ret.valid &= (ret.error_threshold >= 0 && ret.error_threshold <= 255);
    ret.file_number = atoi(argv[5]);
    if (argc == 7) {
      if (string(argv[6]) == "q") {
        ret.q_mers = true;
      } else {
        ret.valid = false;
      }
    }
  }
  return ret;
}

/**
 * This function reads reads from the stream and splits them into
 * k-mers. Then k-mers are written to several file almost
 * uniformly. It is guaranteed that the same k-mers are written to the
 * same files.
 * @param ifs Steam to read reads from.
 * @param ofiles Files to write the result k-mers. They are written
 * one per line.
 */
void SplitToFiles(ireadstream ifs, vector<ofstream *> &ofiles,
                  bool q_mers, uint8_t error_threshold) {
  uint32_t file_number = ofiles.size();
  uint64_t read_number = 0;
  while (!ifs.eof()) {
    ++read_number;
    if (read_number % kStep == 0) {
      LOG("Reading read " << read_number << ".");
    }
    Read r;
    ifs >> r;
    KMer::hash hash_function;
    for (ValidKMerGenerator<kK> gen(r, error_threshold); gen.HasMore(); gen.Next()) {
      KMer kmer = gen.kmer();
      if (KMer::less2()(!kmer, kmer)) {
        kmer = !kmer;
      }
      ofstream &cur_file = *ofiles[hash_function(kmer) % file_number];
      KMer::BinWrite(cur_file, kmer);
      if (q_mers) {
        double correct_probability = gen.correct_probability();
        cur_file.write((const char*) &correct_probability, sizeof(correct_probability));
      }
    }
  }
}

/**
 * This function reads k-mer and calculates number of occurrences for
 * each of them.
 * @param ifile File with k-mer to process. One per line.
 * @param ofile Output file. For each unique k-mer there will be a
 * line with k-mer itself and number of its occurrences.
 */
template<typename KMerStatMap>
void EvalFile(ifstream &ifile, FILE *ofile, bool q_mers) {
  KMerStatMap stat_map;
  char buffer[kK + 1];
  buffer[kK] = 0;
  KMer kmer;
  while  (KMer::BinRead(ifile, &kmer)) {
    KMerFreqInfo &info = stat_map[kmer];
    if (q_mers) {
      double correct_probability = -1;
      ifile.read((char *) &correct_probability, sizeof(correct_probability));
      assert(ifile.fail());
      info.q_count += correct_probability;
    } else {
      info.count += 1;
    }
  }
  for (typename KMerStatMap::iterator it = stat_map.begin();
       it != stat_map.end(); ++it) {
    fprintf(ofile, "%s ", it->first.str().c_str());
    if (q_mers) {
      fprintf(ofile, "%f\n", it->second.q_count);
    } else {
      fprintf(ofile, "%d\n", it->second.count);
    }
  }
}
}

int main(int argc, char *argv[]) {
  Options opts = ParseOptions(argc, argv);
  if (!opts.valid) {
    PrintHelp(argv[0]);
    return 1;
  }
  // BasicConfigurator::configure();
  LOG("Starting preproc: evaluating " << opts.ifile << ".");
  vector<ofstream*> ofiles(opts.file_number);
  for (uint32_t i = 0; i < opts.file_number; ++i) {
    char filename[50];
    snprintf(filename, sizeof(filename), "%u.kmer.part", i);
    ofiles[i] = new ofstream(filename);
    assert(!ofiles[i]->fail() && "Too many files to open");
  }
  SplitToFiles(ireadstream(opts.ifile, opts.qvoffset),
               ofiles, opts.q_mers, opts.error_threshold);
  for (uint32_t i = 0; i < opts.file_number; ++i) {
    delete ofiles[i];
  }
  FILE *ofile = fopen(opts.ofile.c_str(), "w");
  assert(ofile != NULL && "Too many files to open");
  for (uint32_t i = 0; i < opts.file_number; ++i) {
    char ifile_name[50];
    snprintf(ifile_name, sizeof(ifile_name), "%u.kmer.part", i);
    ifstream ifile(ifile_name);
    LOG("Processing " << ifile_name << ".");
    EvalFile<UnorderedMap>(ifile, ofile, opts.q_mers);
    LOG("Processed " << ifile_name << ".");
  }
  fclose(ofile);
  LOG("Preprocessing done. You can find results in " << opts.ofile << ".");
  return 0;
}
