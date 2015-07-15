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
#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include <set>
#include <unordered_map>
#include <vector>
#include "logging.hpp"
#include "io/ireadstream.hpp"
#include "io/read.hpp"
#include "sequence/seq.hpp"
#include "valid_kmer_generator.hpp"
#define SUPPRESS_UNUSED(X) ((void) (X))

using std::string;
using std::set;
using std::vector;
using std::unordered_map;
using std::map;

namespace {

DECL_LOGGER("count")

struct KMerInfo {
  int count;
  double q_count;
  double q_inversed_count;
};

const uint32_t kK = 55;
typedef Seq<kK> KMer;
typedef unordered_map<KMer, KMerInfo, KMer::hash> UnorderedMap;

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
  bool valid;
  Options()
      : qvoffset(0),
        ifile(""),
        ofile(""),
        error_threshold(0),
        file_number(3),
        valid(true) {}
};

void PrintHelp(char *program_name) {
  printf("Usage: %s qvoffset ifile.fastq ofile.[q]cst file_number error_threshold\n",
         program_name);
  printf("Where:\n");
  printf("\tqvoffset\tan offset of fastq quality data\n");
  printf("\tifile.fastq\tan input file with reads in fastq format\n");
  printf("\tofile.[q]cst\ta filename where k-mer statistics will be outputted\n");
  printf("\terror_threshold\tnucliotides with quality lower then threshold will be cut from the ends of reads\n");
  printf("\tfile_number\thow many files will be used when splitting k-mers\n");
}

Options ParseOptions(int argc, char *argv[]) {
  Options ret;
  if (argc != 6) {
    ret.valid =  false;
  } else {
    ret.qvoffset = atoi(argv[1]);
    ret.valid &= (ret.qvoffset >= 0 && ret.qvoffset <= 255);
    ret.ifile = argv[2];
    ret.ofile = argv[3];
    ret.error_threshold = atoi(argv[4]);
    ret.valid &= (ret.error_threshold >= 0 && ret.error_threshold <= 255);
    ret.file_number = atoi(argv[5]);
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
void SplitToFiles(ireadstream ifs, const vector<FILE*> &ofiles,
                  uint8_t error_threshold) {
  uint32_t file_number = ofiles.size();
  uint64_t read_number = 0;
  while (!ifs.eof()) {
    ++read_number;
    if (read_number % kStep == 0) {
      INFO("Reading read " << read_number << ".");
    }
    Read r;
    ifs >> r;
    KMer::hash hash_function;
    for (ValidKMerGenerator<kK> gen(r, error_threshold); gen.HasMore(); gen.Next()) {
      KMer kmer = gen.kmer();
      if (KMer::less2()(!kmer, kmer)) {
        kmer = !kmer;
      }
      FILE *cur_file = ofiles[hash_function(kmer) % file_number];
      KMer::BinWrite(cur_file, kmer);
      double correct_probability = gen.correct_probability();
      fwrite(&correct_probability, sizeof(correct_probability), 1, cur_file);
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
void EvalFile(FILE *ifile, FILE *ofile) {
  UnorderedMap stat_map;
  char buffer[kK + 1];
  buffer[kK] = 0;
  KMer kmer;
  while (KMer::BinRead(ifile, &kmer)) {
    KMerInfo &info = stat_map[kmer];
    double correct_probability = -1;
    bool readed =
        fread(&correct_probability, sizeof(correct_probability),
              1, ifile);
    assert(readed == 1);
    SUPPRESS_UNUSED(readed);
    double inversed_probability = 1 / correct_probability;
    // ToDo 0.5 threshold ==>> command line option
    if (correct_probability < 0.5) {
      inversed_probability = 0;
    }
    info.q_count += correct_probability;
    info.count += 1;
    info.q_inversed_count += inversed_probability;
  }
  for (UnorderedMap::iterator it = stat_map.begin();
       it != stat_map.end(); ++it) {
    const KMerInfo &info = it->second;
    fprintf(ofile, "%s %d %f %f\n", it->first.str().c_str(),
            info.count, info.q_count, info.q_inversed_count);
  }
}

void run(const Options &opts) {
  INFO("Starting preproc: evaluating "
       << opts.ifile << ".");
  vector<FILE*> ofiles(opts.file_number);
  for (uint32_t i = 0; i < opts.file_number; ++i) {
    char filename[50];
    snprintf(filename, sizeof(filename), "%u.kmer.part", i);
    ofiles[i] = fopen(filename, "wb");
    assert(ofiles[i] != NULL && "Too many files to open");
  }
  SplitToFiles(ireadstream(opts.ifile, opts.qvoffset),
               ofiles, opts.error_threshold);
  for (uint32_t i = 0; i < opts.file_number; ++i) {
    fclose(ofiles[i]);
  }
  FILE *ofile = fopen(opts.ofile.c_str(), "w");
  assert(ofile != NULL && "Too many files to open");
  for (uint32_t i = 0; i < opts.file_number; ++i) {
    char ifile_name[50];
    snprintf(ifile_name, sizeof(ifile_name), "%u.kmer.part", i);
    FILE *ifile = fopen(ifile_name, "rb");
    INFO("Processing " << ifile_name << ".");
    EvalFile(ifile, ofile);
    INFO("Processed " << ifile_name << ".");
    fclose(ifile);
  }
  fclose(ofile);
  INFO("Preprocessing done. You can find results in " <<
       opts.ofile << ".");
}
}

int main(int argc, char *argv[]) {
  Options opts = ParseOptions(argc, argv);
  if (!opts.valid) {
    PrintHelp(argv[0]);
    return 1;
  }
  run(opts);
  return 0;
}
