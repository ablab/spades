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
#include <omp.h>
#include <cstdlib>
#include <string>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
#include <map>
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#include "common/read/ireadstream.hpp"
#include "common/read/strobe_read.hpp"
#include "common/read/read.hpp"
#include "hammer/defs.hpp"
#include "hammer/kmer_freq_info.hpp"
#include "hammer/kmer_part_joiner.hpp"
#include "hammer/valid_kmer_generator.hpp"

using std::make_pair;
using std::pair;
using std::string;
using std::set;
using std::vector;
using std::unordered_map;
using std::map;
using log4cxx::LoggerPtr;
using log4cxx::Logger;
using log4cxx::BasicConfigurator;

namespace {

const uint32_t kK = 2;
typedef Seq<kK> KMer;
typedef RCReaderWrapper<ireadstream, Read> ReadStream;
typedef map<KMer, KMerFreqInfo, KMer::less2> OrderedMap;
typedef unordered_map<KMer, KMerFreqInfo, KMer::hash> UnorderedMap;

LoggerPtr logger(Logger::getLogger("preproc"));
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
  uint32_t nthreads;
  /**
   * @variable How many files will be used when splitting k-mers.
   */
  uint32_t file_number;
  bool sorted;
  /**
   * @variable If options provided are valid.
   */
  bool valid;
  Options()
      : qvoffset(0),
        ifile(""),
        ofile(""),
        nthreads(1),
        file_number(3),
        sorted(false),
        valid(true) {}
};

void PrintHelp() {
  printf("Usage: ./preproc qvoffset ifile.fastq ofile.kmer sorted nthreads\n");
  printf("Where:\n");
  printf("\tqvoffset\tan offset of fastq quality data\n");
  printf("\tifile.fastq\tan input file with reads in fastq format\n");
  printf("\tofile.kmer\ta filename where k-mer statistics will be outputted\n");
  printf("\tfile_number\thow many files will be used when splitting k-mers\n");
  printf("\tsorted\t\t'y' if you need sorting, 'n' otherwise\n");
  printf("\tnthreads\ta number of threads (one by default)\n");
}

Options ParseOptions(int argc, char * argv[]) {
  Options ret;
  if (argc != 7) {
    ret.valid =  false;
  } else {
    ret.qvoffset = atoi(argv[1]);
    ret.valid &= (ret.qvoffset >= 0 && ret.qvoffset <= 255);
    ret.ifile = argv[2];
    ret.ofile = argv[3];
    ret.file_number = atoi(argv[4]);
    string sorted(argv[5]);
    if (sorted == "y") {
      ret.sorted = true;
    } else if (sorted == "n") {
      ret.sorted = false;
    } else {
      ret.valid = false;
    }
    ret.nthreads = atoi(argv[6]);
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
void SplitToFiles(ReadStream ifs, const vector<FILE*> &ofiles) {
  uint32_t file_number = ofiles.size();
  uint32_t read_number = 0;
  while (!ifs.eof()) {
    ++read_number;
    if (read_number % kStep == 0) {
      LOG4CXX_INFO(logger, "Reading read " << read_number << ".");
    }
    Read r;
    ifs >> r;
    KMer::hash hash_function;
    for (ValidKMerGenerator<kK> gen(r); gen.HasMore(); gen.Next()) {
      int file_id = hash_function(gen.kmer()) % file_number;
      float p = gen.correct_probability();
      if (p != p) {
        for (int i = 0; i < kK; ++i) {
          printf("%d ", (int)r.getQualityString()[i + gen.pos()]);
        }
        printf(": %d\n", gen.pos());
      }
      fprintf(ofiles[file_id], "%s %f\n", gen.kmer().str().c_str(), 
              static_cast<float>(gen.correct_probability()));
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
void EvalFile(FILE *ifile, FILE *ofile) {
  KMerStatMap stat_map;
  char format[10];
  snprintf(format, sizeof(format), "%%%ds%%f", kK);
  char buffer[kK + 1];
  float correct_probability;
  while (fscanf(ifile, format, buffer, &correct_probability) != EOF) {
    KMer kmer(buffer);
    KMerFreqInfo &info = stat_map[kmer];
    ++info.count;
    info.q_count += correct_probability;
  }
  for (typename KMerStatMap::iterator it = stat_map.begin();
       it != stat_map.end();
       ++it) {
    fprintf(ofile,
            "%s %u %f\n", it->first.str().c_str(),
            it->second.count, it->second.q_count);
  }
}

/**
 * Given a set of sorted files with k-mers and their frequency this
 * function merge them into one sorted file.
 * @param ifiles Files to merge.
 * @param ofile Output file.
 */
void MergeAndSort(const vector<FILE*> &ifiles, FILE *ofile) {
  KMerPartJoiner joiner(ifiles, kK);
  while (!joiner.IsEmpty()) {
    KMerFreqInfo info = joiner.Next();
    fprintf(ofile, "%s %d %f\n", info.kmer.c_str(), info.count, info.q_count);
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
  {
    vector<FILE*> ofiles(opts.file_number);
    for (uint32_t i = 0; i < opts.file_number; ++i) {
      char filename[50];
      snprintf(filename, sizeof(filename), "%u.kmer.part", i);
      ofiles[i] = fopen(filename, "w");
    }
    ireadstream ir(opts.ifile, opts.qvoffset);
    SplitToFiles(ReadStream(ir), ofiles);
    for (uint32_t i = 0; i < opts.file_number; ++i) {
      fclose(ofiles[i]);
    }
  }
  if (opts.sorted) {
    LOG4CXX_INFO(logger, "Reads written to separate files.");
#pragma omp parallel for num_threads(opts.nthreads)
    for (uint32_t i = 0; i < opts.file_number; ++i) {
      char ifile_name[50];
      char ofile_name[50];
      snprintf(ifile_name, sizeof(ifile_name), "%u.kmer.part", i);
      snprintf(ofile_name, sizeof(ofile_name), "%u.result.part", i);
      FILE *ifile = fopen(ifile_name, "r");
      FILE *ofile = fopen(ofile_name, "w");
      LOG4CXX_INFO(logger, "Processing " << ifile_name << ".");
      EvalFile<OrderedMap>(ifile, ofile);
      LOG4CXX_INFO(logger, "Processed " << ifile_name << ". " <<
                   "You can find the result in " << ofile_name <<
                   ".");
      fclose(ifile);
      fclose(ofile);
    }

    LOG4CXX_INFO(logger, "Starting merge.");
      vector<FILE*> ifiles;
      for (uint32_t i = 0; i < opts.file_number; ++i) {
        char ifile_name[50];
        snprintf(ifile_name, sizeof(ifile_name), "%u.result.part", i);
        FILE *ifile = fopen(ifile_name, "r");
      ifiles.push_back(ifile);
      }
      FILE *ofile = fopen(opts.ofile.c_str(), "w");
      MergeAndSort(ifiles, ofile);
      for (uint32_t i = 0; i < opts.file_number; ++i) {
        fclose(ifiles[i]);
      }
      fclose(ofile);
  } else {
    FILE *ofile = fopen(opts.ofile.c_str(), "w");
    for (uint32_t i = 0; i < opts.file_number; ++i) {
      char ifile_name[50];
      snprintf(ifile_name, sizeof(ifile_name), "%u.kmer.part", i);
      FILE *ifile = fopen(ifile_name, "r");
      LOG4CXX_INFO(logger, "Processing " << ifile_name << ".");
      EvalFile<UnorderedMap>(ifile, ofile);
      LOG4CXX_INFO(logger, "Processed " << ifile_name << ".");
      fclose(ifile);
    }
    fclose(ofile);
  }
  LOG4CXX_INFO(logger,
               "Preprocessing done. You can find results in " <<
               opts.ofile << ".");
  return 0;
}
