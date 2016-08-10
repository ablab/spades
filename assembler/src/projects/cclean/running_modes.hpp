//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef RUNNING_MODES_HPP
#define RUNNING_MODES_HPP

#include <unordered_map>
#include <string>
#include <iostream>
#include <iomanip>
#include "additional.cpp"
#include "adapter_index.hpp"

class Cleaner {

  public:
    static void ProcessDataset();
    // Correct reads in a given file
    static void CorrectReadFile(const cclean::AdapterIndex &index,
                                const std::string &fname,
                                std::ofstream *outf_good, std::ofstream *outf_bad,
                                WorkModeType mode);
    // Correct reads in a given pair of files
    static void CorrectPairedReadFiles(const cclean::AdapterIndex &index,
                                       const std::string &fnamel,
                                       const std::string &fnamer,
                                       std::ofstream *ofbadl,
                                       std::ofstream *ofcorl,
                                       std::ofstream *ofbadr,
                                       std::ofstream *ofcorr,
                                       std::ofstream *ofunp,
                                       WorkModeType mode);
    // Parallel correction of batch of reads
    static inline unsigned CorrectReadsBatch(AbstractCclean *cleaner,
                                             std::vector<bool> *results,
                                             std::vector<Read> *reads,
                                             size_t buf_size, unsigned nthreads) {
      unsigned bad = 0;
#     pragma omp parallel for shared(reads, results) num_threads(nthreads)
      for (size_t i = 0; i < buf_size; ++i) {
        bool ok;
        (*reads)[i] = (*cleaner)((*reads)[i], &ok);
        (*results)[i] = ok;
        if (!ok) ++bad;
      }
      return bad;
    }
    // Get pure file name without extension
    inline static std::string getPureFilename(const std::string &fname) {
      std::string tmp = path::filename(fname);
      std::string pure_file_name = "";
      size_t pos = tmp.find(".fastq");
      if (pos == std::string::npos)
        pure_file_name = tmp;
      else
        pure_file_name = tmp.substr(0, pos);
      return pure_file_name;
    }
    // Get filename for reads
    inline static std::string getReadsFilename(const std::string &dirprefix,
                                               const std::string &fname,
                                               const std::string &suffix) {
      const std::string &pure_file_name = getPureFilename(fname);
      return (dirprefix + "/" + pure_file_name + "." + suffix + ".fastq");
    }
    // Define mode depends on config file data
    inline static WorkModeType getMode() {
        WorkModeType mode;
        if (cfg::get().use_bruteforce) {
          if (cfg::get().use_quality) mode = BRUTE_WITH_Q;
          else                        mode = BRUTE_SIMPLE;
        }
        else {
          if (cfg::get().use_quality) mode = SINGLE_END_Q;
          else                        mode = SINGLE_END;
        }
        return mode;
    }
    // Create and return cleaner depends on mode
    inline static AbstractCclean* getCleaner(std::ofstream *outf_alig_debug,
                                             std::ofstream *outf_bad_deb,
                                             const std::string &db,
                                             WorkModeType mode, unsigned mlen,
                                             const cclean::AdapterIndex &index,
                                             bool deb_info);

};

#endif /* RUNNING_MODES_H_ */
