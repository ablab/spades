//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "stdint.h"
#include <string>
#include <vector>
#include "io/reader.hpp"
#include "io/ireadstream.hpp"

// ToDo => to settings.hpp
const uint32_t kK = 31;

namespace quake_enhanced {
class Quake {
 public:
  Quake() : cur_state_(kInitial) {}
  void Count(std::string ifile, std::string ofile,
             std::string hash_file_prefix, uint32_t hash_file_number, 
             uint8_t quality_offset, uint8_t quality_threshold);
  void PrepareHists(std::string hist_file, std::string trusted_hist_file,
                    std::string bad_hist_file, uint32_t top_threshold,
                    double average_min);
  void PrepareLimits(long double threshold, std::string limits_file);
  void FilterTrusted(std::string ifile, std::string ofile, std::string badfile);
 private:
  enum QuakeState {kInitial, kCountDone, kRealHistPrepared, kRealHistPrinted,
                   kTrustedHistPrepared, kLimitsCounted, kLimitsPrinted, 
                   kTrustedFiltered};
  QuakeState cur_state_;
  std::string kmer_count_file_;
  std::vector<uint32_t> trusted_hist_;
  std::vector<uint32_t> real_hist_;
  std::vector<double> limits_;

  // Count
  /**
   * This function reads reads from the stream and splits them into
   * k-mers. Then k-mers are written to several file almost
   * uniformly. It is guaranteed that the same k-mers are written to the
   * same files.
   * @param ifs Steam to read reads from.
   * @param ofiles Files to write the result k-mers. They are written
   * one per line.
   */
  void SplitToFiles(ireadstream ifs, vector<ofstream*> &ofiles,
                    uint8_t error_threshold);
  /**
   * This function reads k-mer and calculates number of occurrences for
   * each of them.
   * @param ifile File with k-mer to process. One per line.
   * @param ofile Output file. For each unique k-mer there will be a
   * line with k-mer itself and number of its occurrences.
   */
  void EvalFile(ifstream &ifile, ofstream &ofile);

  // PrepareHists
  void AddToHist(double freq);
  void PrepareRealHist();
  void PrintRealHist(std::string hist_file);
  void PrepareTrustedHist(std::string trusted_hist_file, 
                          std::string bad_hist_file, uint32_t top_threshold, 
                          double average_min);
  // PrepareLimits
  long double GetLimit(uint32_t x, long double threshold);
  void CountLimits(long double threshold);
  void PrintLimits(std::string limits_file);
  // Filter k-mers
};
}
