//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef QUAKE_ENHANCED_HPP_
#define QUAKE_ENHANCED_HPP_

#include <string>
#include <sstream>

namespace quake_enhanced {
class Options {
 public:
  std::string read_file;
  std::string corrected_read_file;
  std::string help_message;
  // Count options.
  //
  std::string kmer_count_file;
  std::string hash_file_prefix;
  int hash_file_number;
  int quality_offset;
  // nucleotides with quality lower than threshold will be cut from
  // the ends of the read.
  int quality_threshold;
  // PrepareHist options.
  //
  std::string hist_file;
  // CorrectHist options.
  //
  std::string trusted_hist_file;
  std::string bad_hist_file;
  // we will look for maximum which is at least top_threshold times
  // higher than previous minimum
  size_t top_threshold; 
  // trying to find Gauss's average we will go to the left and to the
  // right until we rich coverage average_min * max
  double average_min;
  // GenerateLimits options.
  std::string limits_file;
  // we will consider k-mer trusted if its probability of being
  // correct is at least bad_threshold times higher then its
  // probability of being bad. 
  double bad_threshold;
  // FilterTrusted options.
  std::string trusted_kmer_file;
  std::string bad_kmer_file;

  bool valid;
  Options(int argc, char **argv);
 private:
  std::ostringstream help_builder;
  void Validate();
};
}

#endif  // QUAKE_ENHANCED_HPP_
