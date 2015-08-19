//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef UTILS_HPP
#define UTILS_HPP

#include <ssw/ssw_cpp.h> // Striped Smith-Waterman aligner
#include <io/read.hpp>
#include "additional.cpp"
#include "running_modes.hpp"
#include "adapter_index.hpp"

namespace cclean_utils {

std::string ReverseComplement(const std::string& read);

std::unordered_map<std::string, std::string> ProcessArgs(int argc, char *argv[],
                                                         bool *ok, std::string *error);

double GetScoreWithQuality(const StripedSmithWaterman::Alignment &a,
                                            const Quality &qual);

inline bool is_alignment_good(const StripedSmithWaterman::Alignment& a,
                              const std::string& sequence,
                              const std::string& query,
                              double aligned_part_fraction) {
  //  Сheck that query adjoins or even overlaps the sequence edge
  return (std::min(a.query_end - a.query_begin + 1, a.ref_end - a.ref_begin + 1)
         / (double) query.size() > aligned_part_fraction) /*&&
         (a.ref_begin == 0 || a.ref_end == sequence.size() - 1)*/;
}

// Cut read from start to end position of best aligment with adapter
Read CutRead(const Read &r, int start_pos, int end_pos);
void RestoreFromCigar(const std::string& ref, const std::string& query,
                      std::string& out_ref, std::string& out_query,
                      const StripedSmithWaterman::Alignment& a);

inline double GetMismatches(const std::string &read, const std::string &adapter,
                         const StripedSmithWaterman::Alignment &a)  {
  std::string aligned_read;
  std::string aligned_adapter;
  RestoreFromCigar(read, adapter, aligned_read, aligned_adapter, a);
  int size = (int)std::min(aligned_read.length(), aligned_adapter.length());
  int mismatched_score = 0;
  for (int i = 0; i < size; ++i)  {
    if (aligned_read[i] != aligned_adapter[i])
      ++mismatched_score;
  }
  return static_cast<double>(mismatched_score);
}
// end of namespace
}
#endif /* UTILS_HPP */
