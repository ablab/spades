//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "brute_force_clean.hpp"

#include <string>
#include <vector>
#include <iostream>

#include "adapter_index.hpp"
#include <ssw/ssw_cpp.h> // Striped Smith-Waterman aligner
#include "additional.cpp"
#include "output.hpp"

using std::string;
using std::vector;
using StripedSmithWaterman::Filter;
using StripedSmithWaterman::Aligner;
using StripedSmithWaterman::Alignment;
using cclean_output::print_alignment;
using cclean_output::print_bad;
using cclean_output::print_match;
using cclean_output::print_read;

static inline bool is_alignment_good(const StripedSmithWaterman::Alignment& a,
                              const std::string& sequence,
                              const std::string& query,
                              double aligned_part_fraction) {
  // Сheck that query adjoins or even overlaps the sequence edge
  return (std::min(a.query_end - a.query_begin + 1, a.ref_end - a.ref_begin + 1)
         / (double) query.size() > aligned_part_fraction) &&
         (a.ref_begin == 0 || a.ref_end == sequence.size() - 1);
}

Read BruteForceClean::operator()(const Read &read, bool *ok) {
  const string &read_name = read.getName();
  const string &seq_string = read.getSequenceString();
  Filter filter; // SSW filter
  Aligner aligner; // SSW aligner
  aligner.SetReferenceSequence(seq_string.c_str(),
                               static_cast<int>(seq_string.size()));
  Alignment alignment;

  //  It can be many alignment adaps, so we searching the most probable
  double best_score;
  if (mode_ == BRUTE_SIMPLE)  // so in both mode first overlap will initialize as best
    best_score = mismatch_threshold_;
  if (mode_ == BRUTE_WITH_Q)
    best_score = score_threshold_;
  std::string best_adapter = "";

  //  For each adapter align read and adapter
  for (std::string adapt_string: adap_seqs_) {

    aligner.Align(adapt_string.c_str(), filter, &alignment);
    if((*checker)(read, alignment, aligned_part_fraction_, adapt_string,
                  &best_score)) {
      best_adapter = adapt_string;
    }
  }

  if (!best_adapter.empty())  {
      aligner.Align(best_adapter.c_str(), filter, &alignment);
      aligned_ += 1;
      Read cuted_read = cclean_utils::CutRead(read, alignment.ref_begin,
                                              alignment.ref_end);
      if (full_inform_)  // If user want full output
#       pragma omp critical
        print_alignment(aligned_output_stream_, alignment, seq_string,
                        best_adapter, read_name, db_name_);

      // Cuted read must be >= minimum lenght specified by arg
      if (cuted_read.getSequenceString().size() >= read_mlen_) {
        if (full_inform_)  // If user want full output
#         pragma omp critical
          print_bad(bad_stream_, read_name, alignment.ref_begin, alignment.ref_end);
        (*ok) = true;
        return cuted_read;
      }
      else {
        if (full_inform_)
#         pragma omp critical
          print_bad(bad_stream_, read_name, 0, alignment.ref_end);
        (*ok) = false;
        return cuted_read;
      }
    }
  else {
    // Read was not aligned with any adapter
    (*ok) = true;
    return read;
  }
}
