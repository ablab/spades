//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <set>

#include "job_wrappers.hpp"
#include "utils/logger/log_writers.hpp"
#include "adapter_index.hpp"
#include "valid_kmer_generator.hpp"
#include "adapter_index.hpp"
#include "output.hpp"
#include "ssw/ssw_cpp.h"
#include "utils.hpp"

using cclean_output::print_alignment;
using cclean_output::print_bad;
using cclean_output::print_match;
using cclean_output::print_read;

Read SimpleClean::operator()(const Read &read, bool *ok)
{
  const std::string& name = read.getName();
  const std::string& sequence = read.getSequenceString();

  std::set<size_t> to_check;
  ValidKMerGenerator<cclean::K> gen(sequence.c_str(), NULL, sequence.size());
  while (gen.HasMore()) {
    cclean::KMer kmer = gen.kmer();

    auto it = index_.find(kmer);
    if (it != index_.end())
      to_check.insert(it->second.begin(), it->second.end());

    gen.Next();
  }

  //  Try to align the artifacts for corresponding kmers
  StripedSmithWaterman::Aligner aligner;
  StripedSmithWaterman::Filter filter;
  StripedSmithWaterman::Alignment alignment; //  why it was in for loop?
  aligner.SetReferenceSequence(sequence.c_str(), sequence.size());

  //  Pointer on best match adapter
  const std::string *best_adapter = nullptr;
  double best_score;
  if (mode_ == SINGLE_END)  // so in both mode first overlap will initialize as best
    best_score = mismatch_threshold_;
  if (mode_ == SINGLE_END_Q)
    best_score = score_threshold_;
  best_adapter = nullptr;

  for (auto it = to_check.begin(), et = to_check.end(); it != et; ++it) {
    const std::string &query = index_.seq(*it);
    aligner.Align(query.c_str(), filter, &alignment);
    // Check is this apapter better then previous best
    if((*checker)(read, alignment, aligned_part_fraction_, query,
                  &best_score)) {
      best_adapter = &query;
    }
  }

  if (best_adapter != nullptr)  {
      aligner.Align(best_adapter->c_str(), filter, &alignment);
      aligned_ += 1;
      Read cuted_read = cclean_utils::CutRead(read, alignment.ref_begin,
                                              alignment.ref_end);
      if (full_inform_)  // If user want full output
#       pragma omp critical
        print_alignment(aligned_output_stream_, alignment, sequence,
                        *best_adapter,name, db_name_);

      // Cuted read must be >= minimum lenght specified by arg
      if (cuted_read.getSequenceString().size() >= read_mlen_) {
        if (full_inform_)
#         pragma omp critical
          print_bad(bad_stream_, name, alignment.ref_begin, alignment.ref_end);
        (*ok) = true;
        return cuted_read;
      }
      else {
        if (full_inform_)
#         pragma omp critical
          print_bad(bad_stream_, name, 0, alignment.ref_end);
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
