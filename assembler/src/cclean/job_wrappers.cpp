#include <set>

#include "job_wrappers.hpp"
#include "logger/log_writers.hpp"
#include "adapter_index.hpp"
#include "valid_kmer_generator.hpp"
#include "adapter_index.hpp"
#include "output.hpp"
#include "ssw/ssw_cpp.h"
#include "utils.hpp"

using cclean_output::print_alignment;
using cclean_output::print_bed;
using cclean_output::print_match;
using cclean_output::print_read;

static inline bool is_alignment_good(const StripedSmithWaterman::Alignment& a,
                              const std::string& sequence,
                              const std::string& query,
                              double aligned_part_fraction) {
  //  Сheck that query adjoins or even overlaps the sequence edge
  return (std::min(a.query_end - a.query_begin + 1, a.ref_end - a.ref_begin + 1)
         / (double) query.size() > aligned_part_fraction) &&
         (a.ref_begin == 0 || a.ref_end == sequence.size() - 1);
}

bool SimpleClean::operator()(const Read &r) {
  if(mode_ == additional::SINGLE_END) SingleEndClean(r);
  if(mode_ == additional::SINGLE_END_Q) SingleEndWQualityClean(r);

  return false;
}

void SimpleClean::SingleEndClean(const Read &r)
{
  const std::string& name = r.getName();
  const std::string& sequence = r.getSequenceString();

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
  StripedSmithWaterman::Alignment alignment; // why it was in for loop?
  aligner.SetReferenceSequence(sequence.c_str(), sequence.size());
  //  It can be aligned many adapters, but we want to find most probable
  int best_mismatch = mismatch_threshold_;
  int cur_mismatch = INT_MAX;
  const std::string *best_adapter = nullptr;

  for (auto it = to_check.begin(), et = to_check.end(); it != et; ++it) {
    const std::string& query = index_.seq(*it);
    aligner.Align(query.c_str(), filter, &alignment);
    cur_mismatch = cclean_utils::GetMismatches(r.getSequenceString(),
                                                 query, alignment);
    if (cur_mismatch < best_mismatch &&
        is_alignment_good(alignment, sequence, query,
                          aligned_part_fraction_)) {
      best_mismatch = cur_mismatch;
      best_adapter = &query;
    }
  }

  if (best_adapter != nullptr)  {
#   pragma omp critical
    {
      aligner.Align(best_adapter->c_str(), filter, &alignment);
      aligned_ += 1;
      if (options_["inform"] == "FULL") { // if user want full output
        print_alignment(aligned_output_, alignment, sequence, *best_adapter,
                        name, db_);
        print_bed(bed_, name, alignment.ref_begin, alignment.ref_end);
      }
      Read cuted_read = cclean_utils::CutRead(r, alignment.ref_begin,
                                              alignment.ref_end);
      if (cuted_read.getSequenceString().size() >= read_mlen_) {
        // printing cuted read only if it >= minimum lenght specified by arg
        print_read(output_stream_, cuted_read);
      }
      else { // 0 - because start of read
        if (options_["inform"] == "FULL") print_bed(bed_, name, 0, alignment.ref_end);
      }
    }
  }
  else {
#   pragma omp critical
    print_read(output_stream_, r);
  }
}

void SimpleClean::SingleEndWQualityClean(const Read &r)
{
  const std::string& name = r.getName();
  const std::string& sequence = r.getSequenceString();

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
  StripedSmithWaterman::Alignment alignment; // why it was in for loop?
  aligner.SetReferenceSequence(sequence.c_str(), sequence.size());
  //  It can be aligned many adapters, but we want to find most probable
  const std::string *best_adapter = nullptr;
  double align_score;
  double best_align_score = mismatch_threshold_; // not lower then theshold

  for (auto it = to_check.begin(), et = to_check.end(); it != et; ++it) {
    const std::string& query = index_.seq(*it);
    aligner.Align(query.c_str(), filter, &alignment);
    align_score = cclean_utils::GetScoreWithQuality(alignment, r.getQuality().str());
    if (align_score >= best_align_score &&
        is_alignment_good(alignment, sequence, query,
                          aligned_part_fraction_)) {
      best_align_score = align_score;
      best_adapter = &query;
    }
  }

  if (best_adapter != nullptr)  {
#   pragma omp critical
    {
      aligner.Align(best_adapter->c_str(), filter, &alignment);
      aligned_ += 1;
      if (options_["inform"] == "FULL") { // if user want full output
        print_alignment(aligned_output_, alignment, sequence, *best_adapter,
                        name, db_);
        print_bed(bed_, name, alignment.ref_begin, alignment.ref_end);
      }
      Read cuted_read = cclean_utils::CutRead(r, alignment.ref_begin,
                                              alignment.ref_end);
      if (cuted_read.getSequenceString().size() >= read_mlen_) {
        // printing cuted read only if it >= minimum lenght specified by arg
        print_read(output_stream_, cuted_read);
      }
      else { // 0 - because start of read
        if (options_["inform"] == "FULL") print_bed(bed_, name, 0, alignment.ref_end);
      }
    }
  }
  else {
#   pragma omp critical
    print_read(output_stream_, r);
  }
}
