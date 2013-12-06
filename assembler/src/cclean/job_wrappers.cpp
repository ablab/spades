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
  // Сheck that query adjoins or even overlaps the sequence edge
  return (std::min(a.query_end - a.query_begin + 1, a.ref_end - a.ref_begin + 1)
         / (double) query.size() > aligned_part_fraction) &&
         (a.ref_begin == 0 || a.ref_end == sequence.size() - 1);
}

bool SimpleClean::operator()(const Read &r) {
  try {
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

    // Try to align the artifacts for corresponding kmers
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    aligner.SetReferenceSequence(sequence.c_str(), sequence.size());
    for (auto it = to_check.begin(), et = to_check.end(); it != et; ++it) {
      StripedSmithWaterman::Alignment alignment;
      const std::string& query = index_.seq(*it);
      aligner.Align(query.c_str(), filter, &alignment);

      if (alignment.mismatches < mismatch_threshold_ &&
          is_alignment_good(alignment, sequence, query,
                            aligned_part_fraction_)) {
#       pragma omp critical
        {
          aligned_ += 1;
          print_alignment(aligned_output_, alignment, sequence, query, name, db_);
          print_bed(bed_, name, alignment.ref_begin, alignment.ref_end);
          Read cuted_read = cclean_utils::CutRead(r, alignment.ref_begin,
                                                  alignment.ref_end);
          print_read(output_stream_, cuted_read);
        }
      }
    }

  } catch (std::exception& e) {
    ERROR(e.what());
  }

  return false;
}
