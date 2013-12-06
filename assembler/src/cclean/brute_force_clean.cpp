#include "brute_force_clean.hpp"

#include <string>
#include <vector>
#include <iostream>

#include "adapter_index.hpp"
#include <ssw/ssw_cpp.h> // Striped Smith-Waterman aligner
#include "additional.cpp"
#include "output.hpp"
#include "utils.hpp"

using std::string;
using std::vector;
using StripedSmithWaterman::Filter;
using StripedSmithWaterman::Aligner;
using StripedSmithWaterman::Alignment;
using additional::WorkModeType;
using additional::NONE;
using additional::SIMPLE;
using additional::BRUTE_SIMPLE;
using additional::BRUTE_DEEP;
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


int BruteForceClean::cuted_ = 0;

bool BruteForceClean::operator()(const Read &read) {
  if (brute_ == BRUTE_SIMPLE)  BruteSimple(read);
  if (brute_ == BRUTE_DEEP)    BruteDeep(read);

  return false;
}

void BruteForceClean::BruteSimple(const Read &read) {
  const string &read_name = read.getName();
  const string &seq_string = read.getSequenceString();
//  std::cout << read.getName() << std::endl;
//  std::cout << read.getPhredQualityString(33) << std::endl;
//  std::cin.ignore().get();
  Filter filter; // SSW filter
  Aligner aligner; // SSW aligner
  aligner.SetReferenceSequence(seq_string.c_str(),
                               static_cast<int>(seq_string.size()));
  Alignment alignment; // Struct that will contain align result
  // Simple mode: Best align score <= threeshold in config
  // Deep mode: Every match give +0.6 score, every mismatch -Q/10
  // >= threeshold in config

  // For each adapter align read and adapter
  for (auto adapt_string: adap_seqs_) {
    // If score >= then THEESHOLD, then cut adapter from read
    // If result is empty save read in BAD_READS_FILE
    // Save read in GOOD_READS_FILE
    aligner.Align(adapt_string.c_str(), filter, &alignment);

    if (alignment.mismatches < threshold_
        && is_alignment_good(alignment, seq_string,
                             adapt_string, aligned_part_fraction_)) {
#       pragma omp critical
      {
        cuted_ += 1;
        print_alignment(aligned_output_stream_, alignment,
                        seq_string, adapt_string, read_name, db_name_);
        print_bed(bed_stream_, read_name, alignment.ref_begin,
                  alignment.ref_end);
        Read cuted_read = cclean_utils::CutRead(read, alignment.ref_begin,
                                                alignment.ref_end);
        print_read(output_stream_, cuted_read);
      }
    }

  }
}

void BruteForceClean::BruteDeep(const Read &read) {
  const string &read_name = read.getName();
  const string &seq_string = read.getSequenceString();

  Filter filter; // SSW filter
  Aligner aligner; // SSW aligner
  aligner.SetReferenceSequence(seq_string.c_str(),
                               static_cast<int>(seq_string.size()));
  Alignment alignment; // Struct that will contain align result
  // Simple mode: Best align score <= threeshold in config
  // Deep mode: Every match give +0.6 score, every mismatch -Q/10
  // >= threeshold in config

  double align_score;
  for (auto adapt_string: adap_seqs_) {
    aligner.Align(adapt_string.c_str(), filter, &alignment);

    align_score = cclean_utils::GetScoreWithQuality(alignment, read.getQuality());

    if (align_score > threshold_
        && is_alignment_good(alignment, seq_string,
                             adapt_string, aligned_part_fraction_)) {
#       pragma omp critical
      {
        cuted_ += 1;
        print_alignment(aligned_output_stream_, alignment,
                        seq_string, adapt_string, read_name, db_name_);
        print_bed(bed_stream_, read_name, alignment.ref_begin,
                  alignment.ref_end);
        Read cuted_read = cclean_utils::CutRead(read, alignment.ref_begin,
                                                alignment.ref_end);
        print_read(output_stream_, cuted_read);
      }
    }

  }
}
