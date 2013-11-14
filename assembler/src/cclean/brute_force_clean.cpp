// brute_force_clean.cpp Copyright (c) 10.11.2013 Kuprashevich Maksim
// This file is under MIT license.
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
// documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
// Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS
// OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

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
using additional::WORK_MODE_TYPE;
using additional::NONE;
using additional::SIMPLE;
using additional::BRUTE_SIMPLE;
using additional::BRUTE_DEEP;

const double MATCH_SCORE = 0.6;
const double MISMATCH_SCORE = 10;

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

bool BruteForceClean::operator()(const Read &read)
{
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
  const vector<string> &adapters = adap_gen_.GetSeqs();

  if (brute_ == BRUTE_SIMPLE) {
    // For each adapter align read and adapter
    for (auto adapt_string: adapters) {
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
          print_alignment(output_stream_, alignment,
                          seq_string, adapt_string, read_name, db_name_);
          print_bed(bed_stream_, read_name, alignment.ref_begin,
                    alignment.ref_end);
        }
      }

    }
  }

  if(brute_ == BRUTE_DEEP) {
    double align_score;
    for (auto adapt_string: adapters) {
      aligner.Align(adapt_string.c_str(), filter, &alignment);

      align_score = GetScoreWithQuality(alignment, read.getQuality());

      if (align_score > threshold_
          && is_alignment_good(alignment, seq_string,
                               adapt_string, aligned_part_fraction_)) {
#       pragma omp critical
        {
          cuted_ += 1;
          print_alignment(output_stream_, alignment,
                          seq_string, adapt_string, read_name, db_name_);
          print_bed(bed_stream_, read_name, alignment.ref_begin,
                    alignment.ref_end);
        }
      }

    }
  }

  return false;
}

double BruteForceClean::GetScoreWithQuality(const StripedSmithWaterman::Alignment &a,
                                            const Quality &qual)
{ // Try to get more realistic align score depend on read quality
  // Mathes and mismatches get from cigar alignment string below
  double score = 0.0;
  int ref_pos = 0, query_pos = 0;
  for (std::vector<uint32_t>::const_iterator it = a.cigar.begin();
       it != a.cigar.end(); ++it) {

    int num = (*it & 0xFFFFFFF0) >> 4;
    int op_code = *it & 0x0000000F;

    switch (op_code) {
      case 0: { //match
        for (int i = 0; i < num; ++i, ++ref_pos, ++query_pos)
          score += MATCH_SCORE;
        break;
      }
      case 1: { //insert
        for (int i = 0; i < num; ++i, ++query_pos)
          score += qual[query_pos] / MISMATCH_SCORE;
        break;
      }
      case 2: { //del
        for (int i = 0; i < num; ++i, ++ref_pos)
          score += qual[query_pos] / MISMATCH_SCORE;
        break;
      }
      default:
        break;
    }
  }

  return score;
}
