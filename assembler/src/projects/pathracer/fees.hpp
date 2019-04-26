#pragma once

#include "common/hmm/hmmer_fwd.h"

#include <string>
#include <iostream>
#include <vector>
#include <limits>
#include <cstdlib>

namespace hmm {

class DigitalCodind {
 public:
  DigitalCodind(const ESL_ALPHABET *abc);
  DigitalCodind();

  size_t operator()(char ch) const { return inmap_[ch]; }

 private:
  std::vector<size_t> inmap_;
  size_t k_;
};

struct Fees {
  size_t M;
  size_t k;
  std::vector<std::vector<double>> t;
  std::vector<std::vector<double>> mat;
  std::vector<std::vector<double>> ins;
  DigitalCodind code;
  std::string consensus;
  std::string name;

  struct {size_t l25 = 1000000, l100 = 50000, l500 = 10000; } state_limits;
  size_t max_insertion_length = 30;
  double absolute_threshold;

  double cleavage_cost;
  bool local = false;
  size_t minimal_match_length = 50;

  bool use_experimental_i_loop_processing = false;

  double empty_sequence_score() const;
  double all_matches_score(const std::string &seq) const;

  double all_matches_score() const {
    return all_matches_score(consensus);
  }

  bool is_proteomic() const {
    return k == 20;
  }

  bool is_nucleotide() const {
    return !is_proteomic();
  }

  Fees() {
    state_limits.l25 = 1000000;
    state_limits.l100 = 50000;
    state_limits.l500 = 10000;
    absolute_threshold = 250.0;
    cleavage_cost = 0;
  }

  bool check_i_loop(size_t i) const;
  size_t count_negative_loops() const;
  bool is_i_loop_non_negative(size_t i) const { return check_i_loop(i); }
  bool check_i_negative_loops() const;

  void reverse();
  Fees reversed() const {
    Fees copy{*this};
    copy.reverse();
    return copy;
  }
};

Fees levenshtein_fees(const std::string &s, double mismatch = 1, double gap_open = 1, double gap_ext = 1);
Fees fees_from_hmm(const P7_HMM *hmm, const ESL_ALPHABET *abc, double lambda = 0);
Fees fees_from_file(const std::string &filename);

}  // namespace hmm

// vim: set ts=2 sw=2 et :
