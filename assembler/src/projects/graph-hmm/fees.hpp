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

  bool check_i_loop(size_t i) const;
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
