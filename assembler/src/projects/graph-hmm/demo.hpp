#pragma once

#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"

#include "graph.hpp"
extern "C" {
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
}

class DigitalCodind {
 public:
  DigitalCodind(const ESL_ALPHABET *abc)
      : inmap_(abc->inmap, abc->inmap + 128) {}
  DigitalCodind(const std::string &s)
      : inmap_(s.c_str(), s.c_str() + s.length()) {}
  DigitalCodind() : inmap_(128) {
    std::string alphabet = "ACGT";
    for (size_t i = 0; i < alphabet.length(); ++i) {
      inmap_[alphabet[i]] = i;
    }
  }

  size_t operator()(char ch) const { return inmap_[ch]; }

 private:
  std::vector<size_t> inmap_;
};

struct Fees {
  size_t M;
  std::vector<std::vector<double>> t;
  std::vector<std::vector<double>> mat;
  std::vector<std::vector<double>> ins;
  DigitalCodind code;

  bool check_i_loop(size_t i) const {
    return t[i][p7H_II] + *min_element(ins[i].cbegin(), ins[i].cend()) > 0;
  }

  bool check_i_negative_loops() const {
    for (size_t i = 0; i <= M; ++i) {
      if (!check_i_loop(i)) {
        std::cout << t[i][p7H_II] << " " << i << " "
                  << *min_element(ins[i].cbegin(), ins[i].cend()) << std::endl;
        std::cout << t[i][p7H_II] << " " << ins[i][0] << " " << ins[i][1] << " "
                  << ins[i][2] << " " << ins[i][3] << std::endl;
        return false;
      }
    }
    return true;
  }

/* Some notes:
 *   0. The model might be either in counts or probability form.
 *   1. t[0] is special: t[0][TMM,TMI,TMD] are the begin->M_1,I_0,D_1 entry probabilities,
 *      t[0][TIM,TII] are the I_0 transitions, and delete state 0 doesn't
 *      exist. Therefore D[0] transitions and mat[0] emissions are unused.
 *      To simplify some normalization code, we adopt a convention that these are set
 *      to valid probability distributions: 1.0 for t[0][TDM] and mat[0][0],
 *      and 0 for the rest.
 *   2. t[M] is also special: TMD and TDD are 0 because there is no next delete state;
 *      TDM is therefore 1.0 by definition. TMM and TDM are interpreted as the
 *      M->E and D->E end transitions. t[M][TDM] must be 1.0, therefore.
 */
  void reverse() {
    std::reverse(ins.begin(), ins.end());
    std::reverse(mat.begin() + 1, mat.end());
    for (size_t i = 0, j = M; i < j; ++i, --j) {
      std::swap(t[i][p7H_II], t[j][p7H_II]);
      std::swap(t[i][p7H_MM], t[j][p7H_MM]);
      std::swap(t[i][p7H_DD], t[j][p7H_DD]);
    }
    for (size_t i = 0; i <= M; ++i) {
      std::swap(t[i][p7H_MD], t[M - i][p7H_DM]);
      std::swap(t[i][p7H_IM], t[M - i][p7H_MI]);
    }

    // TODO Does it reaaly needed???
    t[0][p7H_DM] = 0;
    t[0][p7H_DD] = std::numeric_limits<double>::infinity();
    t[M][p7H_MD] = std::numeric_limits<double>::infinity();
    t[M][p7H_DD] = std::numeric_limits<double>::infinity();
  }

  Fees reversed() const {
    Fees copy{*this};
    copy.reverse();
    return copy;
  }
};

Fees levenshtein_fees(const std::string &s, double mismatch = 1, double gap_open = 1, double gap_ext = 1);
Fees fees_from_hmm(const P7_HMM *hmm, const ESL_ALPHABET *abc);
Fees read_hmm_file(const std::string &filename);
std::vector<std::string> read_fasta_edges(const std::string &filename, bool add_rc = false);
std::vector<std::pair<std::string, double>> find_best_path(const Fees &fees, const std::vector<DBGraph::GraphPointer> &initial);
std::vector<std::pair<std::string, double>> find_best_path_rev(const Fees &fees,
                                                               const std::vector<ReversalGraphPointer<DBGraph::GraphPointer>> &initial);
std::vector<std::pair<std::string, double>> find_best_path(const Fees &fees, const std::vector<Graph::GraphPointer> &initial);
std::vector<std::pair<std::string, double>> find_best_path_rev(const Fees &fees,
                                                               const std::vector<ReversalGraphPointer<Graph::GraphPointer>> &initial);
