#include "fees.hpp"

#include "hmmfile.hpp"

#include <algorithm>

extern "C" {
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "hmmer.h"
}

namespace hmm {

DigitalCodind::DigitalCodind(const ESL_ALPHABET *abc) : inmap_(abc->inmap, abc->inmap + 128), k_{abc->K} {
  if (k_ == 20) {
    // Fix map to be consistent with aa.hpp
    inmap_['*'] = 20;
  }
}

DigitalCodind::DigitalCodind()
    : inmap_(128), k_{4} {
  std::string alphabet = "ACGT";
  for (size_t i = 0; i < alphabet.length(); ++i) {
    inmap_[alphabet[i]] = i;
  }
}

bool Fees::check_i_loop(size_t i) const {
  return t[i][p7H_II] + *min_element(ins[i].cbegin(), ins[i].cend()) > 0;
}

bool Fees::check_i_negative_loops() const {
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
void Fees::reverse() {
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

Fees levenshtein_fees(const std::string &s, double mismatch, double gap_open, double gap_ext) {
  size_t M = s.size();
  Fees fees;
  fees.M = M;
  fees.k = 4;
  DigitalCodind encode;
  fees.code = encode;

  fees.t.resize(M + 1);
  fees.mat.resize(M + 1);
  fees.ins.resize(M + 1);
  for (size_t i = 0; i <= M; ++i) {
    fees.mat[i].resize(4, mismatch);
    fees.ins[i].resize(4, 0);
    fees.t[i].resize(p7H_NTRANSITIONS);
  }

  const double inf = std::numeric_limits<double>::infinity();
  fees.mat[0] = {inf, inf, inf, inf};
  for (size_t i = 1; i <= M; ++i) {
    fees.mat[i][encode(s[i - 1])] = 0;
  }

  // Insertions
  for (size_t i = 0; i <= M; ++i) {
    fees.t[i][p7H_IM] = 0;
    fees.t[i][p7H_II] = gap_ext;
    fees.t[i][p7H_MI] = gap_open;
  }

  // Deletions
  fees.t[0][p7H_MD] = gap_open;
  fees.t[0][p7H_DD] = inf;
  fees.t[0][p7H_DM] = inf;
  for (size_t i = 1; i < M; ++i) {
    fees.t[i][p7H_MD] = gap_open;
    fees.t[i][p7H_DD] = gap_ext;
    fees.t[i][p7H_DM] = 0;
  }
  fees.t[M][p7H_MD] = inf;
  fees.t[M][p7H_DD] = inf;
  fees.t[M][p7H_DM] = 0;

  // Matches
  for (size_t i = 0; i <= M; ++i) {
    fees.t[i][p7H_MM] = 0;
  }

  return fees;
}

Fees fees_from_hmm(const P7_HMM *hmm, const ESL_ALPHABET *abc, double lambda) {
  size_t M = hmm->M;
  Fees fees;
  fees.code = DigitalCodind(abc);
  fees.M = M;

  size_t k = abc->K;
  fees.k = k;
  size_t all_k = (k == 4) ? 4 : 21;

  fees.t.resize(M + 1);
  fees.mat.resize(M + 1);
  fees.ins.resize(M + 1);
  for (size_t i = 0; i <= M; ++i) {
    fees.mat[i].resize(all_k);
    fees.ins[i].resize(all_k);
    fees.t[i].resize(p7H_NTRANSITIONS);
  }

  for (size_t i = 0; i <= M; ++i) {
    for (size_t j = 0; j < p7H_NTRANSITIONS; ++j) {
      fees.t[i][j] = -log(hmm->t[i][j]);
    }

    for (size_t j = 0; j < all_k; ++j) {
      fees.mat[i][j] = -log(hmm->mat[i][j]) + log(1. / k) + lambda;
      fees.ins[i][j] = -log(hmm->ins[i][j]) + log(1. / k) + lambda;
    }

    if (all_k == 21) {
      fees.mat[i][20] = std::numeric_limits<double>::infinity();
      fees.ins[i][20] = std::numeric_limits<double>::infinity();
    }
  }

  return fees;
}

Fees fees_from_file(const std::string &filename) {
  hmmer::HMMFile hmm_file(filename);
  auto hmm = hmm_file.read();
  return fees_from_hmm(hmm->get(), hmm->abc());
}

};
