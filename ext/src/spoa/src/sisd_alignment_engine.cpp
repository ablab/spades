// Copyright (c) 2020 Robert Vaser

#include "sisd_alignment_engine.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

#include "spoa/graph.hpp"

namespace spoa {

constexpr std::int32_t kNegativeInfinity =
    std::numeric_limits<std::int32_t>::min() + 1024;

std::unique_ptr<AlignmentEngine> SisdAlignmentEngine::Create(
    AlignmentType type,
    AlignmentSubtype subtype,
    std::int8_t m,
    std::int8_t n,
    std::int8_t g,
    std::int8_t e,
    std::int8_t q,
    std::int8_t c) {
  return std::unique_ptr<AlignmentEngine>(
      new SisdAlignmentEngine(type, subtype, m, n, g, e, q, c));
}

struct SisdAlignmentEngine::Implementation {
  std::vector<std::uint32_t> node_id_to_rank;
  std::vector<std::int32_t> sequence_profile;
  std::vector<std::int32_t> M;
  std::int32_t* H;
  std::int32_t* F;
  std::int32_t* E;
  std::int32_t* O;
  std::int32_t* Q;

  Implementation()
      : node_id_to_rank(),
        sequence_profile(),
        M(),
        H(nullptr),
        F(nullptr),
        E(nullptr),
        O(nullptr),
        Q(nullptr) {
  }
};

SisdAlignmentEngine::SisdAlignmentEngine(
    AlignmentType type,
    AlignmentSubtype subtype,
    std::int8_t m,
    std::int8_t n,
    std::int8_t g,
    std::int8_t e,
    std::int8_t q,
    std::int8_t c)
    : AlignmentEngine(type, subtype, m, n, g, e, q, c),
      pimpl_(new Implementation()) {
}

void SisdAlignmentEngine::Prealloc(
    std::uint32_t max_sequence_len,
    std::uint8_t alphabet_size) {
  if (max_sequence_len > std::numeric_limits<int32_t>::max()) {
    throw std::invalid_argument(
        "[spoa::SisdAlignmentEngine::Prealloc] error: too large sequence!");
  }
  try {
    Realloc(
        static_cast<std::uint64_t>(max_sequence_len) + 1,
        static_cast<std::uint64_t>(max_sequence_len) * alphabet_size + alphabet_size,  // NOLINT
        alphabet_size);
  } catch (std::bad_alloc& ba) {
    throw std::invalid_argument(
        "[spoa::SisdAlignmentEngine::Prealloc] error: insufficient memory!");
  }
}

void SisdAlignmentEngine::Realloc(
    std::uint64_t matrix_width,
    std::uint64_t matrix_height,
    std::uint8_t num_codes) {
  if (pimpl_->node_id_to_rank.size() < matrix_height - 1) {
    pimpl_->node_id_to_rank.resize(matrix_height - 1, 0);
  }
  if (pimpl_->sequence_profile.size() < num_codes * matrix_width) {
    pimpl_->sequence_profile.resize(num_codes * matrix_width, 0);
  }
  if (subtype_ == AlignmentSubtype::kLinear) {
    if (pimpl_->M.size() < matrix_height * matrix_width) {
      pimpl_->M.resize(matrix_width * matrix_height, 0);
      pimpl_->H = pimpl_->M.data();
      pimpl_->F = nullptr;
      pimpl_->E = nullptr;
    }
  } else if (subtype_ == AlignmentSubtype::kAffine) {
    if (pimpl_->M.size() < 3 * matrix_height * matrix_width) {
      pimpl_->M.resize(3 * matrix_width * matrix_height, 0);
      pimpl_->H = pimpl_->M.data();
      pimpl_->F = pimpl_->H + matrix_width * matrix_height;
      pimpl_->E = pimpl_->F + matrix_width * matrix_height;
    }
  } else if (subtype_ == AlignmentSubtype::kConvex) {
    if (pimpl_->M.size() < 5 * matrix_height * matrix_width) {
      pimpl_->M.resize(5 * matrix_width * matrix_height, 0);
      pimpl_->H = pimpl_->M.data();
      pimpl_->F = pimpl_->H + matrix_width * matrix_height;
      pimpl_->E = pimpl_->F + matrix_width * matrix_height;
      pimpl_->O = pimpl_->E + matrix_width * matrix_height;
      pimpl_->Q = pimpl_->O + matrix_width * matrix_height;
    }
  }
}

void SisdAlignmentEngine::Initialize(
    const char* sequence, std::uint32_t sequence_len,
    const Graph& graph) noexcept {
  std::uint32_t matrix_width = sequence_len + 1;
  std::uint32_t matrix_height = graph.nodes().size() + 1;

  for (std::uint32_t i = 0; i < graph.num_codes(); ++i) {
    char c = graph.decoder(i);
    pimpl_->sequence_profile[i * matrix_width] = 0;
    for (std::uint32_t j = 0; j < sequence_len; ++j) {
      pimpl_->sequence_profile[i * matrix_width + (j + 1)] =
          (c == sequence[j] ? m_ : n_);
    }
  }

  const auto& rank_to_node = graph.rank_to_node();
  for (std::uint32_t i = 0; i < rank_to_node.size(); ++i) {
    pimpl_->node_id_to_rank[rank_to_node[i]->id] = i;
  }

  // initialize secondary matrices
  switch (subtype_) {
    case AlignmentSubtype::kConvex:
      pimpl_->O[0] = 0;
      pimpl_->Q[0] = 0;
      for (std::uint32_t j = 1; j < matrix_width; ++j) {
        pimpl_->O[j] = kNegativeInfinity;
        pimpl_->Q[j] = q_ + (j - 1) * c_;
      }
      for (std::uint32_t i = 1; i < matrix_height; ++i) {
        const auto& edges = rank_to_node[i - 1]->inedges;
        std::int32_t penalty = edges.empty() ? q_ - c_ : kNegativeInfinity;
        for (const auto& it : edges) {
          std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;
          penalty = std::max(penalty, pimpl_->O[pred_i * matrix_width]);
        }
        pimpl_->O[i * matrix_width] = penalty + c_;
        pimpl_->Q[i * matrix_width] = kNegativeInfinity;
      }
      // fall through
    case AlignmentSubtype::kAffine:
      pimpl_->F[0] = 0;
      pimpl_->E[0] = 0;
      for (std::uint32_t j = 1; j < matrix_width; ++j) {
        pimpl_->F[j] = kNegativeInfinity;
        pimpl_->E[j] = g_ + (j - 1) * e_;
      }
      for (std::uint32_t i = 1; i < matrix_height; ++i) {
        const auto& edges = rank_to_node[i - 1]->inedges;
        std::int32_t penalty = edges.empty() ? g_ - e_ : kNegativeInfinity;
        for (const auto& it : edges) {
          std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;
          penalty = std::max(penalty, pimpl_->F[pred_i * matrix_width]);
        }
        pimpl_->F[i * matrix_width] = penalty + e_;
        pimpl_->E[i * matrix_width] = kNegativeInfinity;
      }
      // fall through
    case AlignmentSubtype::kLinear:
      pimpl_->H[0] = 0;
      break;
    default:
      break;
  }

  // initialize primary matrix
  switch (type_) {
    case AlignmentType::kSW:
      for (std::uint32_t j = 1; j < matrix_width; ++j) {
        pimpl_->H[j] = 0;
      }
      for (std::uint32_t i = 1; i < matrix_height; ++i) {
        pimpl_->H[i * matrix_width] = 0;
      }
      break;
    case AlignmentType::kNW:
      switch (subtype_) {
        case AlignmentSubtype::kConvex:
          for (std::uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->H[j] = std::max(pimpl_->Q[j], pimpl_->E[j]);
          }
          for (std::uint32_t i = 1; i < matrix_height; ++i) {
            pimpl_->H[i * matrix_width] = std::max(
                pimpl_->O[i * matrix_width],
                pimpl_->F[i * matrix_width]);
          }
          break;
        case AlignmentSubtype::kAffine:
          for (std::uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->H[j] = pimpl_->E[j];
          }
          for (std::uint32_t i = 1; i < matrix_height; ++i) {
            pimpl_->H[i * matrix_width] = pimpl_->F[i * matrix_width];
          }
          break;
        case AlignmentSubtype::kLinear:
          for (std::uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->H[j] = j * g_;
          }
          for (std::uint32_t i = 1; i < matrix_height; ++i) {
            const auto& edges = rank_to_node[i - 1]->inedges;
            std::int32_t penalty = edges.empty() ? 0 : kNegativeInfinity;
            for (const auto& it : edges) {
              std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;
              penalty = std::max(penalty, pimpl_->H[pred_i * matrix_width]);
            }
            pimpl_->H[i * matrix_width] = penalty + g_;
          }
        default:
          break;
      }
      break;
    case AlignmentType::kOV:
      switch (subtype_) {
        case AlignmentSubtype::kConvex:
          for (std::uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->H[j] = std::max(pimpl_->Q[j], pimpl_->E[j]);
          }
          break;
        case AlignmentSubtype::kAffine:
          for (std::uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->H[j] = pimpl_->E[j];
          }
          break;
        case AlignmentSubtype::kLinear:
          for (std::uint32_t j = 1; j < matrix_width; ++j) {
            pimpl_->H[j] = j * g_;
          }
          break;
        default:
          break;
      }
      for (std::uint32_t i = 1; i < matrix_height; ++i) {
        pimpl_->H[i * matrix_width] = 0;
      }
      break;
    default:
      break;
  }
}

Alignment SisdAlignmentEngine::Align(
    const char* sequence, std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score) {
  if (sequence_len > std::numeric_limits<int32_t>::max()) {
    throw std::invalid_argument(
        "[spoa::SisdAlignmentEngine::Align] error: too large sequence!");
  }

  if (graph.nodes().empty() || sequence_len == 0) {
    return Alignment();
  }

  if (WorstCaseAlignmentScore(sequence_len, graph.nodes().size()) < kNegativeInfinity) {  // NOLINT
    throw std::invalid_argument(
        "[spoa::SisdAlignmentEngine::Align] error: possible overflow!");
  }

  try {
    Realloc(sequence_len + 1, graph.nodes().size() + 1, graph.num_codes());
  } catch (std::bad_alloc& ba) {
    throw std::invalid_argument(
        "[spoa::SisdAlignmentEngine::Align] error: insufficient memory!");
  }
  Initialize(sequence, sequence_len, graph);

  if (subtype_ == AlignmentSubtype::kLinear) {
    return Linear(sequence_len, graph, score);
  } else if (subtype_ == AlignmentSubtype::kAffine) {
    return Affine(sequence_len, graph, score);
  } else if (subtype_ == AlignmentSubtype::kConvex) {
    return Convex(sequence_len, graph, score);
  }
  return Alignment();
}

Alignment SisdAlignmentEngine::Linear(
    std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score) noexcept {
  std::uint64_t matrix_width = sequence_len + 1;
  const auto& rank_to_node = graph.rank_to_node();

  std::int32_t max_score = type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;
  std::uint32_t max_i = 0;
  std::uint32_t max_j = 0;
  auto update_max_score = [&max_score, &max_i, &max_j] (
      std::int32_t* H_row,
      std::uint32_t i,
      std::uint32_t j) -> void {
    if (max_score < H_row[j]) {
      max_score = H_row[j];
      max_i = i;
      max_j = j;
    }
    return;
  };

  // alignment
  for (const auto& it : rank_to_node) {
    const auto& char_profile =
        &(pimpl_->sequence_profile[it->code * matrix_width]);

    std::uint32_t i = pimpl_->node_id_to_rank[it->id] + 1;
    std::uint32_t pred_i = it->inedges.empty() ?
        0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

    std::int32_t* H_row = &(pimpl_->H[i * matrix_width]);
    std::int32_t* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

    // update H
    for (std::uint64_t j = 1; j < matrix_width; ++j) {
      H_row[j] = std::max(
          H_pred_row[j - 1] + char_profile[j],
          H_pred_row[j] + g_);
    }
    // check other predeccessors
    for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
      pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

      H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

      for (std::uint64_t j = 1; j < matrix_width; ++j) {
        H_row[j] = std::max(
            H_pred_row[j - 1] + char_profile[j],
            std::max(
                H_row[j],
                H_pred_row[j] + g_));
      }
    }

    for (std::uint64_t j = 1; j < matrix_width; ++j) {
      H_row[j] = std::max(H_row[j - 1] + g_, H_row[j]);

      if (type_ == AlignmentType::kSW) {
        H_row[j] = std::max(H_row[j], 0);
        update_max_score(H_row, i, j);
      } else if (type_ == AlignmentType::kNW &&
          it->outedges.empty() && j == matrix_width - 1) {
        update_max_score(H_row, i, j);
      } else if (type_ == AlignmentType::kOV && it->outedges.empty()) {
        update_max_score(H_row, i, j);
      }
    }
  }

  if (max_i == 0 && max_j == 0) {
    return Alignment();
  }
  if (score) {
    *score = max_score;
  }

  // backtrack
  Alignment alignment;
  std::uint32_t i = max_i;
  std::uint32_t j = max_j;

  auto sw_condition = [this, &i, &j, &matrix_width] () -> bool {
    return (pimpl_->H[i * matrix_width + j] == 0) ? false : true;
  };
  auto nw_condition = [&i, &j] () -> bool {
    return (i == 0 && j == 0) ? false : true;
  };
  auto ov_condition = [&i, &j] () -> bool {
    return (i == 0 || j == 0) ? false : true;
  };

  std::uint32_t prev_i = 0;
  std::uint32_t prev_j = 0;
  while ((type_ == AlignmentType::kSW && sw_condition()) ||
         (type_ == AlignmentType::kNW && nw_condition()) ||
         (type_ == AlignmentType::kOV && ov_condition())) {
    auto H_ij = pimpl_->H[i * matrix_width + j];
    bool predecessor_found = false;

    if (i != 0 && j != 0) {
      const auto& it = rank_to_node[i - 1];
      std::int32_t match_cost =
          pimpl_->sequence_profile[it->code * matrix_width + j];

      std::uint32_t pred_i = it->inedges.empty() ?
          0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
        prev_i = pred_i;
        prev_j = j - 1;
        predecessor_found = true;
      } else {
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          std::uint32_t pred_i =
              pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
            prev_i = pred_i;
            prev_j = j - 1;
            predecessor_found = true;
            break;
          }
        }
      }
    }

    if (!predecessor_found && i != 0) {
      const auto& it = rank_to_node[i - 1];

      std::uint32_t pred_i = it->inedges.empty() ? 0 :
          pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      if (H_ij == pimpl_->H[pred_i * matrix_width + j] + g_) {
        prev_i = pred_i;
        prev_j = j;
        predecessor_found = true;
      } else {
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          std::uint32_t pred_i =
              pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          if (H_ij == pimpl_->H[pred_i * matrix_width + j] + g_) {
            prev_i = pred_i;
            prev_j = j;
            predecessor_found = true;
            break;
          }
        }
      }
    }

    if (!predecessor_found && H_ij == pimpl_->H[i * matrix_width + j - 1] + g_) {  // NOLINT
      prev_i = i;
      prev_j = j - 1;
      predecessor_found = true;
    }

    alignment.emplace_back(
        i == prev_i ? -1 : rank_to_node[i - 1]->id,
        j == prev_j ? -1 : j - 1);

    i = prev_i;
    j = prev_j;
  }

  std::reverse(alignment.begin(), alignment.end());
  return alignment;
}

Alignment SisdAlignmentEngine::Affine(
    std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score) noexcept {
  std::uint64_t matrix_width = sequence_len + 1;
  const auto& rank_to_node = graph.rank_to_node();

  std::int32_t max_score = type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;
  std::uint32_t max_i = 0;
  std::uint32_t max_j = 0;
  auto update_max_score = [&max_score, &max_i, &max_j] (
      std::int32_t* H_row,
      std::uint32_t i,
      std::uint32_t j) -> void {
    if (max_score < H_row[j]) {
      max_score = H_row[j];
      max_i = i;
      max_j = j;
    }
    return;
  };

  // alignment
  for (const auto& it : rank_to_node) {
    const auto& char_profile =
        &(pimpl_->sequence_profile[it->code * matrix_width]);

    std::uint32_t i = pimpl_->node_id_to_rank[it->id] + 1;
    std::uint32_t pred_i = it->inedges.empty() ? 0 :
        pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

    std::int32_t* H_row = &(pimpl_->H[i * matrix_width]);
    std::int32_t* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

    std::int32_t* F_row = &(pimpl_->F[i * matrix_width]);
    std::int32_t* F_pred_row = &(pimpl_->F[pred_i * matrix_width]);

    // update F and H
    for (std::uint64_t j = 1; j < matrix_width; ++j) {
      F_row[j] = std::max(
          H_pred_row[j] + g_,
          F_pred_row[j] + e_);
      H_row[j] = H_pred_row[j - 1] + char_profile[j];
    }
    // check other predeccessors
    for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
      pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

      H_pred_row = &(pimpl_->H[pred_i * matrix_width]);
      F_pred_row = &(pimpl_->F[pred_i * matrix_width]);

      for (std::uint64_t j = 1; j < matrix_width; ++j) {
        F_row[j] = std::max(
            F_row[j],
            std::max(
                H_pred_row[j] + g_,
                F_pred_row[j] + e_));
        H_row[j] = std::max(
            H_row[j],
            H_pred_row[j - 1] + char_profile[j]);
      }
    }

    // update E and H
    std::int32_t* E_row = &(pimpl_->E[i * matrix_width]);
    for (std::uint64_t j = 1; j < matrix_width; ++j) {
      E_row[j] = std::max(H_row[j - 1] + g_, E_row[j - 1] + e_);
      H_row[j] = std::max(H_row[j], std::max(F_row[j], E_row[j]));

      if (type_ == AlignmentType::kSW) {
        H_row[j] = std::max(H_row[j], 0);
        update_max_score(H_row, i, j);
      } else if (type_ == AlignmentType::kNW &&
          (it->outedges.empty() && j == matrix_width - 1)) {
        update_max_score(H_row, i, j);
      } else if (type_ == AlignmentType::kOV && (it->outedges.empty())) {
        update_max_score(H_row, i, j);
      }
    }
  }

  if (max_i == 0 && max_j == 0) {
    return Alignment();
  }
  if (score) {
    *score = max_score;
  }

  // backtrack
  Alignment alignment;
  std::uint32_t i = max_i;
  std::uint32_t j = max_j;

  auto sw_condition = [this, &i, &j, &matrix_width] () -> bool {
    return (pimpl_->H[i * matrix_width + j] == 0) ? false : true;
  };
  auto nw_condition = [&i, &j] () -> bool {
    return (i == 0 && j == 0) ? false : true;
  };
  auto ov_condition = [&i, &j] () -> bool {
    return (i == 0 || j == 0) ? false : true;
  };

  std::uint32_t prev_i = 0;
  std::uint32_t prev_j = 0;

  while ((type_ == AlignmentType::kSW && sw_condition()) ||
         (type_ == AlignmentType::kNW && nw_condition()) ||
         (type_ == AlignmentType::kOV && ov_condition())) {
    auto H_ij = pimpl_->H[i * matrix_width + j];
    bool predecessor_found = false, extend_left = false, extend_up = false;

    if (i != 0 && j != 0) {
      const auto& it = rank_to_node[i - 1];
      std::int32_t match_cost =
          pimpl_->sequence_profile[it->code * matrix_width + j];

      std::uint32_t pred_i = it->inedges.empty() ?
          0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
        prev_i = pred_i;
        prev_j = j - 1;
        predecessor_found = true;
      } else {
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
            prev_i = pred_i;
            prev_j = j - 1;
            predecessor_found = true;
            break;
          }
        }
      }
    }

    if (!predecessor_found && i != 0) {
      const auto& it = rank_to_node[i - 1];

      std::uint32_t pred_i = it->inedges.empty() ?
          0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      if ((extend_up = H_ij == pimpl_->F[pred_i * matrix_width + j] + e_) ||
                       H_ij == pimpl_->H[pred_i * matrix_width + j] + g_) {
        prev_i = pred_i;
        prev_j = j;
        predecessor_found = true;
      } else {
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          if ((extend_up = H_ij == pimpl_->F[pred_i * matrix_width + j] + e_) ||
                           H_ij == pimpl_->H[pred_i * matrix_width + j] + g_) {
            prev_i = pred_i;
            prev_j = j;
            predecessor_found = true;
            break;
          }
        }
      }
    }

    if (!predecessor_found && j != 0) {
      if ((extend_left = H_ij == pimpl_->E[i * matrix_width + j - 1] + e_) ||
                         H_ij == pimpl_->H[i * matrix_width + j - 1] + g_) {
        prev_i = i;
        prev_j = j - 1;
        predecessor_found = true;
      }
    }

    alignment.emplace_back(
        i == prev_i ? -1 : rank_to_node[i - 1]->id,
        j == prev_j ? -1 : j - 1);

    i = prev_i;
    j = prev_j;

    if (extend_left) {
      while (true) {
        alignment.emplace_back(-1, j - 1);
        --j;
        if (pimpl_->E[i * matrix_width + j] + e_ !=
            pimpl_->E[i * matrix_width + j + 1]) {
          break;
        }
      }
    } else if (extend_up) {
      while (true) {
        bool stop = false;
        prev_i = 0;
        for (const auto& it : rank_to_node[i - 1]->inedges) {
          std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;

          if ((stop = pimpl_->F[i * matrix_width + j] == pimpl_->H[pred_i * matrix_width + j] + g_) ||  // NOLINT
                      pimpl_->F[i * matrix_width + j] == pimpl_->F[pred_i * matrix_width + j] + e_) {  // NOLINT
            prev_i = pred_i;
            break;
          }
        }

        alignment.emplace_back(rank_to_node[i - 1]->id, -1);
        i = prev_i;
        if (stop || i == 0) {
          break;
        }
      }
    }
  }

  std::reverse(alignment.begin(), alignment.end());
  return alignment;
}

Alignment SisdAlignmentEngine::Convex(
    std::uint32_t sequence_len,
    const Graph& graph,
    std::int32_t* score) noexcept {
  std::uint64_t matrix_width = sequence_len + 1;
  const auto& rank_to_node = graph.rank_to_node();

  std::int32_t max_score = type_ == AlignmentType::kSW ? 0 : kNegativeInfinity;
  std::uint32_t max_i = 0;
  std::uint32_t max_j = 0;
  auto update_max_score = [&max_score, &max_i, &max_j] (
      std::int32_t* H_row,
      std::uint32_t i,
      std::uint32_t j) -> void {
    if (max_score < H_row[j]) {
      max_score = H_row[j];
      max_i = i;
      max_j = j;
    }
    return;
  };

  // alignment
  for (const auto& it : rank_to_node) {
    const auto& char_profile =
        &(pimpl_->sequence_profile[it->code * matrix_width]);

    std::uint32_t i = pimpl_->node_id_to_rank[it->id] + 1;
    std::uint32_t pred_i = it->inedges.empty() ? 0 :
        pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

    std::int32_t* H_row = &(pimpl_->H[i * matrix_width]);
    std::int32_t* H_pred_row = &(pimpl_->H[pred_i * matrix_width]);

    std::int32_t* F_row = &(pimpl_->F[i * matrix_width]);
    std::int32_t* F_pred_row = &(pimpl_->F[pred_i * matrix_width]);

    std::int32_t* O_row = &(pimpl_->O[i * matrix_width]);
    std::int32_t* O_pred_row = &(pimpl_->O[pred_i * matrix_width]);

    // update F, O and H
    for (std::uint64_t j = 1; j < matrix_width; ++j) {
      F_row[j] = std::max(H_pred_row[j] + g_, F_pred_row[j] + e_);
      O_row[j] = std::max(H_pred_row[j] + q_, O_pred_row[j] + c_);
      H_row[j] = H_pred_row[j - 1] + char_profile[j];
    }
    // check other predeccessors
    for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
      pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

      H_pred_row = &(pimpl_->H[pred_i * matrix_width]);
      F_pred_row = &(pimpl_->F[pred_i * matrix_width]);
      O_pred_row = &(pimpl_->O[pred_i * matrix_width]);

      for (std::uint64_t j = 1; j < matrix_width; ++j) {
        F_row[j] = std::max(
            F_row[j],
            std::max(
                H_pred_row[j] + g_,
                F_pred_row[j] + e_));
        O_row[j] = std::max(
            O_row[j],
            std::max(
                H_pred_row[j] + q_,
                O_pred_row[j] + c_));
        H_row[j] = std::max(H_row[j], H_pred_row[j - 1] + char_profile[j]);
      }
    }

    // update E, Q and H
    std::int32_t* E_row = &(pimpl_->E[i * matrix_width]);
    std::int32_t* Q_row = &(pimpl_->Q[i * matrix_width]);
    for (std::uint64_t j = 1; j < matrix_width; ++j) {
      E_row[j] = std::max(H_row[j - 1] + g_, E_row[j - 1] + e_);
      Q_row[j] = std::max(H_row[j - 1] + q_, Q_row[j - 1] + c_);
      H_row[j] = std::max(
        H_row[j],
        std::max(
            std::max(F_row[j], E_row[j]),
            std::max(O_row[j], Q_row[j])));

      if (type_ == AlignmentType::kSW) {
        H_row[j] = std::max(H_row[j], 0);
        update_max_score(H_row, i, j);
      } else if (type_ == AlignmentType::kNW &&
          (it->outedges.empty() && j == matrix_width - 1)) {
        update_max_score(H_row, i, j);
      } else if (type_ == AlignmentType::kOV && it->outedges.empty()) {
        update_max_score(H_row, i, j);
      }
    }
  }

  if (max_i == 0 && max_j == 0) {
    return Alignment();
  }
  if (score) {
    *score = max_score;
  }

  // backtrack
  Alignment alignment;
  std::uint32_t i = max_i;
  std::uint32_t j = max_j;

  auto sw_condition = [this, &i, &j, &matrix_width] () -> bool {
    return (pimpl_->H[i * matrix_width + j] == 0) ? false : true;
  };
  auto nw_condition = [&i, &j] () -> bool {
    return (i == 0 && j == 0) ? false : true;
  };
  auto ov_condition = [&i, &j] () -> bool {
    return (i == 0 || j == 0) ? false : true;
  };

  std::uint32_t prev_i = 0;
  std::uint32_t prev_j = 0;

  while ((type_ == AlignmentType::kSW && sw_condition()) ||
         (type_ == AlignmentType::kNW && nw_condition()) ||
         (type_ == AlignmentType::kOV && ov_condition())) {
    auto H_ij = pimpl_->H[i * matrix_width + j];
    bool predecessor_found = false, extend_left = false, extend_up = false;

    if (i != 0 && j != 0) {
      const auto& it = rank_to_node[i - 1];
      std::int32_t match_cost =
        pimpl_->sequence_profile[it->code * matrix_width + j];

      std::uint32_t pred_i = it->inedges.empty() ?
          0 : pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
        prev_i = pred_i;
        prev_j = j - 1;
        predecessor_found = true;
      } else {
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          if (H_ij == pimpl_->H[pred_i * matrix_width + (j - 1)] + match_cost) {
            prev_i = pred_i;
            prev_j = j - 1;
            predecessor_found = true;
            break;
          }
        }
      }
    }

    if (!predecessor_found && i != 0) {
      const auto& it = rank_to_node[i - 1];

      std::uint32_t pred_i = it->inedges.empty() ? 0 :
        pimpl_->node_id_to_rank[it->inedges[0]->tail->id] + 1;

      if ((extend_up |= H_ij == pimpl_->F[pred_i * matrix_width + j] + e_) ||
                        H_ij == pimpl_->H[pred_i * matrix_width + j] + g_  ||
          (extend_up |= H_ij == pimpl_->O[pred_i * matrix_width + j] + c_) ||
                        H_ij == pimpl_->H[pred_i * matrix_width + j] + q_) {
        prev_i = pred_i;
        prev_j = j;
        predecessor_found = true;
      } else {
        for (std::uint32_t p = 1; p < it->inedges.size(); ++p) {
          pred_i = pimpl_->node_id_to_rank[it->inedges[p]->tail->id] + 1;

          if ((extend_up |= H_ij == pimpl_->F[pred_i * matrix_width + j] + e_) ||  // NOLINT
                            H_ij == pimpl_->H[pred_i * matrix_width + j] + g_  ||  // NOLINT
              (extend_up |= H_ij == pimpl_->O[pred_i * matrix_width + j] + c_) ||  // NOLINT
                            H_ij == pimpl_->H[pred_i * matrix_width + j] + q_) {
            prev_i = pred_i;
            prev_j = j;
            predecessor_found = true;
            break;
          }
        }
      }
    }

    if (!predecessor_found && j != 0) {
      if ((extend_left |= H_ij == pimpl_->E[i * matrix_width + j - 1] + e_) ||
                          H_ij == pimpl_->H[i * matrix_width + j - 1] + g_  ||
          (extend_left |= H_ij == pimpl_->Q[i * matrix_width + j - 1] + c_) ||
                          H_ij == pimpl_->H[i * matrix_width + j - 1] + q_) {
        prev_i = i;
        prev_j = j - 1;
        predecessor_found = true;
      }
    }

    alignment.emplace_back(
        i == prev_i ? -1 : rank_to_node[i - 1]->id,
        j == prev_j ? -1 : j - 1);

    i = prev_i;
    j = prev_j;

    if (extend_left) {
      while (true) {
        alignment.emplace_back(-1, j - 1);
        --j;
        if (pimpl_->E[i * matrix_width + j] + e_ != pimpl_->E[i * matrix_width + j + 1] &&  // NOLINT
            pimpl_->Q[i * matrix_width + j] + c_ != pimpl_->Q[i * matrix_width + j + 1]) {  // NOLINT
          break;
        }
      }
    } else if (extend_up) {
      while (true) {
        bool stop = true;
        prev_i = 0;
        for (const auto& it : rank_to_node[i - 1]->inedges) {
          std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;

          if (pimpl_->F[i * matrix_width + j] == pimpl_->F[pred_i * matrix_width + j] + e_ ||  // NOLINT
              pimpl_->O[i * matrix_width + j] == pimpl_->O[pred_i * matrix_width + j] + c_) {  // NOLINT
            prev_i = pred_i;
            stop = false;
            break;
          }
        }
        if (stop == true) {
          for (const auto& it : rank_to_node[i - 1]->inedges) {
            std::uint32_t pred_i = pimpl_->node_id_to_rank[it->tail->id] + 1;

            if (pimpl_->F[i * matrix_width + j] == pimpl_->H[pred_i * matrix_width + j] + g_ ||  // NOLINT
                pimpl_->O[i * matrix_width + j] == pimpl_->H[pred_i * matrix_width + j] + q_) {  // NOLINT
              prev_i = pred_i;
              break;
            }
          }
        }

        alignment.emplace_back(rank_to_node[i - 1]->id, -1);
        i = prev_i;

        if (stop || i == 0) {
          break;
        }
      }
    }
  }

  std::reverse(alignment.begin(), alignment.end());
  return alignment;
}

}  // namespace spoa
