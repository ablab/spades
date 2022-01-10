// Copyright (c) 2020 Robert Vaser

#include "spoa/alignment_engine.hpp"

#include <algorithm>
#include <exception>
#include <limits>
#include <stdexcept>

#include "sisd_alignment_engine.hpp"
#include "simd_alignment_engine.hpp"

namespace spoa {

std::unique_ptr<AlignmentEngine> AlignmentEngine::Create(
    AlignmentType type,
    std::int8_t m,
    std::int8_t n,
    std::int8_t g) {
  return Create(type, m, n, g, g);
}

std::unique_ptr<AlignmentEngine> AlignmentEngine::Create(
    AlignmentType type,
    std::int8_t m,
    std::int8_t n,
    std::int8_t g,
    std::int8_t e) {
  return Create(type, m, n, g, e, g, e);
}

std::unique_ptr<AlignmentEngine> AlignmentEngine::Create(
    AlignmentType type,
    std::int8_t m,
    std::int8_t n,
    std::int8_t g,
    std::int8_t e,
    std::int8_t q,
    std::int8_t c) {
  if (type != AlignmentType::kSW &&
      type != AlignmentType::kNW &&
      type != AlignmentType::kOV) {
    throw std::invalid_argument(
        "[spoa::AlignmentEngine::Create] error: invalid alignment type!");
  }
  if (g > 0 || q > 0) {
    throw std::invalid_argument(
        "[spoa::AlignmentEngine::Create] error: "
        "gap opening penalty must be non-positive!");
  }
  if (e > 0 || c > 0) {
    throw std::invalid_argument(
        "[spoa::AlignmentEngine::Create] error: "
        "gap extension penalty must be non-positive!");
  }

  AlignmentSubtype subtype = g >= e ?
      AlignmentSubtype::kLinear : (g <= q || e >= c ?
      AlignmentSubtype::kAffine : AlignmentSubtype::kConvex);

  if (subtype == AlignmentSubtype::kLinear) {
    e = g;
  } else if (subtype == AlignmentSubtype::kAffine) {
    q = g;
    c = e;
  }

  auto dst = CreateSimdAlignmentEngine(type, subtype, m, n, g, e, q, c);
  if (!dst) {
    return SisdAlignmentEngine::Create(type, subtype, m, n, g, e, q, c);
  }
  return dst;
}

AlignmentEngine::AlignmentEngine(
    AlignmentType type,
    AlignmentSubtype subtype,
    std::int8_t m,
    std::int8_t n,
    std::int8_t g,
    std::int8_t e,
    std::int8_t q,
    std::int8_t c)
    : type_(type),
      subtype_(subtype),
      m_(m),
      n_(n),
      g_(g),
      e_(e),
      q_(q),
      c_(c) {
}

Alignment AlignmentEngine::Align(
    const std::string& sequence,
    const Graph& graph,
    std::int32_t* score) {
  return Align(sequence.c_str(), sequence.size(), graph, score);
}

std::int64_t AlignmentEngine::WorstCaseAlignmentScore(
    std::int64_t i,
    std::int64_t j) const {
  auto gap_score = [&] (std::int64_t len) -> std::int64_t {
    return len == 0 ? 0 : std::min(g_ + (len - 1) * e_, q_ + (len - 1) * c_);
  };
  return std::min(
      -1 * (m_ * std::min(i, j) + gap_score(std::abs(i - j))),
      gap_score(i) + gap_score(j));
}

}  // namespace spoa
