// Copyright (c) 2020 Robert Vaser

#ifndef SPOA_ALIGNMENT_ENGINE_HPP_
#define SPOA_ALIGNMENT_ENGINE_HPP_

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace spoa {

enum class AlignmentType {
  kSW,  // Smith Waterman
  kNW,  // Needleman Wunsch
  kOV   // Overlap
};

enum class AlignmentSubtype {
  kLinear,  // g * i
  kAffine,  // g + (i - 1) * e
  kConvex   // min(g1 + (i - 1) * e1, g2 + (i - 1) * e2)
};

class Graph;
using Alignment = std::vector<std::pair<std::int32_t, std::int32_t>>;

class AlignmentEngine {
 public:
  virtual ~AlignmentEngine() = default;

  static std::unique_ptr<AlignmentEngine> Create(
      AlignmentType type,
      std::int8_t m,   // match
      std::int8_t n,   // mismatch
      std::int8_t g);  // gap

  static std::unique_ptr<AlignmentEngine> Create(
      AlignmentType type,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,   // gap open
      std::int8_t e);  // gap extend

  static std::unique_ptr<AlignmentEngine> Create(
      AlignmentType type,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,
      std::int8_t e,
      std::int8_t q,   // gap open of second affine
      std::int8_t c);  // gap extend of second affine

  virtual void Prealloc(
      std::uint32_t max_sequence_len,
      std::uint8_t alphabet_size) = 0;

  Alignment Align(
      const std::string& sequence,
      const Graph& graph,
      std::int32_t* score = nullptr);

  virtual Alignment Align(
      const char* sequence, std::uint32_t sequence_len,
      const Graph& graph,
      std::int32_t* score = nullptr) = 0;

 protected:
  AlignmentEngine(
      AlignmentType type,
      AlignmentSubtype subtype,
      std::int8_t m,
      std::int8_t n,
      std::int8_t g,
      std::int8_t e,
      std::int8_t q,
      std::int8_t c);

  std::int64_t WorstCaseAlignmentScore(
      std::int64_t sequence_len,
      std::int64_t graph_len) const;

  AlignmentType type_;
  AlignmentSubtype subtype_;
  std::int8_t m_;
  std::int8_t n_;
  std::int8_t g_;
  std::int8_t e_;
  std::int8_t q_;
  std::int8_t c_;
};

}  // namespace spoa

#endif  // SPOA_ALIGNMENT_ENGINE_HPP_
