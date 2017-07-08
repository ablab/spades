//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_IT_READ_CORRECTOR_HPP__
#define __HAMMER_IT_READ_CORRECTOR_HPP__

#include "HSeq.hpp"
#include "config_struct.hpp"
#include "consensus.hpp"
#include "flow_space_read.hpp"
#include "hkmer_distance.hpp"
#include "io/reads/single_read.hpp"
#include "valid_hkmer_generator.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/optional.hpp>

#include <bamtools/api/BamAlignment.h>
#include <bamtools/api/SamHeader.h>
#include "seqeval/BaseHypothesisEvaluator.h"

#include <algorithm>
#include <cassert>
#include <deque>
#include <fstream>
#include <iterator>
#include <limits>
#include <list>
#include <string>
#include <vector>

#if 1
#include <iomanip>
#include <iostream>
#include "sequence/nucl.hpp"
#endif

namespace hammer {
namespace correction {

namespace numeric = boost::numeric::ublas;

typedef numeric::matrix<double> ScoreMatrix;
typedef std::vector<ScoreMatrix> ScoreStorage;

template <typename It1, typename It2>
static bool exactAlignH(It1 a_begin, It1 a_initial_pos, It1 a_end,
                        It2 b_initial_pos, It2 /*b_end*/, uint8_t max_offset,
                        uint8_t n_cmp, int *p_offset) {
  int M = max_offset * 2 + 1;
  for (int i = 0; i < M; i++) {
    int offset = (i / 2) * ((i & 1) ? 1 : -1);  // 0, -1, 1, -2, 2, ...
    auto a_it = a_initial_pos + offset;
    auto b_it = b_initial_pos;
    if (a_it < a_begin || a_it + n_cmp > a_end) continue;
    bool match = true;
    for (size_t j = 0; j < n_cmp; j++)
      if ((a_it + j)->raw != (b_it + j)->raw) {
        match = false;
        break;
      }
    if (match) {
      *p_offset = offset;
      return true;
    }
  }
  return false;
}

template <typename It1, typename It2>
static int overlapAlignH(It1 a_begin, It1 a_end, It2 b_begin, It2 b_end,
                         uint8_t max_offset) {
  // TODO: use dynamic programming
  int M = max_offset * 2 + 1;
  int best_offset = 0;
  int best_score = 0;
  for (int i = 0; i < M; i++) {
    int offset = (i / 2) * ((i & 1) ? 1 : -1);  // 0, -1, 1, -2, 2, ...
    auto a_it = offset < 0 ? a_begin : a_begin + offset;
    auto b_it = offset < 0 ? b_begin - offset : b_begin;
    if (b_it < b_begin || a_it >= a_end) continue;
    int score = 0;
    for (; a_it != a_end && b_it != b_end; ++a_it, ++b_it)
      if (a_it->nucl == b_it->nucl) score += std::min(a_it->len, b_it->len);
    score -= i / 4;
    if (score > best_score) {
      best_offset = offset;
      best_score = score;
    }
  }
  return best_offset;
}

struct Score {
  short value;
  short dir;
  Score(short v, short d) : value(v), dir(d) {}
};

#if 1
template <typename It1, typename It2>
static void dump(boost::numeric::ublas::matrix<Score> &scores, It1 x_begin,
                 It1 x_end, It2 y_begin, It2 y_end) {
  std::cerr << "        ";
  for (auto it = x_begin; it != x_end; ++it)
    std::cerr << std::setw(3) << int(it->len) << nucl(it->nucl);
  std::cerr << "\n    ";
  auto m = x_end - x_begin;
  auto n = y_end - y_begin;
  for (int i = 0; i <= m; i++) std::cerr << std::setw(4) << scores(i, 0).value;
  std::cerr << '\n';
  for (int i = 1; i <= n; i++) {
    auto run = *(y_begin + i - 1);
    std::cerr << std::setw(2) << int(run.len) << nucl(run.nucl) << ' ';
    for (int j = 0; j <= m; j++)
      std::cerr << std::setw(4) << scores(j, i).value;
    std::cerr << '\n';
  }
}
#endif

template <typename It1, typename It2>
static int alignH(It1 read_begin, It1 read_end, It2 consensus_begin,
                  It2 consensus_end, int approx_read_offset,
                  size_t n_skip_consensus, uint8_t n_side = 5,
                  uint8_t n_cmp = 8) {
  int left_offset = n_side;
  int read_len = int(read_end - read_begin);
  int consensus_len = int(consensus_end - consensus_begin);

  It1 x_begin = read_begin + std::max(approx_read_offset - n_side, 0);
  if (x_begin == read_begin) left_offset = approx_read_offset;

  if (approx_read_offset - n_side + n_cmp >= read_len) {
    x_begin = read_end - std::min(n_cmp + 2 * n_side, read_len);
    left_offset = int(read_begin + approx_read_offset - x_begin);
  }

  auto x_end =
      x_begin + std::min(int(2 * n_side + n_cmp), int(read_end - x_begin));

  auto y_begin =
      consensus_begin + std::min(int(n_skip_consensus), consensus_len);
  if (y_begin == consensus_end) return 0;  // weird situation
  auto y_end = y_begin + std::min(int(n_cmp), int(consensus_end - y_begin));

  // glocal alignment of homopolymer runs
  const short kDirUpLeft = 0;
  const short kDirUp = 1;
  const short kDirLeft = 2;

  const short kBaseDiff = -3;
  const short kRunInsertionStart = -4;
  const short kRunInsertionExtend = -5;
  const short kRunDeletionStart = -4;
  const short kRunDeletionExtend = -5;
  const short kNuclMismatch = -5;
  const short kNuclMatch = 1;
  const short kFullMatch = 5;

  int m = int(x_end - x_begin);
  int n = int(y_end - y_begin);

  using namespace boost::numeric::ublas;
  matrix<Score> scores(m + 1, n + 1, Score(0, 0));

  size_t highest_x = 0, highest_y = 0;
  int highest_entry = std::numeric_limits<int>::min();

  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= n; j++) {
      int best_score = std::numeric_limits<int>::min();
      short best_dir = 0;

      auto run_x = *(x_begin + i - 1);
      auto run_y = *(y_begin + j - 1);

      int score;
      if (run_x.raw == run_y.raw) {
        score = kNuclMatch * run_x.len + scores(i - 1, j - 1).value;
        score += kFullMatch;
        if (score > best_score) {
          best_score = score;
          best_dir = kDirUpLeft;
        }
      } else if (run_x.nucl == run_y.nucl) {
        score = kBaseDiff * std::abs(run_x.len - run_y.len);
        score += kNuclMatch * std::min(run_x.len, run_y.len);
        score += scores(i - 1, j - 1).value;
        if (score > best_score) {
          best_score = score;
          best_dir = kDirUpLeft;
        }
      } else {
        score = scores(i - 1, j - 1).value;
        score += kNuclMismatch * std::max(run_x.len, run_y.len);

        if (score > best_score) {
          best_score = score;
          best_dir = kDirUpLeft;
        }
      }

      int multiplier;

      if (scores(i - 1, j).dir == kDirUp)
        multiplier = kRunDeletionExtend;
      else
        multiplier = kRunDeletionStart;
      score = scores(i - 1, j).value + multiplier * run_x.len;
      if (score > best_score) {
        best_score = score;
        best_dir = kDirUp;
      }

      if (scores(i, j - 1).dir == kDirLeft)
        multiplier = kRunInsertionStart;
      else
        multiplier = kRunInsertionExtend;
      score = scores(i, j - 1).value + multiplier * run_y.len;
      if (score > best_score) {
        best_score = score;
        best_dir = kDirLeft;
      }

      scores(i, j) = Score(static_cast<short>(best_score), best_dir);

      if (i == m || j == n) {
        const int kOffset = 4;
        int approx_offset = i - j - left_offset;
        int offset_penalty = std::abs(approx_offset) * kOffset;
        if (best_score - offset_penalty > highest_entry) {
          highest_entry = best_score - offset_penalty;
          highest_x = i;
          highest_y = j;
        }
      }
    }
  }

  int min_acceptable_score = ((kNuclMatch + kFullMatch) * n_cmp * 4) / 5;
  if (scores(highest_x, highest_y).value < min_acceptable_score && n_cmp < 16U)
    return alignH(read_begin, read_end, consensus_begin, consensus_end,
                  approx_read_offset, n_skip_consensus, n_side,
                  uint8_t(n_cmp * 2));

  int x = int(highest_x);
  int y = int(highest_y);
  while (x > 0 && y > 0) {
    int dir = scores(x, y).dir;
    switch (dir) {
      case kDirUp:
        --x;
        break;
      case kDirLeft:
        --y;
        break;
      case kDirUpLeft:
        --x, --y;
        break;
      default:
        break;
    }
  }

#if 0
  if (std::abs(x - y - left_offset) >= 4)
    dump(scores, x_begin, x_end, y_begin, y_end);
#endif

  return x - y - left_offset;
}

// Not used now
class HKMerProlonger {
  const KMerData &kmer_data_;

 public:
  struct RightSide {
    static size_t changingPosition() { return hammer::K - 1; }
    static hammer::HKMer shift(const hammer::HKMer &kmer) {
      hammer::HKMer res;
      for (size_t i = 1; i < hammer::K; ++i) res[i - 1] = kmer[i];
      return res;
    }
    template <typename T, typename U>
    static void append(T &cont, U obj) {
      cont.push_back(obj);
    }
  };

  struct LeftSide {
    static size_t changingPosition() { return 0; }
    static hammer::HKMer shift(const hammer::HKMer &kmer) {
      hammer::HKMer res;
      for (size_t i = 1; i < hammer::K; ++i) res[i] = kmer[i - 1];
      return res;
    }
    template <typename T, typename U>
    static void append(T &cont, U obj) {
      cont.push_front(obj);
    }
  };

 public:
  /// @param[in] seed               kmer to prolong
  /// @param[in] bases_to_recover   maximum number of bases to recover
  template <typename Side>
  std::deque<hammer::HomopolymerRun> prolong(const hammer::HKMer &seed,
                                             size_t bases_to_recover) {
    std::deque<hammer::HomopolymerRun> good_runs(hammer::K);
    for (size_t i = 0; i < hammer::K; ++i) good_runs[i] = seed[i];

    auto good_kmer = seed;
    auto changing_pos = Side::changingPosition();

    for (size_t recov = 0; recov < bases_to_recover; ++recov) {
      double inf = -std::numeric_limits<double>::infinity();
      double best_qual = inf;
      int best_nucl = -1;
      int best_len = -1;
      double next_best_qual = inf;

      auto kmer = Side::shift(good_kmer);

      for (size_t nucl = 0; nucl < 4; ++nucl) {
        if (nucl == good_kmer[changing_pos].nucl) continue;
        for (size_t len = 1; len <= 4; ++len) {
          kmer[changing_pos] = hammer::HomopolymerRun(nucl, len);
          auto &k = kmer_data_[kmer];
          auto qual = k.count * (1 - k.qual);
          if (qual > best_qual) {
            next_best_qual = best_qual;
            best_qual = qual;
            best_nucl = nucl;
            best_len = len;
          }
        }
      }

      // stop if high-quality kmer is not unique
      if (best_nucl == -1 || best_qual - next_best_qual < 0.8 * best_qual)
        break;

      kmer[changing_pos] = hammer::HomopolymerRun(best_nucl, best_len);
      Side::append(good_runs, kmer[changing_pos]);
      good_kmer = kmer;
    }

    return good_runs;
  }

 public:
  HKMerProlonger(const KMerData &kmer_data) : kmer_data_(kmer_data) {}
};

static const double kLowScoreThreshold = 1.0;

class CorrectedRead {
  FlowSpaceRead raw_read_;  // Uncorrected read
  const KMerData &kmer_data_;
  bool debug_mode_;

  // Stores runs after joining chunks
  std::vector<hammer::HomopolymerRun> corrected_runs_;

  // Contiguous part of read with strong consensus
  struct ConsensusChunk {
    int approx_read_offset;  // in the vector of raw read runs
    int approx_end_read_offset_;
    unsigned rollback_end;  // remove if don't align well

    int initial_read_offset_;

    enum { kChunkLeftAligned, kChunkRightAligned, kChunkNotAligned } alignment;

    const FlowSpaceRead &raw_read;
    size_t trimmed_left;
    size_t trimmed_right;
    bool debug_mode;

    std::vector<hammer::HomopolymerRun> consensus;
    std::vector<double> consensus_scores;

    int raw_start_offset() const { return initial_read_offset_; }

    ConsensusChunk(int initial_read_offset, int approximate_read_offset,
                   int approximate_end_read_offset, const ScoreStorage &scores,
                   unsigned rollback_end, const FlowSpaceRead &read,
                   bool debug_mode)
        : approx_read_offset(approximate_read_offset),
          approx_end_read_offset_(approximate_end_read_offset),
          rollback_end(rollback_end),
          initial_read_offset_(initial_read_offset),
          alignment(kChunkNotAligned),
          raw_read(read),
          trimmed_left(0),
          trimmed_right(0),
          debug_mode(debug_mode) {
      bool left_trim = true;
      for (size_t i = 0; i < scores.size(); ++i) {
        auto run = hammer::iontorrent::consensus(scores[i]);

        // trim low-quality runs from the left side
        if (run.second <= kLowScoreThreshold && left_trim) {
          approx_read_offset += 1;
          trimmed_left += 1;
          continue;
        }

        if (debug_mode && left_trim) {
          std::cerr << "[ConsensusChunk] trimmed from left: " << trimmed_left
                    << std::endl;
          std::cerr << "[ConsensusChunk] approx. read offset: "
                    << approx_read_offset << std::endl;
        }

        left_trim = false;
        VERIFY(run.first.len > 0);
        consensus.push_back(run.first);
        consensus_scores.push_back(run.second);
      }

      size_t right_end = consensus_scores.size();
      if (right_end == 0) return;

      while (consensus_scores[right_end - 1] <= kLowScoreThreshold) {
        --right_end;
        if (right_end == 0) break;
      }

      trimmed_right = consensus.size() - right_end;
      consensus.resize(right_end);
      consensus_scores.resize(right_end);
    }

    void AlignLeftEndAgainstRead(size_t skip = 0) {
      const auto &data = raw_read.data();

      int offset = alignH(data.begin(), data.end(), consensus.begin(),
                          consensus.end(), approx_read_offset, skip);

      if (debug_mode) {
        std::cerr << "[approx. read offset (left)] before: "
                  << approx_read_offset
                  << "; after: " << approx_read_offset + offset << std::endl;
      }

      approx_read_offset += offset;
      alignment = kChunkLeftAligned;
    }

    void AlignRightEndAgainstRead(size_t skip = 0) {
      const auto &data = raw_read.data();
      int position_on_read = approx_end_read_offset_ - 1;
      int offset = alignH(data.rbegin(), data.rend(), consensus.rbegin(),
                          consensus.rend(),
                          int(data.size()) - 1 - position_on_read, skip);
      if (debug_mode) {
        std::cerr << "[approx. read offset (right)] before: "
                  << approx_read_offset
                  << "; after: " << approx_read_offset - offset << std::endl;
      }
      approx_read_offset -= offset;
      alignment = kChunkRightAligned;
    }

    int approx_end_read_offset() const { return approx_end_read_offset_; }

    int approx_end_read_offset_untrimmed() const {
      return approx_end_read_offset() + int(trimmed_right);
    }

   private:
    void RollBack() {
      trimmed_right += rollback_end;
      auto old_size = consensus.size();
      VERIFY(old_size >= rollback_end);
      consensus.resize(old_size - rollback_end);
      approx_end_read_offset_ -= rollback_end;
      consensus_scores.resize(old_size - rollback_end);
      rollback_end = 0;
    }

    bool DoMerge(ConsensusChunk &chunk) {
      int right_end_offset = approx_end_read_offset();

      if (debug_mode) {
        std::cerr << "============== Merging chunks ==============="
                  << std::endl;
        std::cerr << "(" << approx_read_offset << " .. " << right_end_offset
                  << ")";
        std::cerr << " -- (" << chunk.approx_read_offset << " .. "
                  << chunk.approx_end_read_offset() << ")" << std::endl;

        int white_l = 0;
        for (int i = right_end_offset - 1; i >= 0; --i)
          white_l += raw_read[i].len;
        for (size_t i = 0; i < consensus.size(); ++i)
          white_l -= consensus[i].len;
        for (int i = 0; i < white_l; ++i) std::cerr << ' ';
        for (size_t i = std::max(-white_l, 0); i < consensus.size(); ++i)
          std::cerr << consensus[i].str();
        std::cerr << std::endl;

        for (int i = 0; i < chunk.approx_read_offset; ++i)
          for (int j = 0; j < raw_read[i].len; ++j) std::cerr << ' ';
        for (size_t i = 0; i < chunk.consensus.size(); ++i)
          std::cerr << chunk.consensus[i].str();
        std::cerr << std::endl;
      }

      if (right_end_offset <= chunk.approx_read_offset) {
        for (int i = right_end_offset; i < chunk.approx_read_offset; ++i) {
          if (i >= static_cast<int>(raw_read.size())) return false;
          consensus.push_back(raw_read[i]);
          alignment = kChunkNotAligned;

          // TODO: maintain quality scores in raw_read_
          consensus_scores.push_back(0);
        }

        consensus.insert(consensus.end(), chunk.consensus.begin(),
                         chunk.consensus.end());

        consensus_scores.insert(consensus_scores.end(),
                                chunk.consensus_scores.begin(),
                                chunk.consensus_scores.end());

      } else {
        int overlap = right_end_offset - chunk.approx_read_offset;
        overlap -= overlapAlignH(consensus.end() - overlap, consensus.end(),
                                 chunk.consensus.begin(),
                                 chunk.consensus.begin() + overlap, 5);

        if (overlap > static_cast<int>(chunk.consensus.size())) return false;

        if (overlap < 0) {
          chunk.approx_read_offset = right_end_offset - overlap;
          return DoMerge(chunk);
        }

        int n_trim = 0;
        int n_runs = int(consensus.size());

        // FIXME
        if (overlap > 0 && rollback_end > 0) {
          for (int i = 0; i < overlap; i++) {
            if (n_runs - overlap + i < 0 ||
                n_runs - overlap + i >= consensus.size())
              continue;
            auto left_run = consensus[n_runs - overlap + i];
            auto right_run = chunk.consensus[i];
            if (left_run != right_run) {
              RollBack();
              AlignRightEndAgainstRead();
              return DoMerge(chunk);
            }
          }
        }

        if (overlap >= 3 && n_runs > overlap) {
          for (; n_trim < overlap / 3; ++n_trim) {
            auto score1 = consensus_scores[n_runs - n_trim - 1];
            auto score2 = chunk.consensus_scores[overlap - n_trim - 1];
            if (score1 > score2) break;
          }

          consensus.resize(consensus.size() - n_trim);
          consensus_scores.resize(consensus_scores.size() - n_trim);
        }

        consensus.insert(consensus.end(),
                         chunk.consensus.begin() + overlap - n_trim,
                         chunk.consensus.end());

        consensus_scores.insert(
            consensus_scores.end(),
            chunk.consensus_scores.begin() + overlap - n_trim,
            chunk.consensus_scores.end());
      }

      approx_end_read_offset_ = chunk.approx_end_read_offset();
      return true;
    }

    bool MergeWithDisjointChunk(ConsensusChunk &chunk) {
      if (debug_mode) std::cerr << "[MergeWithDisjointChunk]" << std::endl;
      AlignRightEndAgainstRead();
      if (chunk.alignment != kChunkLeftAligned) chunk.AlignLeftEndAgainstRead();
      return DoMerge(chunk);
    }

    bool MergeWithOverlappingChunk(ConsensusChunk &chunk) {
      if (debug_mode) std::cerr << "[MergeWithOverlappingChunk]" << std::endl;
      int right_end_offset = approx_end_read_offset_;
      size_t overlap = right_end_offset - chunk.approx_read_offset;
      if (overlap > chunk.consensus_scores.size()) return false;

      AlignRightEndAgainstRead();
      if (chunk.alignment != kChunkLeftAligned) chunk.AlignLeftEndAgainstRead();
      return DoMerge(chunk);
    }

   public:
    bool TryMergeWith(ConsensusChunk &chunk) {
      if (chunk.consensus.empty()) return true;

      alignment = kChunkNotAligned;
      int right_end_offset = approx_end_read_offset_;

      if (right_end_offset <= chunk.approx_read_offset)
        return MergeWithDisjointChunk(chunk);
      else
        return MergeWithOverlappingChunk(chunk);
    }
  };

  // Chunks where strong consensus was obtained
  std::list<ConsensusChunk> chunks_;
  int trimmed_by_gen_;

  void PushChunk(const ScoreStorage &scores, int initial_read_offset,
                 int approx_read_offset, int approx_end_read_offset,
                 unsigned rollback_end) {
    chunks_.push_back(ConsensusChunk(initial_read_offset, approx_read_offset,
                                     approx_end_read_offset, scores,
                                     rollback_end, raw_read_, debug_mode_));
    if (debug_mode_) {
      auto &consensus = chunks_.back().consensus;
      size_t len = consensus.size();
      size_t nucl_len = 0;
      for (size_t i = 0; i < len; ++i) nucl_len += consensus[i].len;
    }

    chunks_.back().AlignLeftEndAgainstRead();
    if (chunks_.size() == 1)
      trimmed_by_gen_ = chunks_.back().raw_start_offset();
  }

  const ConsensusChunk &LastChunk() const { return chunks_.back(); }

  class ChunkCollector {
    CorrectedRead &cread_;
    const KMerData &kmer_data_;
    bool debug_mode_;

    ValidHKMerGenerator<hammer::K> gen;
    int pos;
    unsigned skipped;
    int raw_pos;

    struct Center {
      hammer::HKMer seq;
      int end_offset;
    };

    Center last_good_center;
    bool last_good_center_is_defined;
    bool is_first_center;
    bool replacing;
    int rollback_size;

    bool need_to_align;

    int approx_read_offset;
    int approx_end_read_offset;
    ScoreStorage scores;
    int chunk_pos;
    int raw_chunk_start_pos;

    unsigned approx_n_insertions;

    Center GetCenterOfCluster(const hammer::HKMer &seq, int start_pos) const {
      hammer::KMerStat k[2];
      k[0] = kmer_data_[kmer_data_[seq].changeto];
      k[1] = kmer_data_[kmer_data_[!seq].changeto];
      k[1].kmer = !k[1].kmer;

      if (k[0].qual > k[1].qual) std::swap(k[0], k[1]);
      using namespace hammer;
      for (size_t i = 0; i < 2; ++i) {
        auto &kmer = k[i].kmer;
        int end_diff;
        auto dist = distanceHKMer(kmer.begin(), kmer.end(), seq.begin(),
                                  seq.end(), 3, &end_diff);
        if (debug_mode_) {
          std::cerr << "[GetCenterOfCluster] distance(" << seq << ", " << kmer
                    << ") = " << dist << std::endl;
        }
        if (dist <= 2) {
          return Center{kmer, start_pos + int(hammer::K) + end_diff};
        }
      }
      return Center{seq, start_pos + int(hammer::K)};
    }

    bool IsInconsistent(const Center &center) const {
      if (!last_good_center_is_defined) return false;

      for (size_t i = 0; i < hammer::K - skipped - 1; ++i)
        if (last_good_center.seq[i + skipped + 1].nucl != center.seq[i].nucl)
          return true;

      return false;
    }

    void FlushCurrentChunk() {
      unsigned rollback_end = 0;

      if (replacing) {
        if (rollback_size < 0) rollback_size = 0;
        if (rollback_size < int(scores.size()))
          rollback_end = int(scores.size()) - rollback_size;
        replacing = false;
        rollback_size = 0;
      }

      if (scores.size() > hammer::K) {
        cread_.PushChunk(scores, raw_chunk_start_pos, approx_read_offset,
                         approx_end_read_offset, rollback_end);
        pos = cread_.LastChunk().approx_end_read_offset_untrimmed() - hammer::K;
        pos += skipped;
      } else {
        pos -= approx_n_insertions;
      }

      scores.clear();
      need_to_align = false;
      chunk_pos = 0;
      skipped = 0;
      approx_n_insertions = 0;
      approx_read_offset = pos;

      last_good_center_is_defined = false;
    }

    // side effect: changes chunk_pos, pos, and approx_n_insertions
    bool TryToAlignCurrentCenter(const Center &center) {
      if (!last_good_center_is_defined) return true;

      if (debug_mode_) {
        std::cerr << "[TryToAlignCurrentCenter] " << center.seq.str()
                  << " (previous good center is " << last_good_center.seq.str()
                  << ","
                  << " skipped " << skipped << " centers)" << std::endl;
      }

      // offset is how many positions the center should be shifted
      // in order to agree with last_good_center
      int offset;
      bool aligned = exactAlignH(last_good_center.seq.begin(),
                                 last_good_center.seq.begin() + skipped + 1,
                                 last_good_center.seq.end(), center.seq.begin(),
                                 center.seq.end(), 3, 8, &offset);

      bool result = aligned && chunk_pos + offset >= 0;
      if (result) {
        if (debug_mode_)
          std::cerr << "[TryToAlignCurrentCenter] offset = " << offset
                    << std::endl;
        if (offset < 0) approx_n_insertions -= offset;
        pos += offset;
        chunk_pos += offset;
      }

      return result;
    }

    void IncludeIntoConsensus(const Center &center) {
      VERIFY(chunk_pos >= 0);
      VERIFY(chunk_pos < (1 << 16));
      is_first_center = false;

      if (chunk_pos + hammer::K > scores.size())
        scores.resize(chunk_pos + hammer::K, ScoreMatrix(4, 64, 0));

      auto k = kmer_data_[center.seq];

      for (size_t i = 0; i < hammer::K; ++i)
        scores[chunk_pos + i](center.seq[i].nucl, center.seq[i].len) +=
            double(k.count) * (1.0 - k.qual);

      last_good_center = center;
      last_good_center_is_defined = true;
      if (raw_chunk_start_pos == -1) raw_chunk_start_pos = raw_pos;
      approx_end_read_offset = center.end_offset;
      if (debug_mode_) {
        std::cerr << "e.o. = " << approx_end_read_offset << std::endl;
      }
      need_to_align = false;
      skipped = 0;
    }

   public:
    ChunkCollector(const io::SingleRead &r, CorrectedRead &cread,
                   const KMerData &kmer_data, bool debug_mode)
        : cread_(cread),
          kmer_data_(kmer_data),
          debug_mode_(debug_mode),
          gen(r),
          pos(int(gen.trimmed_left())),
          skipped(0),
          last_good_center(),
          last_good_center_is_defined(false),
          is_first_center(true),
          replacing(false),
          rollback_size(0),
          need_to_align(false),
          approx_read_offset(0),
          approx_end_read_offset(0),
          scores(),
          chunk_pos(0),
          raw_chunk_start_pos(-1),
          approx_n_insertions(0) {
      --pos;
      --chunk_pos;
    }

    void Run() {
      double lowQualThreshold = cfg::get().kmer_qual_threshold;

      raw_pos = int(gen.trimmed_left()) - 1;

      if (debug_mode_) {
        std::cerr << "gen. trimmed = " << gen.trimmed_left() << std::endl;
      }

      while (gen.HasMore()) {
        auto prev_chunk_pos = chunk_pos;
        auto seq = gen.kmer();
        gen.Next();
        ++pos;
        ++raw_pos;
        if (debug_mode_) {
          std::cerr << "=================================" << std::endl;
          std::cerr << "pos = " << pos << ", raw_pos = " << raw_pos
                    << ", last_good_center_is_defined = "
                    << last_good_center_is_defined << ", skipped = " << skipped
                    << std::endl;
        }
        ++chunk_pos;

        auto center = Center{seq, raw_pos + int(hammer::K)};
        auto qual = kmer_data_[seq].qual;

        bool can_be_changed = last_good_center_is_defined || is_first_center;
        if (qual > lowQualThreshold && can_be_changed) {
          center = GetCenterOfCluster(seq, raw_pos);
          qual = kmer_data_[center.seq].qual;
        }

        if (qual > lowQualThreshold && last_good_center_is_defined &&
            skipped == 0) {
          if (debug_mode_) {
            std::cerr << "raw_pos + hammer::K = " << raw_pos + hammer::K
                      << std::endl;
            std::cerr << "last_good_center.end_offset + 1 = "
                      << last_good_center.end_offset + 1 << std::endl;
          }
          // Finding a center by means of clustering failed.
          // Let's try the following: take last good center and make a new one
          // from it by appending next homopolymer run; if its quality is high,
          // we use it.
          if (raw_pos + hammer::K < last_good_center.end_offset + 1) {
            --pos;
            --chunk_pos;
            if (debug_mode_) {
              std::cerr << "skipping low-quality hk-mer" << std::endl;
            }
            continue;  // move to next hk-mer
          } else if (raw_pos + hammer::K == last_good_center.end_offset + 1) {
            auto seq_corr = last_good_center.seq;
            for (size_t i = 0; i < hammer::K - 1; ++i)
              seq_corr[i] = seq_corr[i + 1];
            seq_corr[hammer::K - 1] = seq[hammer::K - 1];
            center = Center{seq_corr, last_good_center.end_offset + 1};
            qual = kmer_data_[center.seq].qual;
            if (debug_mode_) {
              std::cerr << "seq_corr = " << seq_corr.str()
                        << " , qual = " << qual << std::endl;
            }

            if (qual > lowQualThreshold && can_be_changed) {
              // our last resort...
              center = GetCenterOfCluster(seq_corr, raw_pos);
              qual = kmer_data_[center.seq].qual;
            }
          }
        }

        bool low_qual = qual > lowQualThreshold;
        bool inconsistent = IsInconsistent(center);

        if (debug_mode_ && !low_qual && seq != center.seq) {
          std::cerr << "replaced " << seq.str() << " (quality "
                    << kmer_data_[seq].qual << ", count "
                    << kmer_data_[seq].count << ")"
                    << " with " << center.seq.str() << std::endl;
        }

        if (debug_mode_) {
          std::cerr << "quality of " << center.seq.str() << " is " << qual
                    << " (count " << kmer_data_[center.seq].count << ") "
                    << (inconsistent ? " INCONSISTENT" : "") << std::endl;
        }

        if (low_qual) {
          ++skipped;
        } else if (inconsistent) {
          if (!TryToAlignCurrentCenter(center)) {
            low_qual = true;
            ++skipped;
          }
        }

        if (skipped > hammer::K / 4) {
          FlushCurrentChunk();
        } else if (!low_qual) {
          if (seq != center.seq && !replacing) {
            rollback_size = prev_chunk_pos + hammer::K;
            replacing = true;
          } else if (seq == center.seq && replacing) {
            replacing = false;
          }

          if (debug_mode_) {
            std::cerr << "[include into consensus] raw_pos = " << raw_pos
                      << std::endl;
          }
          IncludeIntoConsensus(center);
        }
      }

      FlushCurrentChunk();
    }
  };

  void CollectChunks(const io::SingleRead &r) {
    ChunkCollector chunk_collector(r, *this, kmer_data_, debug_mode_);
    chunk_collector.Run();
  }

 public:
  CorrectedRead(const io::SingleRead &read, const KMerData &kmer_data,
                bool debug_mode = false)
      : raw_read_(read), kmer_data_(kmer_data), debug_mode_(debug_mode) {
    CollectChunks(read);
  }

  void MergeChunks() {
    if (chunks_.empty()) return;

    auto iter = chunks_.begin();
    ConsensusChunk &merged = *iter;

    if (debug_mode_) {
      if (chunks_.size() == 1) {
        iter->AlignLeftEndAgainstRead();
        for (int i = 0; i < iter->approx_read_offset; ++i)
          for (int j = 0; j < raw_read_[i].len; ++j) std::cerr << ' ';
        for (size_t i = 0; i < iter->consensus.size(); ++i)
          std::cerr << iter->consensus[i].str();
        std::cerr << std::endl;
      }
    }

    ++iter;
    while (iter != chunks_.end()) {
      if (iter->consensus.size() > hammer::K) merged.TryMergeWith(*iter);
      iter = chunks_.erase(iter);
    }

    corrected_runs_ = std::move(merged.consensus);
  }

  void AttachUncorrectedRuns() {
    // attach runs from the right
    const auto &data = raw_read_.data();
    int n_raw = int(raw_read_.size());
    int end_read_offset = LastChunk().approx_end_read_offset();
    if (end_read_offset < n_raw && end_read_offset >= 0) {
      corrected_runs_.insert(corrected_runs_.end(),
                             data.begin() + end_read_offset, data.end());
    }
    if (debug_mode_) {
      std::cerr << "n_raw = " << n_raw
                << ", end_read_offset = " << end_read_offset << std::endl;
    }

    // attach runs from the left
    if (trimmed_by_gen_ > 0 && size_t(trimmed_by_gen_) <= data.size()) {
      std::vector<HomopolymerRun> runs;
      runs.reserve(corrected_runs_.size() + trimmed_by_gen_);
      runs.insert(runs.end(), data.begin(), data.begin() + trimmed_by_gen_);
      runs.insert(runs.end(), corrected_runs_.begin(), corrected_runs_.end());
      std::swap(runs, corrected_runs_);
    }
  }

  std::string GetSequenceString() const {
    if (chunks_.empty() && corrected_runs_.empty()) return "";
    std::string res;
    if (!corrected_runs_.empty()) {
      for (auto it = corrected_runs_.begin(); it != corrected_runs_.end(); ++it)
        res += it->str();
    } else {
      auto &runs = chunks_.front().consensus;
      for (auto it = runs.begin(); it != runs.end(); ++it) res += it->str();
    }
    return res;
  }
};

class SingleReadCorrector {
  const KMerData &kmer_data_;

 public:
  struct ReadSelectionPredicate {
    virtual bool operator()(const io::SingleRead &read) = 0;
  };

  struct DebugOutputPredicate : public ReadSelectionPredicate {};

  struct NoDebug : public DebugOutputPredicate {
    virtual bool operator()(const io::SingleRead &) { return false; }
  };

  struct FullDebug : public DebugOutputPredicate {
    virtual bool operator()(const io::SingleRead &) { return true; }
  };

  class DebugIfContains : public DebugOutputPredicate {
    Sequence needle_;
    Sequence needle_rc_;

   public:
    DebugIfContains(const Sequence &seq) : needle_(seq), needle_rc_(!seq) {}

    virtual bool operator()(const io::SingleRead &read) {
      auto read_seq = read.sequence();
      if (read_seq.size() < needle_.size()) return false;
      if (read_seq.find(needle_, 0) != -1ULL) return true;
      if (read_seq.find(needle_rc_, 0) != -1ULL) return true;
      return false;
    }
  };

  struct SelectPredicate : public ReadSelectionPredicate {};
  struct SelectAll : public SelectPredicate {
    virtual bool operator()(const io::SingleRead &) { return true; }
  };

  class SelectByName : public SelectPredicate {
    std::set<std::string> names_;

   public:
    SelectByName(const std::set<std::string> &names) : names_(names) {}
    virtual bool operator()(const io::SingleRead &r) {
      return names_.find(r.name()) != names_.end();
    }
  };

 private:
  BamTools::SamHeader *sam_header_;
  DebugOutputPredicate &debug_pred_;
  SelectPredicate &select_pred_;

 public:
  SingleReadCorrector(const KMerData &kmer_data,
                      BamTools::SamHeader *sam_header,
                      DebugOutputPredicate &debug, SelectPredicate &select)
      : kmer_data_(kmer_data),
        sam_header_(sam_header),
        debug_pred_(debug),
        select_pred_(select) {}

  SingleReadCorrector(const KMerData &kmer_data, DebugOutputPredicate &debug,
                      SelectPredicate &select)
      : kmer_data_(kmer_data),
        sam_header_(NULL),
        debug_pred_(debug),
        select_pred_(select) {}

  std::unique_ptr<io::SingleRead> operator()(
      std::unique_ptr<io::SingleRead> r) {
    return operator()(*r);
  }

  std::unique_ptr<io::SingleRead> operator()(const io::SingleRead &r) {
    if (!select_pred_(r)) return nullptr;
    bool debug_mode = debug_pred_(r);
    if (debug_mode) {
      std::cerr << "=============================================" << std::endl;

      std::cerr << '>' << r.name() << '\n'
                << r.GetSequenceString() << std::endl;
    }

    CorrectedRead read(r, kmer_data_, debug_mode);
    read.MergeChunks();
    if (cfg::get().keep_uncorrected_ends) read.AttachUncorrectedRuns();

    if (debug_mode) {
      std::cerr << "final result: " << read.GetSequenceString() << std::endl;
    }

    auto seq = read.GetSequenceString();
    if (seq.empty()) return nullptr;

    return std::unique_ptr<io::SingleRead>(new io::SingleRead(r.name(), seq));
  }

  std::unique_ptr<io::BamRead> operator()(
      std::unique_ptr<BamTools::BamAlignment> alignment) {
    VERIFY(sam_header_);
    io::SingleRead r(alignment->Name, alignment->QueryBases);
    // reverse strand means we're working with a mapped BAM, might be
    // the case for datasets downloaded from IonCommunity
    if (alignment->IsReverseStrand()) r = !r;
    auto corrected_r = operator()(r);
    std::string rg;
    if (!alignment->GetTag("RG", rg) || !corrected_r) return nullptr;
    auto flow_order = sam_header_->ReadGroups[rg].FlowOrder;

    float delta_score, fit_score;
    auto seq = corrected_r->GetSequenceString();
    if (alignment->IsReverseStrand()) {
      std::reverse(seq.begin(), seq.end());
      for (auto it = seq.begin(); it != seq.end(); ++it) {
        switch (*it) {
          case 'A':
            *it = 'T';
            break;
          case 'C':
            *it = 'G';
            break;
          case 'G':
            *it = 'C';
            break;
          case 'T':
            *it = 'A';
            break;
          default:
            break;
        }
      }
    }

    BaseHypothesisEvaluator(*alignment, flow_order, seq, delta_score, fit_score,
                            0);
    std::stringstream ss;
    ss << alignment->Name << "_" << delta_score << "_" << fit_score;
    alignment->Name = ss.str();
    if (delta_score >= cfg::get().delta_score_threshold)
      return std::unique_ptr<io::BamRead>(new io::BamRead(*alignment));

    BamTools::BamAlignment corrected(*alignment);
    corrected.QueryBases = corrected_r->GetSequenceString();
    return std::unique_ptr<io::BamRead>(new io::BamRead(corrected));
  }
};

class PairedReadCorrector : public SingleReadCorrector {
 public:
  PairedReadCorrector(const KMerData &kmer_data, DebugOutputPredicate &debug,
                      SelectPredicate &select)
      : SingleReadCorrector(kmer_data, debug, select) {}

  std::unique_ptr<io::PairedRead> operator()(
      std::unique_ptr<io::PairedRead> r) {
    auto corrected_r = SingleReadCorrector::operator()(r->first());
    auto corrected_l = SingleReadCorrector::operator()(r->second());

    if (!corrected_r || !corrected_l) return nullptr;

    return std::unique_ptr<io::PairedRead>(
        new io::PairedRead(*corrected_r, *corrected_l, 0));
  }
};

};      // namespace correction
};      // namespace hammer
#endif  // __HAMMER_IT_READ_CORRECTOR_HPP__
