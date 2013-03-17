#ifndef __HAMMER_IT_READ_CORRECTOR_HPP__
#define __HAMMER_IT_READ_CORRECTOR_HPP__

#include "HSeq.hpp"
#include "flow_space_read.hpp"
#include "hkmer_distance.hpp"
#include "consensus.hpp"
#include "valid_hkmer_generator.hpp"
#include "io/single_read.hpp"

#include <boost/optional.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <deque>
#include <vector>
#include <cassert>
#include <list>
#include <string>
#include <algorithm>

namespace hammer {
namespace correction {

namespace numeric = boost::numeric::ublas;

typedef numeric::matrix<double> ScoreMatrix;
typedef std::vector<ScoreMatrix> ScoreStorage;

template <typename It1, typename It2>
static bool exactAlignH(It1 a_begin, It1 a_initial_pos, It1 a_end,
                        It2 b_initial_pos, It2 b_end,
                        size_t max_offset, size_t n_cmp, int* p_offset)
{
  int M = max_offset * 2 + 1;
  for (int i = 0; i <= M; i++) {
    int offset = (i / 2) * ((i & 1) ? 1 : -1); // 0, -1, 1, -2, 2, ...
    auto a_it = a_initial_pos + offset;
    auto b_it = b_initial_pos;
    if (a_it < a_begin || a_it + n_cmp > a_end)
      continue;
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
static int alignH(It1 read_begin, It1 read_end,
                  It2 consensus_begin, It2 consensus_end,
                  int approx_read_offset, size_t n_skip_consensus,
                  size_t n_side = 5) {

  const int kRunsToCompare = 8;

  int left_offset = n_side;

  auto x_begin = read_begin + approx_read_offset - n_side;
  if (x_begin + kRunsToCompare >= read_end) {
    x_begin = read_end - kRunsToCompare - 2 * n_side;
    left_offset = read_begin + approx_read_offset - x_begin;
  }

  if (x_begin < read_begin) {
    x_begin = read_begin;
    left_offset = approx_read_offset;
  }
 
  auto x_end = x_begin + 2 * n_side + kRunsToCompare;
  if (x_end > read_end)
    x_end = read_end;

  auto y_begin = consensus_begin + n_skip_consensus;
  auto y_end = y_begin + kRunsToCompare;
  if (y_end > consensus_end)
    y_end = consensus_end;

  // Smith-Waterman for homopolymer runs
  const int kDirOffset = 16;
  const int kScoreMask = (1 << kDirOffset) - 1;

  const int kDirUpLeft = 0;
  const int kDirUp = 1;
  const int kDirLeft = 2;

  const int kBaseDiff = -1;
  const int kRunInsertion = -5;
  const int kRunDeletion = -5;
  const int kMismatch = -4;
  const int kMatch = 3;

  size_t m = x_end - x_begin;
  size_t n = y_end - y_begin;

  using namespace boost::numeric::ublas;
  matrix<int> scores(m + 1, n + 1, 0);

  size_t highest_x = 0, highest_y = 0;
  int highest_entry = 0;

  for (size_t i = 1; i <= m; i++) {
    for (size_t j = 1; j <= n; j++) {
      int best_score = 0;
      int best_dir = 0;

      auto run_x = *(x_begin + i - 1);
      auto run_y = *(y_begin + j - 1);

      int score;
      if (run_x.raw == run_y.raw) {
        score = kMatch + (scores(i - 1, j - 1) & kScoreMask);
        if (score > best_score) {
          best_score = score;
          best_dir = kDirUpLeft;
        }
      } else if (run_x.nucl == run_y.nucl) {
        score = kBaseDiff * std::abs(run_x.len - run_y.len);
        score += scores(i - 1, j - 1) & kScoreMask;
        if (score > best_score) {
          best_score = score;
          best_dir = kDirUpLeft;
        }
      } else {
        score = kMismatch + (scores(i - 1, j - 1) & kScoreMask);
        if (score > best_score) {
          best_score = score;
          best_dir = kDirUpLeft;
        }
      }
 
      score = (scores(i - 1, j) & kScoreMask) + kRunDeletion;
      if (score > best_score) {
        best_score = score;
        best_dir = kDirUp;
      }

      score = (scores(i, j - 1) & kScoreMask) + kRunInsertion;
      if (score > best_score) {
        best_score = score;
        best_dir = kDirLeft;
      }

      scores(i, j) = best_score + (best_dir << kDirOffset);

      if (best_score > highest_entry) {
        highest_entry = best_score;
        highest_x = i;
        highest_y = j;
      }
    }
  }

  // if score is low, search on larger portion of the read
  if (highest_entry < (kRunsToCompare * kMatch * 3) / 4 &&
      (x_begin != read_begin || x_end != read_end))
    return alignH(read_begin, read_end,
                  consensus_begin, consensus_end,
                  approx_read_offset, n_skip_consensus,
                  n_side * 3);

  int x = highest_x;
  int y = highest_y;
  while ((x >= 0 || y >= 0) && (scores(x, y) & kScoreMask) != 0) {
    int dir = scores(x, y) >> kDirOffset;
    switch (dir) {
      case kDirUp:
        --x; break;
      case kDirLeft:
        --y; break;
      case kDirUpLeft:
        --x, --y; break;
      default:
        break;
    }
  }

  return x - y - left_offset;
}

// Not used now
class HKMerProlonger {
  const KMerData& kmer_data_;

 public:
  struct RightSide {
    static size_t changingPosition() { return hammer::K - 1; }
    static hammer::HKMer shift(const hammer::HKMer &kmer) {
      hammer::HKMer res;
      for (size_t i = 1; i < hammer::K; ++i)
        res[i - 1] = kmer[i];
      return res;
    }
    template <typename T, typename U>
    static void append(T& cont, U obj) { cont.push_back(obj); }
  }; 

  struct LeftSide {
    static size_t changingPosition() { return 0; }
    static hammer::HKMer shift(const hammer::HKMer &kmer) {
      hammer::HKMer res;
      for (size_t i = 1; i < hammer::K; ++i)
        res[i] = kmer[i - 1];
      return res;
    }
    template <typename T, typename U>
    static void append(T& cont, U obj) { cont.push_front(obj); }
  };

 public:

  /// @param[in] seed               kmer to prolong
  /// @param[in] bases_to_recover   maximum number of bases to recover
  /// @param[in] side               side to prolong to (RightSide/LeftSide)
  template <typename Side>
  std::deque<hammer::HomopolymerRun> prolong(const hammer::HKMer &seed, 
                                             size_t bases_to_recover,
                                             Side side) {
    std::deque<hammer::HomopolymerRun> good_runs(hammer::K);
    for (size_t i = 0; i < hammer::K; ++i)
      good_runs[i] = seed[i];

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
        if (nucl == good_kmer[changing_pos].nucl)
          continue;
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
  HKMerProlonger(const KMerData& kmer_data) : kmer_data_(kmer_data) {}
};

static const double kLowScoreThreshold = 1.0;

class CorrectedRead {
  FlowSpaceRead raw_read_; // Uncorrected read
  const KMerData& kmer_data_;

  // Stores runs after joining chunks
  std::vector<hammer::HomopolymerRun> corrected_runs_;

  // Contiguous part of read with strong consensus
  struct ConsensusChunk {
    int approx_read_offset; // in the vector of raw read runs

    enum {
      kChunkLeftAligned,
      kChunkRightAligned,
      kChunkNotAligned
    } alignment;

    const FlowSpaceRead& raw_read;
    size_t trimmed_left;
    size_t trimmed_right;

    std::vector<hammer::HomopolymerRun> consensus;
    std::vector<double> consensus_scores;

    ConsensusChunk(int approximate_read_offset,
                   const ScoreStorage &scores,
                   const FlowSpaceRead &read,
                   double trim_threshold=1.0)
        : approx_read_offset(approximate_read_offset), 
          alignment(kChunkNotAligned), raw_read(read),
          trimmed_left(0), trimmed_right(0)
    {
      bool left_trim = true;
      for (size_t i = 0; i < scores.size(); ++i) {
        auto run = hammer::iontorrent::consensus(scores[i]);

        // trim low-quality runs from the left side
        if (run.second <= kLowScoreThreshold && left_trim) {
          approx_read_offset += 1;
          trimmed_left += 1;
          continue;
        }

        left_trim = false;
        VERIFY(run.first.len > 0);
        consensus.push_back(run.first);
        consensus_scores.push_back(run.second);
      }

      size_t right_end = consensus_scores.size();
      if (right_end == 0)
        return;

      while (consensus_scores[right_end - 1] <= kLowScoreThreshold) {
        --right_end;
        if (right_end == 0)
          break;
      }

      trimmed_right = consensus.size() - right_end;
      consensus.resize(right_end);
      consensus_scores.resize(right_end);
    }

    void AlignLeftEndAgainstRead(size_t skip=0) {
      const auto& data = raw_read.data();

#if 0
      std::cerr << "[approx. read offset (left)] before: " << approx_read_offset << "; after: ";
#endif
      int offset = alignH(data.begin(), data.end(),
                          consensus.begin(), consensus.end(),
                          approx_read_offset, skip);

      approx_read_offset += offset;
#if 0
      std::cerr << approx_read_offset << std::endl;
#endif
      alignment = kChunkLeftAligned;
    }

    void AlignRightEndAgainstRead(size_t skip=0) {
      const auto& data = raw_read.data();
#if 0
      std::cerr << "[approx. read offset (right)] before: " << approx_read_offset << "; after: ";
#endif
      int position_on_read = approx_read_offset + consensus.size() - 1;
      int offset = alignH(data.rbegin(), data.rend(),
                          consensus.rbegin(), consensus.rend(),
                          data.size() - 1 - position_on_read, skip);

      approx_read_offset -= offset;
#if 0
      std::cerr << approx_read_offset << std::endl;
#endif
      alignment = kChunkRightAligned;
    }

    int approx_end_read_offset() const {
      return approx_read_offset + consensus.size();
    }
 
    int approx_end_read_offset_untrimmed() const {
      return approx_end_read_offset() + trimmed_right;
    }


   private:
    bool DoMerge(ConsensusChunk& chunk) {
      int right_end_offset = approx_read_offset + consensus.size();

#if 0
      std::cerr << "============== Merging chunks ===============" << std::endl;

      int white_l = 0;
      for (int i = right_end_offset - 1; i >= 0; --i)
        white_l += raw_read[i].len;
      for (int i = 0; i < consensus.size(); ++i)
        white_l -= consensus[i].len;
      for (int i = 0; i < white_l; ++i)
        std::cerr << ' ';
      for (int i = std::max(-white_l, 0); i < consensus.size(); ++i)
        std::cerr << consensus[i].str();
      std::cerr << std::endl;

      for (int i = 0; i < chunk.approx_read_offset; ++i)
        for (int j = 0; j < raw_read[i].len; ++j)
          std::cerr << ' ';
      for (int i = 0; i < chunk.consensus.size(); ++i)
        std::cerr << chunk.consensus[i].str();
      std::cerr << std::endl;
#endif

      if (right_end_offset <= chunk.approx_read_offset) {

        for (int i = right_end_offset; i < chunk.approx_read_offset; ++i) {
          if (i >= static_cast<int>(raw_read.size()))
            return false;
          consensus.push_back(raw_read[i]);
          alignment = kChunkNotAligned;

          // TODO: maintain quality scores in raw_read_
          consensus_scores.push_back(0); 
        }

        consensus.insert(consensus.end(), 
                         chunk.consensus.begin(), chunk.consensus.end());

        consensus_scores.insert(consensus_scores.end(),
                                chunk.consensus_scores.begin(),
                                chunk.consensus_scores.end());

      } else {
        int overlap = right_end_offset - chunk.approx_read_offset;
        if (overlap > static_cast<int>(chunk.consensus.size()))
          return false;

        consensus.insert(consensus.end(),
                         chunk.consensus.begin() + overlap,
                         chunk.consensus.end());

        consensus_scores.insert(consensus_scores.end(),
                                chunk.consensus_scores.begin() + overlap,
                                chunk.consensus_scores.end());
      }

      return true;
    }

    bool MergeWithDisjointChunk(ConsensusChunk& chunk) {
      AlignRightEndAgainstRead();
      if (chunk.alignment != kChunkLeftAligned)
        chunk.AlignLeftEndAgainstRead();
      return DoMerge(chunk);
    }

    bool MergeWithOverlappingChunk(ConsensusChunk& chunk) {
      int right_end_offset = approx_read_offset + consensus.size();
      size_t overlap = right_end_offset - chunk.approx_read_offset;
      if (overlap > chunk.consensus_scores.size())
        return false;

      // FIXME
      size_t trim_r = 0, trim_l = 0;
      /*
      // remove low scored bases from the two ends
      while (overlap != 0) {
        size_t pos_r = consensus.size() - 1 - trim_r;
        size_t pos_l = chunk.consensus.size() + trim_l;
        if (pos_r == 0 || pos_l == chunk.consensus.size())
          return false;

        double score_r = consensus_scores[pos_r];
        double score_l = chunk.consensus_scores[pos_l];
        if (score_r > score_l)
          ++trim_l;
        else
          ++trim_r;

        if (std::min(score_r, score_l) / std::max(score_r, score_l) >= 0.1)
          break;

        --overlap;
      }*/

      AlignRightEndAgainstRead(trim_r);
      chunk.AlignLeftEndAgainstRead(trim_l);
      return DoMerge(chunk);
    }

   public:

    bool TryMergeWith(ConsensusChunk& chunk) {
      if (chunk.consensus.empty())
        return true;

      alignment = kChunkNotAligned;
      int right_end_offset = approx_read_offset + consensus.size();

      if (right_end_offset <= chunk.approx_read_offset)
        return MergeWithDisjointChunk(chunk);
      else
        return MergeWithOverlappingChunk(chunk);
    }

  };

  // Chunks where strong consensus was obtained
  std::list<ConsensusChunk> chunks_;

  void PushChunk(const ScoreStorage &scores, int approx_read_offset) {
    chunks_.push_back(ConsensusChunk(approx_read_offset, scores, raw_read_));
    chunks_.back().AlignLeftEndAgainstRead();
  }

  const ConsensusChunk& LastChunk() const {
    return chunks_.back();
  }

  void CollectChunks(const io::SingleRead& r) {
    ScoreStorage scores;

    ValidHKMerGenerator<hammer::K> gen(r);
    size_t pos = gen.trimmed_left(); // current position
    size_t chunk_pos = 0;
    unsigned skipped = 0; // number of skipped centers of low quality
    hammer::HKMer last_good_center;
    bool last_good_center_is_defined = false;
    bool start_new_chunk = true;
    bool need_to_align = false;
    int approx_read_offset = 0;

    // approximate number of insertions in the current chunk
    unsigned approx_n_insertions = 0;

    while (gen.HasMore()) {
      hammer::HKMer seq = gen.kmer();
      gen.Next();

      hammer::KMerStat k = kmer_data_[kmer_data_[seq].changeto];
      hammer::HKMer center = k.kmer;

      const double LOW_QUALITY_THRESHOLD = 1e-11;

      bool low_qual = k.qual > LOW_QUALITY_THRESHOLD;

      // if too many centers are skipped, start new chunk
      if (skipped > 4 && !start_new_chunk) {
        PushChunk(scores, approx_read_offset);
        start_new_chunk = true;
        pos = LastChunk().approx_end_read_offset_untrimmed() - hammer::K;
        pos += skipped + 1;
        last_good_center_is_defined = false;
      } else {
        // otherwise, compare current center with the last 'good' center.
        // if they differ too much or the center has low quality, skip it.
        if (last_good_center_is_defined && !low_qual) {
          for (size_t i = 0; i < hammer::K - skipped - 1; ++i) {
            if (last_good_center[i + skipped + 1].nucl != center[i].nucl) {
              low_qual = true;
              break;
            }
          }
        }
      }

      if (low_qual) {
        need_to_align = true;
        skipped += 1;
        goto process_next_kmer;
      }

      if (need_to_align && last_good_center_is_defined) {
        VERIFY(skipped + 1 < hammer::K);
        int offset;
        bool aligned = exactAlignH(last_good_center.begin(),
                                   last_good_center.begin() + skipped + 1,
                                   last_good_center.end(),
                                   center.begin(), center.end(), 3, 5, &offset);
        if (!aligned) {
          if (!start_new_chunk) {
            PushChunk(scores, approx_read_offset);
            pos = LastChunk().approx_end_read_offset_untrimmed() - hammer::K;
            pos += skipped + 1;
          }
          start_new_chunk = true;
        } else {
          if (offset < 0)
            approx_n_insertions -= offset;
          pos += skipped + offset;
          chunk_pos += skipped + offset;
        }
      }

      if (start_new_chunk) {
        chunk_pos = 0;
        approx_read_offset = pos + approx_n_insertions;
        scores.clear();
        approx_n_insertions = 0;
        start_new_chunk = false;
      }

      if (chunk_pos + hammer::K > scores.size())
        scores.resize(chunk_pos + hammer::K, ScoreMatrix(4, 64, 0));

      for (size_t i = 0; i < hammer::K; ++i)
        scores[chunk_pos + i](center[i].nucl, center[i].len) += k.count * (1.0 - k.qual);

      last_good_center = center;
      last_good_center_is_defined = true;
      need_to_align = false;
      skipped = 0;
      ++pos;
      ++chunk_pos;

process_next_kmer: ;
    }

    if (!start_new_chunk)
      PushChunk(scores, approx_read_offset);
  }
  
 public:
  CorrectedRead(const io::SingleRead& read, const KMerData& kmer_data) :
    raw_read_(read), 
    kmer_data_(kmer_data)
  {
    CollectChunks(read);
  }

  void MergeChunks() {
    if (chunks_.empty())
      return;

    auto iter = chunks_.begin();
    ConsensusChunk& merged = *iter;

#if 0
    if (chunks_.size() == 1) {
      iter->AlignLeftEndAgainstRead();
      for (int i = 0; i < iter->approx_read_offset; ++i)
        for (int j = 0; j < raw_read_[i].len; ++j)
          std::cerr << ' ';
      for (int i = 0; i < iter->consensus.size(); ++i)
        std::cerr << iter->consensus[i].str();
      std::cerr << std::endl;
    }
#endif

    ++iter;
    while (iter != chunks_.end()) {
      merged.TryMergeWith(*iter);
      iter = chunks_.erase(iter);
    }

    corrected_runs_ = std::move(merged.consensus);
  }
  
  std::string GetSequenceString() const {
    if (chunks_.empty())
      return "";
    std::string res;
    if (!corrected_runs_.empty()) {
      for (auto it = corrected_runs_.begin(); it != corrected_runs_.end(); ++it)
        res += it->str();
    } else {
      auto& runs = chunks_.front().consensus;
      for (auto it = runs.begin(); it != runs.end(); ++it)
        res += it->str();
    }
    return res;
  }
};


class SingleReadCorrector {
  const KMerData &kmer_data_;

 public:
  SingleReadCorrector(const KMerData &kmer_data) :
    kmer_data_(kmer_data) {}

  boost::optional<io::SingleRead> operator()(const io::SingleRead &r) {

#if 0
    std::cerr << "=============================================" << std::endl;

    std::cerr << '>' << r.name() << '\n'
              << r.GetSequenceString() << std::endl;
#endif

    CorrectedRead read(r, kmer_data_);
    read.MergeChunks();

    auto seq = read.GetSequenceString();
    if (seq.empty())
      return boost::optional<io::SingleRead>();

    return io::SingleRead(r.name(), seq);
  } 
};

}; // namespace correction
}; // namespace hammer
#endif // __HAMMER_IT_READ_CORRECTOR_HPP__
