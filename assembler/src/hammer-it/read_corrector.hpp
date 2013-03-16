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
static int alignH(It1 a_begin, It1 a_initial_pos, It1 a_end, 
                  It2 b_initial_pos, It2 b_end,
                  size_t max_offset=3, size_t n_cmp=5, unsigned* p_score=NULL) {

  VERIFY(a_begin <= a_end);

  // compute most probable offset
  int best_offset = 0;
  unsigned best_score = -1U;

  int offset = 0;
  for (size_t i = 1; i <= max_offset * 2 + 1; ++i) {
    offset = (i / 2) * ((i & 1) ? 1 : -1); // 0, -1, 1, -2, 2, ...

    if (a_initial_pos + offset >= a_end)
      continue;

    if (a_initial_pos + offset < a_begin)
      continue;

    auto a_it = a_initial_pos + offset;
    auto b_it = b_initial_pos;

    unsigned score = 0;
    hammer::IonPairAligner<It1, It2> aligner(a_it, a_end, b_it, b_end);
    while (!aligner.empty()) {
      auto event = aligner.front();
      size_t x_offset = event.x_iter - a_it;
      size_t y_offset = event.y_iter - b_it;
      if (x_offset > n_cmp || y_offset > n_cmp)
        break;

      switch (event.type) {
        case kIonEventMismatch:
          score += 9 * event.length; break;
        case kIonEventBaseInsertion:
        case kIonEventBaseDeletion:
          score += 3 * event.length; break;
        case kIonEventRunInsertion:
        case kIonEventRunDeletion:
          score += 5 * event.length; break;
      }

      if (score > best_score)
        break;
      aligner.popFront();
    }

    if (score < best_score) {
      best_score = score;
      best_offset = offset;
    }

    if (best_score == 0)
      break;
  }

  if (p_score != NULL)
    *p_score = best_score;

  return best_offset;
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

class CorrectedRead {
  FlowSpaceRead raw_read_; // Uncorrected read
  const KMerData& kmer_data_;

  // Stores runs after joining chunks
  std::vector<hammer::HomopolymerRun> corrected_runs_;

  // Contiguous part of read with strong consensus
  struct ConsensusChunk {
    int approx_read_offset; // in the vector of raw read runs

    static const double kLowScoreThreshold = 1.0;

    enum {
      kChunkLeftAligned,
      kChunkRightAligned,
      kChunkNotAligned
    } alignment;

    const FlowSpaceRead& raw_read;

    std::vector<hammer::HomopolymerRun> consensus;
    std::vector<double> consensus_scores;

    ConsensusChunk(int approximate_read_offset,
                   const ScoreStorage &scores,
                   const FlowSpaceRead &read,
                   double trim_threshold=1.0)
        : approx_read_offset(approximate_read_offset), 
          alignment(kChunkNotAligned),
          raw_read(read)
    {
      bool left_trim = true;
      for (size_t i = 0; i < scores.size(); ++i) {
        auto run = hammer::iontorrent::consensus(scores[i]);

        // trim low-quality runs from the left side
        if (run.second <= kLowScoreThreshold && left_trim) {
          approx_read_offset += 1;
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

      consensus.resize(right_end);
      consensus_scores.resize(right_end);
    }

    unsigned AlignLeftEndAgainstRead() {
      unsigned score;
      const auto& data = raw_read.data();
      int offset = alignH(data.begin(), 
                          data.begin() + approx_read_offset,
                          data.end(),
                          consensus.begin(), consensus.end(), 5, 8, &score);
      approx_read_offset += offset;
      alignment = kChunkLeftAligned;
      return score;
    }

    unsigned AlignRightEndAgainstRead() {
      unsigned score;
      const auto& data = raw_read.data();

      //    approx_read_offset        approx_read_offset + consensus.size() - 1
      // ...--------------*----------------*-------------...                    <- read
      //                 ------------------  <- consensus
      int offset = alignH(data.rbegin(), 
                          data.rend() - (approx_read_offset + consensus.size() - 1) - 1,
                          data.rend(),
                          consensus.rbegin(), consensus.rend(), 5, 8, &score);
      approx_read_offset -= offset;
      alignment = kChunkRightAligned;
      return score;
    }

    bool TryMergeWith(ConsensusChunk& chunk) {
      if (chunk.consensus.empty())
        return true;
     
      if (alignment != kChunkRightAligned)
        AlignRightEndAgainstRead();
      if (chunk.alignment != kChunkLeftAligned)
        chunk.AlignLeftEndAgainstRead();

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
        size_t overlap = right_end_offset - chunk.approx_read_offset;
        if (overlap > chunk.consensus_scores.size())
          return false;

#if 0
        for (size_t i = 0; i < overlap; i++)
          std::cerr << *(consensus_scores.end() - i - 1) << ' ';
        std::cerr << " | ";
        for (size_t i = 0; i < overlap ; i++) 
          std::cerr << chunk.consensus_scores[i] << ' ';
        std::cerr << std::endl;
#endif

        consensus.insert(consensus.end(),
                         chunk.consensus.begin() + overlap,
                         chunk.consensus.end());

        consensus_scores.insert(consensus_scores.end(),
                                chunk.consensus_scores.begin() + overlap,
                                chunk.consensus_scores.end());
      }

      // right end needs to be aligned again
      alignment = kChunkNotAligned;
      return true;
    }
  };

  // Chunks where strong consensus was obtained
  std::list<ConsensusChunk> chunks_;

  void PushChunk(const ScoreStorage &scores, int approx_read_offset) {
    chunks_.push_back(ConsensusChunk(approx_read_offset, scores, raw_read_));
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

    unsigned approx_n_insertions = 0; // approximate number of insertions

    while (gen.HasMore()) {
      hammer::HKMer seq = gen.kmer();
      gen.Next();

      hammer::KMerStat k = kmer_data_[kmer_data_[seq].changeto];
      hammer::HKMer center = k.kmer;

      // if too many centers are skipped, start new chunk
      if (skipped > 4 && !start_new_chunk) {
        PushChunk(scores, approx_read_offset);
        start_new_chunk = true;
        pos += skipped;
        last_good_center_is_defined = false;
      } else {
        // otherwise, compare current center with the last 'good' center.
        // if they differ too much or the center has low quality, skip it.
        const double LOW_QUALITY_THRESHOLD = 1e-11;

        bool low_qual = k.qual > LOW_QUALITY_THRESHOLD;

        if (last_good_center_is_defined && !low_qual) {
          auto dist = distanceHKMer<2, 0, 1, 0, 1>
                                   (last_good_center.begin() + skipped + 1,
                                    last_good_center.end(),
                                    center.begin(),
                                    center.end() - 1 - skipped, 1);

          if (dist > 0) {
            low_qual = true;
#if 0
            std::cerr << "Skipped " << center 
                      << " (distance = " << dist << ")" << std::endl;
#endif
          }
        }

        if (low_qual) {
          need_to_align = true;
          skipped += 1;
          goto process_next_kmer;
        }
      }

      if (need_to_align && last_good_center_is_defined) {
        VERIFY(skipped + 1 < hammer::K);
        unsigned score;
        int offset = alignH(last_good_center.begin(), 
                            last_good_center.begin() + skipped + 1, 
                            last_good_center.end(),
                            center.begin(), center.end(), 3, 5, &score);
        if (score != 0) {
          if (!start_new_chunk)
            PushChunk(scores, approx_read_offset);
          start_new_chunk = true;
          pos += skipped;
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

    if (chunks_.size() == 1) {
      iter->AlignLeftEndAgainstRead();
#if 0
      for (int i = 0; i < iter->approx_read_offset; ++i)
        for (int j = 0; j < raw_read_[i].len; ++j)
          std::cerr << ' ';
      for (int i = 0; i < iter->consensus.size(); ++i)
        std::cerr << iter->consensus[i].str();
      std::cerr << std::endl;
#endif
    }

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
