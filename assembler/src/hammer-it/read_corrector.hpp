#ifndef __HAMMER_IT_READ_CORRECTOR_HPP__
#define __HAMMER_IT_READ_CORRECTOR_HPP__

#include "HSeq.hpp"
#include "hkmer_distance.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <deque>
#include <vector>
#include <algorithm>

#ifdef DEBUG_ION_CONSENSUS
#include <iomanip>
#include <iostream>
#endif

namespace hammer {
namespace correction {

namespace numeric = boost::numeric::ublas;

typedef numeric::matrix<double> ScoreMatrix;
typedef std::vector<ScoreMatrix> ScoreStorage;

struct EndsTrimmer {
  EndsTrimmer(unsigned trim_left, unsigned trim_right)
    : trim_left_(trim_left), trim_right_(trim_right) {}

  ScoreStorage::const_iterator trimmedBegin(const ScoreStorage& scores) const {
    return scores.begin() + trim_left_;
  }

  ScoreStorage::const_iterator trimmedEnd(const ScoreStorage& scores) const {
    return scores.end() - trim_right_;
  }

 private:
  unsigned trim_left_;
  unsigned trim_right_;
};

struct ClipTrimmedEnds {
  ClipTrimmedEnds(const io::SingleRead &r) {}

  std::string buildStart(size_t n) const { 
    return ""; 
  }

  std::string buildEnd(size_t n) const {
    return "";
  }
};

struct KeepTrimmedEnds {
  KeepTrimmedEnds(const io::SingleRead &r) : 
    runs_(hammer::iontorrent::toHomopolymerRuns(r.GetSequenceString())) {
  }

  // n - number of nucs trimmed from the left end
  std::string buildStart(size_t n) const {
    std::string res;
    for (size_t i = 0; i < n; ++i)
      res += runs_[i].str();
    return res;
  }

  // start_pos - offset of the start in the array of runs
  std::string buildEnd(size_t start_pos) const {
    size_t n = start_pos;
    if (runs_.size() < n)
      return "";
    std::string res;
    for (auto it = runs_.cbegin() + n; it != runs_.cend(); ++it)
      res += it->str();
    return res;
  }

 private:
  std::vector<hammer::HomopolymerRun> runs_;
};

// finds such offset D that sequences a[a_offset + D ..] and b[b_offset .. ]
// have longest common prefix
//
// parameters:
//    a, b, a_offset, b_offset - sequences and offsets on them
//    max_offset - maximum absolute value of offset that will be checked
//    n_cmp - maximum number of compared elements
//    p_score - address at which length of longest common prefix will be stored
//              (this length can't be more than n_cmp)
template <typename T1, typename T2>
static int alignH(const T1 &a, size_t a_offset, const T2 &b, size_t b_offset,
                  size_t max_offset=3, size_t n_cmp=5, int* p_score=NULL) {

  int a_len = hammer::internal::getSize(a);
  int b_len = hammer::internal::getSize(b);

  // compute most probable offset
  int best_offset = 0;
  int best_score = -1;

  int offset = 0;
  for (size_t i = 1; i <= max_offset * 2 + 1; ++i) {
    offset = (i / 2) * ((i & 1) ? 1 : -1); // 0, -1, 1, -2, 2, ...
    int score = 0;
    for (size_t j = 0; j < n_cmp; ++j) {
      int a_index = a_offset + offset + j;
      int b_index = b_offset + j;
      if (a_index < 0 || a_index >= a_len || b_index >= b_len)
        break;
      if (a[a_index].raw != b[b_index].raw) {
        // allow first run to differ
        if (j > 0)
          break;
      } else {
        ++score;
      }
    }

    if (score > best_score) {
      best_score = score;
      best_offset = offset;
    }
  }

  if (p_score != NULL)
    *p_score = best_score;

  return best_offset;
}

template <typename ReadTrimmer,
          typename TrimPolicy = KeepTrimmedEnds>
class SingleReadCorrector {
  const KMerData &kmer_data_;
  const ReadTrimmer &trimmer_;

#ifdef DEBUG_ION_CONSENSUS
  const std::string READ_NAME;
  const int START_POS;
  const int END_POS;
  std::string NUCS;
#endif

  struct Right {
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

  struct Left {
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

  template <typename Side>
  std::deque<hammer::HomopolymerRun> prolong(const hammer::HKMer &seed, 
                                             size_t bases_to_recover,
                                             Side) {
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
  SingleReadCorrector(const KMerData &kmer_data, const ReadTrimmer &trimmer) :
    kmer_data_(kmer_data), trimmer_(trimmer) {}

#ifdef DEBUG_ION_CONSENSUS
  SingleReadCorrector(const KMerData &kmer_data, const ReadTrimmer &trimmer,
                      std::string rname, int start_pos) :
    kmer_data_(kmer_data), trimmer_(trimmer), 
    READ_NAME(rname), START_POS(start_pos), END_POS(start_pos + 20), NUCS("ACGT") {}
#endif

  boost::optional<io::SingleRead> operator()(const io::SingleRead &r) {

    ScoreStorage scores(r.size() * 3 / 2, ScoreMatrix(4, 64, 0));

    std::string read_seq = r.GetSequenceString();
    auto original_runs = hammer::iontorrent::toHomopolymerRuns(read_seq);

    ValidHKMerGenerator<hammer::K> gen(r);
    int pos = -1; // offset on corrected read of 1st base of current center
    TrimPolicy trim_policy(r);

    int n_beg = gen.trimmed_left(); // number of runs skipped at the start
    unsigned skipped = 0; // low-quality centers
    bool need_to_align = false; // this is done after low-quality hkmers
    int n_insertions = 0; // APPROXIMATE number of insertions on the read
    hammer::HKMer last_good_center; // last high-quality center

    while (gen.HasMore()) {
      hammer::HKMer seq = gen.kmer();
      hammer::KMerStat k = kmer_data_[kmer_data_[seq].changeto];
      hammer::HKMer center = k.kmer;

      static const double QUALITY_THRESHOLD = 1e-11;

      // TODO: better detection of low-quality centers,
      //       based on _relative_ difference in quality
      if (k.qual < QUALITY_THRESHOLD) {

        if (need_to_align) {
          if (pos >= 0) {

            int offset; // of current good center relative to the previous one

            if (skipped > hammer::K - 4) { 
              // case of possible 'hole' -  try to recover a few extra bases
              auto to_recover = skipped - hammer::K + 6;

              auto left_runs = prolong(last_good_center, to_recover, Right());
              auto recovered_left = left_runs.size() - hammer::K;
              to_recover -= recovered_left;

              auto right_runs = prolong(center, 
                                        std::min(skipped - 3, to_recover), 
                                        Left());
              auto recovered_right = right_runs.size() - hammer::K;

              int score;
              offset = alignH(left_runs, skipped + 1 - recovered_right,
                              right_runs, 0, 3, 5, &score);

              if (score < 3) { // failed to get significant intersection
                int score, d, d1;
                // try to find right end of previous center on the read
                // (looking for last bases of left_runs doesn't make
                // much sense because the read is likely to have errors there)
                size_t prev_offset = pos + n_beg + n_insertions;
                d = alignH(original_runs, prev_offset + hammer::K - 5,
                           last_good_center, hammer::K - 5, 2, 5, &score);
                if (score < 4) {
                  pos += recovered_left;
                  break;
                }
                prev_offset += d;

                // try to find left end of current center on the read
                size_t curr_offset = prev_offset + skipped + 1;
                d1 = alignH(original_runs, curr_offset, center, 0, 2, 5, 
                            &score);
                if (score < 4) {
                  pos += recovered_left;
                  break;
                }
                curr_offset += d1;

                size_t read_offset_l = prev_offset + hammer::K;
                size_t read_offset_r = curr_offset;
                const double kCopiedScore = 1;

                pos += hammer::K;
                for (size_t j = read_offset_l; j < read_offset_r; ++j) {
                  auto run = original_runs[j];
                  scores[pos++](run.nucl, run.len) = kCopiedScore;
                }
                offset = skipped + 1 - curr_offset + prev_offset;
                pos += offset;

#ifdef DEBUG_ION_CONSENSUS
                std::cerr << r.name() << std::endl;
                for (size_t i = 0; i < 5; i++)
                  std::cerr << last_good_center[hammer::K - 5 + i].str() << " ";
                std::cerr << " | ";
                for (size_t j = read_offset_l; j < read_offset_r; ++j) {
                  auto run = original_runs[j];
                  std::cerr << run.str() << " ";
                }
                std::cerr << "| ";
                for (size_t i = 0; i < 5; i++)
                  std::cerr << center[i].str() << " ";
                std::cerr << std::endl << "offset: " << offset << std::endl;
                std::cerr << "skipped: " << skipped << std::endl;
#endif
              } else { // the two centers have in intersection >= 3 runs
                pos += skipped + offset + 1;
              }
            } else {
              offset = alignH(last_good_center, skipped + 1, center, 0);
              pos += skipped + offset + 1;
            }

            if (offset < 0) {
              n_insertions -= offset;
            }
          } else {
            pos = 0;
          }

          need_to_align = false;
          skipped = 0;
        } else {
          ++pos;
        }
       
#ifdef DEBUG_ION_CONSENSUS
        if (r.name() == READ_NAME && pos >= START_POS && pos < END_POS) {
          std::cerr << std::setprecision(5) << std::setw(8) << std::log(k.qual); 
          for (int i = 0; i < pos - START_POS; ++i)
            std::cerr << "    ";
          for (size_t i = 0; i < hammer::K; ++i)
            std::cerr << std::setw(3) << static_cast<unsigned>(center[i].len)
                      << NUCS[center[i].nucl];
          std::cerr << std::endl;
        }
#endif

        for (size_t i = 0; i < hammer::K; ++i)
          scores[pos + i](center[i].nucl, center[i].len) += k.count * (1 - k.qual);

        last_good_center = center;
      } else {

        need_to_align = true; 
        ++skipped;
        if (pos == -1)
          ++n_beg;
      }

      gen.Next();
    }

    if (pos == -1)
      return boost::optional<io::SingleRead>();

    scores.resize(pos + hammer::K);

    auto it_begin = trimmer_.trimmedBegin(scores);
    auto it_end = trimmer_.trimmedEnd(scores);

    if (it_begin >= it_end)
      return boost::optional<io::SingleRead>();

    // initialize sequence with bases trimmed from the left
    n_beg += it_begin - scores.begin();
    std::string res = trim_policy.buildStart(n_beg);

    std::array<hammer::HomopolymerRun, 5> last_corrected_runs;
    const int n_tail = last_corrected_runs.size();

    // take consensus for the middle of the read
    for (auto it = it_begin; it != it_end; ++it) {
      hammer::HomopolymerRun run = hammer::iontorrent::consensus(*it);
      res += run.str();

      if (it_end - it <= n_tail) {
        last_corrected_runs[n_tail - (it_end - it)] = run;
      }
    }

    // try to find trimmed bases
    if (it_end - it_begin >= n_tail && pos + n_beg + n_insertions >= n_tail) {
      int next_base_pos = pos + n_beg + n_insertions;
      int score;
      next_base_pos -= alignH(original_runs, pos + n_beg + n_insertions - n_tail,
                              last_corrected_runs, 0,
                              1, 5, &score);

      if (score == 5) {
        res += trim_policy.buildEnd(next_base_pos);
      }
    }

    return io::SingleRead(r.name(), res);
  }
};

}; // namespace correction
}; // namespace hammer
#endif // __HAMMER_IT_READ_CORRECTOR_HPP__
