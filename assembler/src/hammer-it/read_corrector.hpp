#ifndef __HAMMER_IT_READ_CORRECTOR_HPP__
#define __HAMMER_IT_READ_CORRECTOR_HPP__

#include "HSeq.hpp"
#include "hkmer_distance.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>

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

    ValidHKMerGenerator<hammer::K> gen(r);
    int pos = -1;
    TrimPolicy trim_policy(r);

    unsigned skipped_hkmers = 0; // number of consecutive hkmers with low quality
    bool need_to_align = false; // this is done after low-quality hkmers
    int skipped_at_start = 0;
    int n_insertions = 0; // APPROXIMATE number of insertions on the read

    hammer::HKMer last_good_center; // last high-quality center

    while (gen.HasMore()) {
      hammer::HKMer seq = gen.kmer();
      hammer::KMerStat k = kmer_data_[kmer_data_[seq].changeto];
      hammer::HKMer center = k.kmer;

      static const double QUALITY_THRESHOLD = 1e-11;

      if (k.qual < QUALITY_THRESHOLD) {

        if (need_to_align) {
          if (pos >= 0) {
            int offset = alignH(last_good_center, skipped_hkmers + 1,
                                center, 0);

            if (offset < 0)
              n_insertions -= offset;
            pos += skipped_hkmers + offset + 1;
          } else {
            pos = 0;
          }

          need_to_align = false;
          skipped_hkmers = 0;
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
        if (skipped_hkmers >= hammer::K / 2) {
          size_t start = pos + skipped_hkmers + 1;
          for (size_t i = 0; i < hammer::K; ++i)
            scores[start + i](center[i].nucl, center[i].len) += k.count * (1 - k.qual);
        }

        need_to_align = true; 
        ++skipped_hkmers;
        if (pos == -1)
          ++skipped_at_start;
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

    size_t triml = gen.trimmed_left();

    // initialize sequence with bases trimmed from the left
    int n_beg = it_begin - scores.begin() + triml + skipped_at_start;
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
      std::string seq = r.GetSequenceString();
      auto original_runs = hammer::iontorrent::toHomopolymerRuns(seq);
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
