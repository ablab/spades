#ifndef __HAMMER_IT_READ_CORRECTOR_HPP__
#define __HAMMER_IT_READ_CORRECTOR_HPP__

#include "HSeq.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>

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

template <typename ReadTrimmer>
class SingleReadCorrector {
  const KMerData &kmer_data_;
  const ReadTrimmer &trimmer_;

 public:
  SingleReadCorrector(const KMerData &kmer_data, const ReadTrimmer &trimmer) :
      kmer_data_(kmer_data), trimmer_(trimmer) {}

  boost::optional<io::SingleRead> operator()(const io::SingleRead &r) {

    ScoreStorage scores(r.size(), ScoreMatrix(4, 64, 0));

    ValidHKMerGenerator<hammer::K> gen(r);
    size_t pos = 0;
    while (gen.HasMore()) {
      hammer::HKMer seq = gen.kmer();
      hammer::KMerStat k = kmer_data_[kmer_data_[seq].changeto];
      hammer::HKMer center = k.kmer;

      for (size_t i = 0; i < hammer::K; ++i)
        scores[pos + i](center[i].nucl, center[i].len) += k.count * (1 - k.qual);

      gen.Next();
      pos += 1;
    }

    if (pos == 0)
      return boost::optional<io::SingleRead>();

    scores.resize(pos + hammer::K);

    std::string res;
    auto it_begin = trimmer_.trimmedBegin(scores);
    auto it_end = trimmer_.trimmedEnd(scores);

    if (it_begin >= it_end)
      return boost::optional<io::SingleRead>();

    for (auto it = it_begin; it < it_end; ++it) {
      hammer::HomopolymerRun run = hammer::iontorrent::consensus(*it);
      res += run.str();
    }

    return io::SingleRead(r.name(), res);
  }
};

}; // namespace correction
}; // namespace hammer
#endif // __HAMMER_IT_READ_CORRECTOR_HPP__
