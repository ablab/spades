#include "read_corrector.hpp"

#include "kmer_data.hpp"
#include "kmer_stat.hpp"
#include "valid_kmer_generator.hpp"

#include <string>
#include <vector>

using namespace hammer;

static bool update(const std::string &seq, const KMerData &data,
                   size_t pos, const KMerStat & stat,
                   std::vector<std::vector<unsigned> > & v,
                   int & left, int & right, bool & isGood,
                   bool correct_threshold, bool discard_singletons) {
  bool res = false;
  if (stat.isGoodForIterative() || (correct_threshold && stat.isGood())) {
    isGood = true;
    KMer kmer = stat.kmer();

    for (size_t j = 0; j < K; ++j)
      v[kmer[j]][pos + j]++;
    if ((int) pos < left)
      left = pos;
    if ((int) pos > right)
      right = pos;
  } else {
    // if discard_only_singletons = true, we always use centers of clusters that do not coincide with the current center
    if (stat.change() &&
        (discard_singletons || data[stat.changeto].isGoodForIterative() ||
         (correct_threshold && stat.isGood()))) {

      isGood = true;
      if ((int) pos < left)
        left = pos;
      if ((int) pos > right)
        right = pos;
      KMer newkmer = data[stat.changeto].kmer();

      for (size_t j = 0; j < K; ++j) {
        v[newkmer[j]][pos + j]++;
      }

      res = true;
    }
  }

  return res;
}

bool ReadCorrector::CorrectOneRead(Read & r,
                                   bool correct_threshold, bool discard_singletons, bool discard_bad) {
  bool isGood = false;
  std::string seq = r.getSequenceString();
  const std::string &qual = r.getQualityString();
  
  unsigned read_size = seq.size();

  VERIFY(read_size >= K);

  // create auxiliary structures for consensus
  std::vector<std::vector<unsigned> > v(4, std::vector<unsigned>(read_size, 0));  // A=0, C=1, G=2, T=3
  isGood = false;

  // getting the leftmost and rightmost positions of a solid kmer
  int left = read_size; int right = -1;
  bool changedRead = false;
  ValidKMerGenerator<K> gen(seq.data(), qual.data(), read_size);
  while (gen.HasMore()) {
    size_t read_pos = gen.pos() - 1;
    const KMerStat &kmer_data = data_[gen.kmer()];
 
    changedRead = changedRead ||
                  update(seq, data_, read_pos, kmer_data, v,
                         left, right, isGood, correct_threshold, discard_singletons);

    gen.Next();
  }

  int left_rev = 0; int right_rev = read_size-(int)K;

  if (left <= right && left_rev <= right_rev) {
    left = std::min(left, (int)read_size - left_rev - (int)K);
    right = std::max(right, (int)read_size - right_rev - (int)K);
  } else if ( left > right && left_rev <= right_rev ) {
    left = (int)read_size - left_rev - (int)K;
    right = (int)read_size - right_rev - (int)K;
  }

  // at this point the array v contains votes for consensus

  size_t res = 0; // how many nucleotides have really changed?
  // find max consensus element
  for (size_t j=0; j < read_size; ++j) {
    char cmax = seq[j]; unsigned nummax = 0;
    for (size_t k=0; k<4; ++k) {
      if (v[k][j] > nummax) {
        cmax = nucl(k); nummax = v[k][j];
      }
    }
    if (seq[j] != cmax)
      ++res;
    seq[j] = cmax;
  }

  r.setSequence(seq.data(), /* preserve_trimming */ true);

  // if discard_bad=false, we retain original sequences when needed
  if (discard_bad) {
    r.trimLeftRight(left, right+K-1);
    if (left > 0 || right + K -1 < read_size)
      changedRead = true;
  }

  if (res) {
#   pragma omp atomic
    changed_nucleotides_ += res;

#   pragma omp atomic
    changed_reads_ += 1;
  }

  return isGood;
}
