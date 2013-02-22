#include "expander.hpp"

#include "config_struct_hammer.hpp"
#include "globals.hpp"
#include "kmer_data.hpp"
#include "valid_kmer_generator.hpp"

#include "io/read.hpp"

#include <vector>
#include <cstring>

bool Expander::operator()(const Read &r) {
  int trim_quality = cfg::get().input_trim_quality;

  // FIXME: Get rid of this
  Read cr = r;
  size_t sz = cr.trimNsAndBadQuality(trim_quality);

  if (sz < hammer::K)
    return false;

  std::vector<unsigned> covered_by_solid(sz, false);
  std::vector<size_t> kmer_indices(sz, -1ull);
    
  ValidKMerGenerator<hammer::K> gen(cr);
  while (gen.HasMore()) {
    hammer::KMer kmer = gen.kmer();
    size_t idx = data_.seq_idx(kmer);
    size_t read_pos = gen.pos() - 1;

    kmer_indices[read_pos] = idx;
    if (data_[idx].isGoodForIterative()) {
      for (size_t j = read_pos; j < read_pos + hammer::K; ++j)
        covered_by_solid[j] = true;
    }
    gen.Next();
  }

  for (size_t j = 0; j < sz; ++j)
    if (!covered_by_solid[j])
      return false;

  for (size_t j = 0; j < sz; ++j) {
    if (kmer_indices[j] == -1ull)
      continue;

    // FIXME: Do not lock everything
    KMerStat &kmer_data = data_[kmer_indices[j]];
    if (!kmer_data.isGoodForIterative() &&
        !kmer_data.isMarkedGoodForIterative()) {
#     pragma omp atomic
      changed_ += 1;

      kmer_data.lock();
      kmer_data.makeGoodForIterative();
      kmer_data.unlock();
    }
  }
    
  return false;
}
