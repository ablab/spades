//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "expander.hpp"

#include "config_struct_hammer.hpp"
#include "globals.hpp"
#include "kmer_data.hpp"
#include "valid_kmer_generator.hpp"

#include "io/reads/read.hpp"

#include <vector>
#include <cstring>

bool Expander::operator()(std::unique_ptr<Read> r) {
  uint8_t trim_quality = (uint8_t)cfg::get().input_trim_quality;

  // FIXME: Get rid of this
  Read cr = *r;
  size_t sz = cr.trimNsAndBadQuality(trim_quality);

  if (sz < hammer::K)
    return false;

  std::vector<unsigned> covered_by_solid(sz, false);
  std::vector<size_t> kmer_indices(sz, -1ull);

  ValidKMerGenerator<hammer::K> gen(cr);
  while (gen.HasMore()) {
    hammer::KMer kmer = gen.kmer();
    size_t idx = data_.checking_seq_idx(kmer);
    if (idx != -1ULL) {
      size_t read_pos = gen.pos() - 1;

      kmer_indices[read_pos] = idx;
      if (data_[idx].good()) {
        for (size_t j = read_pos; j < read_pos + hammer::K; ++j)
          covered_by_solid[j] = true;
      }
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
    if (!kmer_data.good()) {
#     pragma omp atomic
      changed_ += 1;

      kmer_data.lock();
      kmer_data.mark_good();
      kmer_data.unlock();
    }
  }
    
  return false;
}
