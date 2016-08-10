//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "hamcluster.hpp"

#include "common/adt/concurrent_dsu.hpp"
#include "io/kmers/mmapped_reader.hpp"
#include "parallel_radix_sort.hpp"

#include "config_struct_hammer.hpp"
#include "globals.hpp"

#include <iostream>
#include <sstream>

class EncoderKMer {
public:
  inline static size_t extract(const SubKMer &x, unsigned shift, unsigned Base) {
    size_t idx = shift / SubKMer::TBits;
    size_t ishift = shift - idx * SubKMer::TBits;
    return (x.data()[idx] >> ishift) & ((1 << Base) - 1);
  }
};

struct SubKMerComparator {
    bool operator()(const SubKMer &l, const SubKMer &r) const {
      for (size_t i = 0; i < SubKMer::DataSize ; ++i) {
        if (l.data()[i] != r.data()[i]) {
          return (l.data()[i] < r.data()[i]);
        }
      }

      return false;
    }
};

template<class Op>
std::pair<size_t, size_t> SubKMerSplitter::split(Op &&op) {
  std::vector<SubKMer> data; std::vector<size_t> blocks;

  MMappedReader bifs(bifname_, /* unlink */ true);
  MMappedReader kifs(kifname_, /* unlink */ true);
  size_t icnt = 0, ocnt = 0;
  while (bifs.good()) {
    deserialize(blocks, data, bifs, kifs);

    using PairSort = parallel_radix_sort::PairSort<SubKMer, size_t, SubKMer, EncoderKMer>;
    // PairSort::InitAndSort(data.data(), blocks.data(), data.size());
    PairSort::InitAndSort(data.data(), blocks.data(), data.size(), data.size() > 1000*16 ? -1 : 1);

    for (auto start = data.begin(), end = data.end(); start != end;) {
      auto chunk_end = std::upper_bound(start + 1, data.end(), *start, SubKMerComparator());
      op(blocks.begin() + (start - data.begin()), chunk_end - start);
      start = chunk_end;
      ocnt += 1;
    }
    icnt += 1;
  }

  return std::make_pair(icnt, ocnt);
}

#if 1
static bool canMerge(const ConcurrentDSU &uf, size_t x, size_t y) {
  size_t szx = uf.set_size(x), szy = uf.set_size(y);
  const size_t hardthr = 2500;

  // Global threshold - no cluster larger than hard threshold
  if (szx + szy > hardthr)
    return false;

  // If one of the clusters is moderately large, than attach "almost" singletons
  // only.
  if ((szx > hardthr * 3 / 4 && szy > 50) ||
      (szy > hardthr * 3 / 4 && szx > 50))
    return false;

  return true;
}
#else
static bool canMerge(const ConcurrentDSU &uf, size_t x, size_t y) {
  return (uf.set_size(x) + uf.set_size(y)) < 10000;
}
#endif


static void processBlockQuadratic(ConcurrentDSU  &uf,
                                  const std::vector<size_t>::iterator &block,
                                  size_t block_size,
                                  const KMerData &data,
                                  unsigned tau) {
  for (size_t i = 0; i < block_size; ++i) {
    size_t x = block[i];
    hammer::KMer kmerx = data.kmer(x);
    for (size_t j = i + 1; j < block_size; j++) {
      size_t y = block[j];
      hammer::KMer kmery = data.kmer(y);
      if (!uf.same(x, y) &&
          canMerge(uf, x, y) &&
          hamdistKMer(kmerx, kmery, tau) <= tau) {
        uf.unite(x, y);
      }
    }
  }
}

void KMerHamClusterer::cluster(const std::string &prefix,
                               const KMerData &data,
                               ConcurrentDSU &uf) {
  // First pass - split & sort the k-mers
  std::string fname = prefix + ".first", bfname = fname + ".blocks", kfname = fname + ".kmers";
  std::ofstream bfs(bfname, std::ios::out | std::ios::binary);
  std::ofstream kfs(kfname, std::ios::out | std::ios::binary);
  VERIFY(bfs.good()); VERIFY(kfs.good());

  INFO("Serializing sub-kmers.");
  for (unsigned i = 0; i < tau_ + 1; ++i) {
    size_t from = (*Globals::subKMerPositions)[i];
    size_t to = (*Globals::subKMerPositions)[i+1];

    INFO("Serializing: [" << from << ", " << to << ")");
    serialize(bfs, kfs,
              data, NULL, 0,
              SubKMerPartSerializer(from, to));
  }
  VERIFY(!bfs.fail()); VERIFY(!kfs.fail());
  bfs.close(); kfs.close();

  size_t big_blocks1 = 0;
  {
    unsigned block_thr = cfg::get().hamming_blocksize_quadratic_threshold;

    INFO("Splitting sub-kmers, pass 1.");
    SubKMerSplitter Splitter(bfname, kfname);

    fname = prefix + ".second", bfname = fname + ".blocks", kfname = fname + ".kmers";
    bfs.open(bfname, std::ios::out | std::ios::binary);
    kfs.open(kfname, std::ios::out | std::ios::binary);
    VERIFY(bfs.good()); VERIFY(kfs.good());

    std::pair<size_t, size_t> stat =
      Splitter.split([&] (const std::vector<size_t>::iterator &start, size_t sz) {
        if (sz < block_thr) {
          // Merge small blocks.
          processBlockQuadratic(uf, start, sz, data, tau_);
        } else {
          big_blocks1 += 1;
          // Otherwise - dump for next iteration.
          for (unsigned i = 0; i < tau_ + 1; ++i) {
            serialize(bfs, kfs,
                      data, &start, sz,
                      SubKMerStridedSerializer(i, tau_ + 1));
          }
        }
    });
    INFO("Splitting done."
         " Processed " << stat.first << " blocks."
         " Produced " << stat.second << " blocks.");

    // Sanity check - there cannot be more blocks than tau + 1 times of total
    // kmer number. And on the first pass we have only tau + 1 input blocks!
    VERIFY(stat.first == tau_ + 1);
    VERIFY(stat.second <= (tau_ + 1) * data.size());

    VERIFY(!bfs.fail()); VERIFY(!kfs.fail());
    bfs.close(); kfs.close();
    INFO("Merge done, total " << big_blocks1 << " new blocks generated.");
  }

  size_t big_blocks2 = 0;
  {
    INFO("Spliting sub-kmers, pass 2.");
    SubKMerSplitter Splitter(bfname, kfname);
    size_t nblocks = 0;
    std::pair<size_t, size_t> stat =
      Splitter.split([&] (const std::vector<size_t>::iterator &start, size_t sz) {
        if (sz > 50) {
          big_blocks2 += 1;
#if 0
          for (size_t i = 0; i < block.size(); ++i) {
            std::string s(Globals::blob + data[block[i]], K);
            INFO("" << block[i] << ": " << s);
          }
#endif
        }
        processBlockQuadratic(uf, start, sz, data, tau_);
        nblocks += 1;
    });
    INFO("Splitting done."
            " Processed " << stat.first << " blocks."
            " Produced " << stat.second << " blocks.");

    // Sanity check - there cannot be more blocks than tau + 1 times of total
    // kmer number. And there should be tau + 1 times big_blocks input blocks.
    VERIFY(stat.first == (tau_ + 1)*big_blocks1);
    VERIFY(stat.second <= (tau_ + 1) * (tau_ + 1) * data.size());

    INFO("Merge done, saw " << big_blocks2 << " big blocks out of " << nblocks << " processed.");
  }
}

enum {
  UNLOCKED = 0,
  PARTIALLY_LOCKED = 1,
  FULLY_LOCKED = 3
};

static bool canMerge2(const ConcurrentDSU &uf, size_t kidx, size_t cidx) {
    // If either of indices is fully locked - bail out
    uint64_t kaux = uf.root_aux(kidx), caux = uf.root_aux(cidx);
    if (kaux == FULLY_LOCKED || caux == FULLY_LOCKED)
        return false;

    // Otherwise there is a possibility to merge stuff.
    if (0 && (kaux == PARTIALLY_LOCKED || caux == PARTIALLY_LOCKED)) {
        // We cannot merge two partially locked clusters.
        return kaux != caux;
    }

    return true;
}

static void ClusterChunk(size_t start_idx, size_t end_idx, const KMerData &data, ConcurrentDSU &uf) {
    unsigned nthreads = cfg::get().general_max_nthreads;

    // INFO("Cluster: " << start_idx << ":" << end_idx);
#   pragma omp parallel num_threads(nthreads)
    {
#       pragma omp for
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
            hammer::KMer kmer = data.kmer(idx);

            if (kmer.GetHash() > (!kmer).GetHash())
                continue;

            size_t kidx = data.seq_idx(kmer);
            size_t rckidx = -1ULL;
            // INFO("" << kmer << ":" << kidx);

            for (size_t k = 0; k < hammer::K; ++k) {
                hammer::KMer candidate = kmer;
                char c = candidate[k];
                for (char nc = 0; nc < 4; ++nc) {
                    if (nc == c)
                        continue;
                    candidate.set(k, nc);
                    size_t cidx = data.checking_seq_idx(candidate);
                    // INFO("" << candidate << ":" << cidx);
                    if (cidx != -1ULL && canMerge2(uf, kidx, cidx)) {
                        uf.unite(kidx, cidx);

                        size_t rccidx = data.seq_idx(!candidate);
                        if (rckidx == -1ULL)
                            rckidx = data.seq_idx(!kmer);
                        uf.unite(rckidx, rccidx);
                    }
                }
            }
        }
#       pragma omp barrier
        //INFO("Lock: " << start_idx << ":" << end_idx);
#       pragma omp for
        for (size_t idx = start_idx; idx < end_idx; ++idx) {
            if (uf.set_size(idx) < 2500)
                continue;

            if (uf.root_aux(idx) != FULLY_LOCKED)
                uf.set_root_aux(idx, FULLY_LOCKED);
        }
    }
}

void TauOneKMerHamClusterer::cluster(const std::string &, const KMerData &data, ConcurrentDSU &uf) {
    size_t start_idx = 0;
    while (start_idx < data.size()) {
        size_t end_idx = start_idx + 64*1024;
        if (end_idx > data.size())
            end_idx = data.size();

        ClusterChunk(start_idx, end_idx, data, uf);

        start_idx = end_idx;
    }
}
