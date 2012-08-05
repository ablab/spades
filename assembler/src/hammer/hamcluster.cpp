//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "hamcluster.hpp"

#include "globals.hpp"
#include "mmapped_reader.hpp"
#include "union.hpp"

#include <iostream>
#include <sstream>

#ifdef USE_GLIBCXX_PARALLEL
#include <parallel/algorithm>
#endif

struct SubKMerComparator {
  bool operator()(const SubKMer &lhs, const SubKMer &rhs) {
    return Seq<K>::less2()(lhs.data, rhs.data);
  }
};

std::pair<size_t, size_t> SubKMerSplitter::split() {
  std::vector<SubKMer> data;

  MMappedReader ifs(ifname_, /* unlink */ true);
  std::ofstream ofs(ofname_, std::ios::out | std::ios::binary);
  VERIFY(ofs.good());
  size_t icnt = 0, ocnt = 0;
  while (ifs.good()) {
    SubKMerComparator comp;

    deserialize(data, ifs);

#ifdef USE_GLIBCXX_PARALLEL
    // Explicitly force a call to parallel sort routine.
    __gnu_parallel::sort(data.begin(), data.end(), comp);
#else
    std::sort(data.begin(), data.end(), comp);
#endif
    for (auto start = data.begin(), end = data.end(); start != end;) {
      auto chunk_end = std::upper_bound(start + 1, data.end(), *start, comp);
      serialize(ofs, start, chunk_end);
      start = chunk_end;
      ocnt += 1;
    }
    icnt += 1;
  }
  VERIFY(!ofs.fail());

  ofs.close();

  return make_pair(icnt, ocnt);
}

static unsigned hamdistKMer(hint_t x, hint_t y, unsigned tau) {
  unsigned dist = 0;
  for (unsigned i = 0; i < K; ++i) {
    if (Globals::blob[x + i] != Globals::blob[y + i]) {
      ++dist; if (dist > tau) return dist;
    }
  }
  return dist;
}

#if 0
static bool canMerge(const unionFindClass &uf, int x, int y) {
  size_t szx = uf.set_size(x), szy = uf.set_size(y);
  const size_t hardthr = 15000;

  // Global threshold - no cluster larger than hard threshold
  if (szx + szy > hardthr)
    return false;

  // If one of the clusters is moderately large, than attach "almost" singletons
  // only.
  if ((szx > hardthr * 3 / 4 && szy > 10) ||
      (szy > hardthr * 3 / 4 && szx > 10))
    return false;

  return true;
}
#else
static bool canMerge(const unionFindClass &uf, int x, int y) {
  return (uf.set_size(x) + uf.set_size(y)) < 10000;
}
#endif


static void processBlockQuadratic(unionFindClass &uf,
                                  const std::vector<size_t> &block,
                                  const KMerIndex &index,
                                  unsigned tau) {
  size_t blockSize = block.size();
  for (size_t i = 0; i < blockSize; ++i) {
    hint_t kmerx = index[block[i]].first.start();
    uf.find_set(block[i]);
    for (uint32_t j = i + 1; j < blockSize; j++) {
      hint_t kmery = index[block[j]].first.start();
      if (hamdistKMer(kmerx, kmery, tau) <= tau &&
          canMerge(uf, block[i], block[j])) {
        uf.unionn(block[i], block[j]);
      }
    }
  }
}

void KMerHamClusterer::cluster(const std::string &prefix,
                               const KMerIndex &index,
                               unionFindClass &uf) {
  // First pass - split & sort the k-mers
  std::ostringstream tmp;
  tmp << prefix << ".first";
  std::string fname(tmp.str());
  std::ofstream ofs(fname, std::ios::out | std::ios::binary);
  VERIFY(ofs.good());

  INFO("Serializing sub-kmers.");
  for (unsigned i = 0; i < tau_ + 1; ++i) {
    size_t from = (*Globals::subKMerPositions)[i];
    size_t to = (*Globals::subKMerPositions)[i+1];

    INFO("Serializing: [" << from << ", " << to << ")");
    serialize(ofs, index, NULL,
              SubKMerPartSerializer(from, to));
  }
  VERIFY(!ofs.fail());
  ofs.close();

  size_t big_blocks1 = 0;
  {
    INFO("Splitting sub-kmers, pass 1.");
    SubKMerSplitter Splitter(fname, fname + ".blocks");
    std::pair<size_t, size_t> stat = Splitter.split();
    INFO("Splitting done."
            " Processed " << stat.first << " blocks."
            " Produced " << stat.second << " blocks.");

    // Sanity check - there cannot be more blocks than tau + 1 times of total
    // kmer number. And on the first pass we have only tau + 1 input blocks!
    VERIFY(stat.first == tau_ + 1);
    VERIFY(stat.second <= (tau_ + 1) * index.size());

    // Ok, now in the files we have everything grouped in blocks in the output files.

    std::vector<size_t> block;

    INFO("Merge sub-kmers, pass 1");
    SubKMerBlockFile blocks(fname + ".blocks", /* unlink */ true);

    std::ostringstream tmp;
    tmp << prefix << ".second";
    fname = tmp.str();

    ofs.open(fname, std::ios::out | std::ios::binary);
    VERIFY(ofs.good());
    while (blocks.get_block(block)) {
      unsigned block_thr = cfg::get().subvectors_blocksize_quadratic_threshold;
      if (block.size() < block_thr) {
        // Merge small blocks.
        processBlockQuadratic(uf, block, index, tau_);
      } else {
        big_blocks1 += 1;
        // Otherwise - dump for next iteration.
        for (unsigned i = 0; i < tau_ + 1; ++i) {
          serialize(ofs, index, &block,
                    SubKMerStridedSerializer(i, tau_ + 1));
        }
      }
    }
    VERIFY(!ofs.fail());
    ofs.close();
    INFO("Merge done, total " << big_blocks1 << " new blocks generated.");
  }

  size_t big_blocks2 = 0;
  {
    INFO("Spliting sub-kmers, pass 2.");
    SubKMerSplitter Splitter(fname, fname + ".blocks");
    std::pair<size_t, size_t> stat = Splitter.split();
    INFO("Splitting done."
            " Processed " << stat.first << " blocks."
            " Produced " << stat.second << " blocks.");

    // Sanity check - there cannot be more blocks than tau + 1 times of total
    // kmer number. And there should be tau + 1 times big_blocks input blocks.
    VERIFY(stat.first == (tau_ + 1)*big_blocks1);
    VERIFY(stat.second <= (tau_ + 1) * (tau_ + 1) * index.size());

    INFO("Merge sub-kmers, pass 2");
    SubKMerBlockFile blocks(fname + ".blocks", /* unlink */ true);
    std::vector<size_t> block;

    size_t nblocks = 0;
    while (blocks.get_block(block)) {
      if (block.size() > 50) {
        big_blocks2 += 1;
#if 0
        for (size_t i = 0; i < block.size(); ++i) {
          std::string s(Globals::blob + index[block[i]], K);
          INFO("" << block[i] << ": " << s);
        }
#endif
      }
      processBlockQuadratic(uf, block, index, tau_);
      nblocks += 1;
    }
    INFO("Merge done, saw " << big_blocks2 << " big blocks out of " << nblocks << " processed.");
  }
}
