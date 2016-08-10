//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "subcluster.hpp"
#include "config_struct.hpp"
#include "consensus.hpp"
#include "hkmer_distance.hpp"
#include "kmer_data.hpp"
#include "utils/logger/log_writers.hpp"

#include <boost/numeric/ublas/matrix.hpp>

#include <vector>
#include <iostream>

hammer::HKMer center(const KMerData &data, const std::vector<size_t>& kmers) {
  hammer::HKMer res;
  namespace numeric = boost::numeric::ublas;

  for (unsigned i = 0; i < hammer::K; ++i) {
    numeric::matrix<double> scores(4, 64, 0);
    for (size_t j = 0; j < kmers.size(); ++j) {
      const hammer::KMerStat &k = data[kmers[j]];
      // FIXME: switch to MLE when we'll have use per-run quality values
#if 1
      scores(k.kmer[i].nucl, k.kmer[i].len) += double(k.count) * (1 - k.qual);
#else
      for (unsigned n = 0; n < 4; ++n)
        for (unsigned l = 1; l < 64; ++l)
          scores(n, l) += k.count * (n == k.kmer[i].nucl && l == k.kmer[i].len ?
                                     log(1 - k.qual) : log(k.qual) - log(4*63 - 1));
#endif
    }

    res[i] = hammer::iontorrent::consensus(scores).first;
  }

  return res;
}

bool assign(KMerData &kmer_data, const std::vector<size_t> &cluster) {
  hammer::HKMer c = center(kmer_data, cluster);
  bool nonread = false;

  size_t idx = kmer_data.seq_idx(c);
  if (kmer_data[idx].kmer != c) {
#   pragma omp critical
    {
      idx = kmer_data.push_back(hammer::KMerStat(0, c, 1.0));
    }
    nonread = true;
  }

  for (size_t j = 0; j < cluster.size(); ++j)
    kmer_data[cluster[j]].changeto = unsigned(idx);

  return nonread;
}

void dump(const KMerData &kmer_data, const std::vector<size_t> &cluster) {
  std::cerr << "{ \n\"kmers\": {";
  for (size_t j = 0; j < cluster.size(); ++j) {
    if (j > 0) std::cerr << ", ";
    std::cerr << '"' << kmer_data[cluster[j]].kmer << "\": ["
              << kmer_data[cluster[j]].count << ", " 
              << 1 - kmer_data[cluster[j]].qual << "] \n";
  }
  std::cerr << "}, \"center\": { \"status\": ";
  hammer::HKMer c = center(kmer_data, cluster);
  size_t idx = kmer_data.seq_idx(c);
  if (kmer_data[idx].kmer == c) {
    std::cerr << "\"ok\", \"center\": \"" << c << "\"}\n";
  } else {
    std::cerr << "\"not\", \"kmer\": \"" << kmer_data[idx].kmer 
              << "\", \"center\": \"" << c << "\"}\n";
  }
  std::cerr << "}" << std::endl;
}

size_t subcluster(KMerData &kmer_data, std::vector<size_t> &cluster) {
  size_t nonread = 0;

  // First, sort the kmer indicies wrt count
  std::sort(cluster.begin(), cluster.end(), CountCmp(kmer_data));

  // The number of subclusters for now is really dumb: we assume that the quality should be 1.
  size_t k = 0;
  for (size_t i = 0; i < cluster.size(); ++i)
      k += kmer_data[cluster[i]].qual < cfg::get().center_qual_threshold;

  if (k <= 1) {
#if 0
    dump(kmer_data, cluster);
#endif
    return assign(kmer_data, cluster);
  }

  // Find the closest center
  std::vector<std::vector<size_t> > idx(k, std::vector<size_t>());
  for (size_t i = 0; i < k; ++i)
    idx[i].push_back(cluster[i]);
  for (size_t i = k; i < cluster.size(); ++i) {
    unsigned dist = std::numeric_limits<unsigned>::max();
    size_t cidx = k;
    hammer::HKMer kmerx = kmer_data[cluster[i]].kmer;
    for (size_t j = 0; j < k; ++j) {
      hammer::HKMer kmery = kmer_data[cluster[j]].kmer;
      unsigned cdist = hammer::distanceHKMer(kmerx.begin(), kmerx.end(),
                                             kmery.begin(), kmery.end());
      if (cdist < dist) {
        cidx = j;
        dist = cdist;
      }
    }
    VERIFY(cidx < k);
    idx[cidx].push_back(cluster[i]);
  }

  for (auto it = idx.begin(), et = idx.end(); it != et; ++it) {
    const std::vector<size_t> &subcluster = *it;

    if (assign(kmer_data, subcluster)) {
      nonread += 1;
#if 0
      dump(kmer_data, subcluster);
#endif
    }
  }

  return nonread;
}
