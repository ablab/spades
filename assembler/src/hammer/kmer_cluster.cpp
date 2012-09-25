//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * kmer_cluster.cpp
 *
 *  Created on: 16.07.2011
 *      Author: snikolenko
 */

#include "standard.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>

#include "read/ireadstream.hpp"
#include "mathfunctions.hpp"
#include "hammer_tools.hpp"
#include "hamcluster.hpp"
#include "kmer_cluster.hpp"
#include "config_struct_hammer.hpp"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using std::max_element;
using std::min_element;

namespace numeric = boost::numeric::ublas;

using namespace hammer;

std::string KMerClustering::GetGoodKMersFname() const {
  // FIXME: This is ugly!
  std::ostringstream tmp;
  tmp.str("");
  tmp << workdir_ << "/" << "kmers.solid";

  return tmp.str();
}

std::string KMerClustering::GetBadKMersFname() const {
  // FIXME: This is ugly!
  std::ostringstream tmp;
  tmp.str("");
  tmp << workdir_ << "/" << "kmers.bad";

  return tmp.str();
}

static double logLikelihoodKMer(const KMer center, const KMerStat &x) {
  double res = 0;
  KMer kmer = x.kmer();
  for (unsigned i = 0; i < K; ++i) {
    if (center[i] != kmer[i]) {
      res += getRevProb(x, i, /* log */ true) - log(3);
    } else {
      res += getProb(x, i, /* log */ true);
    }
  }

  return res;
}

KMer KMerClustering::ConsensusWithMask(const vector<unsigned> & cl, const vector<unsigned> & mask, unsigned maskVal) const {
  size_t blockSize = cl.size();

  // consensus of a single string is trivial
  if (blockSize == 1)
    return data_[cl[0]].kmer();

  uint64_t scores[4*K] = {0};
  for (size_t j = 0; j < blockSize; ++j) {
    if (mask[j] != maskVal)
      continue;

    const KMerStat &kms = data_[cl[j]];
    const KMer kmer = kms.kmer();

    for (unsigned i = 0; i < K; ++i) {
      scores[4*i + kmer[i]] += kms.count;
    }
  }

  std::string res(K, 'A');
  for (unsigned i = 0; i < K; ++i) {
    res[i] = nucl(std::max_element(scores + 4*i, scores + 4*i + 4) - (scores + 4*i));
  }

  return KMer(res);
}

KMer KMerClustering::Consensus(const vector<unsigned> & cl) const {
  size_t blockSize = cl.size();

  // consensus of a single string is trivial
  if (blockSize == 1)
    return data_[cl[0]].kmer();

  uint64_t scores[4*K] = {0};
  for (size_t j = 0; j < blockSize; ++j) {
    const KMerStat &kms = data_[cl[j]];
    const KMer kmer = kms.kmer();

    for (unsigned i = 0; i < K; ++i) {
      scores[4*i + kmer[i]] += kms.count;
    }
  }

  std::string res(K, 'A');
  for (unsigned i = 0; i < K; ++i) {
    res[i] = nucl(std::max_element(scores + 4*i, scores + 4*i + 4) - (scores + 4*i));
  }

  return KMer(res);
}

double KMerClustering::ClusterBIC(const vector<unsigned> & cl, const vector<StringCount> & centers, const vector<unsigned> & indices) const {
  size_t blockSize = cl.size();
  size_t clusters = centers.size();
  if (blockSize == 0)
    return -std::numeric_limits<double>::infinity();
  assert(blockSize == indices.size());
  assert(centers.size() > 0);

  double loglik = 0;
  unsigned total = 0;
  for (size_t i=0; i<blockSize; ++i) {
    const KMerStat &kms = data_[cl[i]];
    unsigned cnt = kms.count;
    loglik += logLikelihoodKMer(centers[indices[i]].first, kms);
    total += cnt;
  }

  unsigned nparams = (clusters - 1) + clusters*K + 2*clusters*K;

  return loglik - nparams*log(blockSize) / 2;
}

double KMerClustering::lMeansClustering(unsigned l, const std::vector<unsigned> &kmerinds,
                                        std::vector<unsigned> & indices, std::vector<StringCount> & centers) {
  centers.resize(l); // there are l centers

  // if l==1 then clustering is trivial
  if (l == 1) {
    centers[0].first = Consensus(kmerinds);
    centers[0].second.first = kmerinds.size();
    centers[0].second.second = 1;
    for (size_t i = 0; i < kmerinds.size(); ++i)
      indices[i] = 0;
    return ClusterBIC(kmerinds, centers, indices);
  }

  // We assume that kmerinds are sorted wrt the count.
  for (unsigned j = 0; j < l; ++j) {
    const KMerStat &kms = data_[kmerinds[j]];
    centers[j].first = kms.kmer();
    centers[j].second.first = kms.count;
    centers[j].second.second = kms.totalQual;
  }

  if (cfg::get().bayes_debug_output > 1) {
#   pragma omp critical
    {
      std::cout << "    centers:\n";
      for (size_t i=0; i < centers.size(); ++i) {
        std::cout << "    " << centers[i].first << "\n";
      }
    }
  }

  // Provide the initial approximation.
  double totalLikelihood = 0.0;
  for (size_t i = 0; i < kmerinds.size(); ++i) {
    const KMerStat &kms = data_[kmerinds[i]];
    const KMer &kmer = kms.kmer();
    unsigned mdist = K;
    unsigned cidx = 0;
    for (unsigned j = 0; j < l; ++j) {
      unsigned cdist = hamdistKMer(kmer, centers[j].first);
      if (cdist < mdist) {
        mdist = cdist;
        cidx = j;
      }
    }
    indices[i] = cidx;
    totalLikelihood += logLikelihoodKMer(centers[cidx].first, kms);
  }

  // Main loop
  bool changed = true, improved = true;

  // auxiliary variables
  std::vector<unsigned> dists(l);
  std::vector<double> loglike(l);
  std::vector<bool> changedCenter(l);

  while (changed && improved) {
    // fill everything with zeros
    changed = false;
    std::fill(changedCenter.begin(), changedCenter.end(), false);
    for (unsigned j = 0; j < l; ++j)
      centers[j].second.first = 0;

    double curlik = 0;

    // E step: find which clusters we belong to
    for (size_t i = 0; i < kmerinds.size(); ++i) {
      const KMerStat &kms = data_[kmerinds[i]];
      const KMer &kmer = kms.kmer();

      unsigned newInd = 0;
      for (unsigned j=0; j < l; ++j)
        loglike[j] = logLikelihoodKMer(centers[j].first, kms);
      newInd = std::max_element(loglike.begin(), loglike.end()) - loglike.begin();

      if (cfg::get().bayes_use_hamming_dist) {
        for (unsigned j=0; j < l; ++j)
          dists[j] = hamdistKMer(kmer, centers[j].first);

        newInd = std::min_element(dists.begin(), dists.end()) - dists.begin();
      }

      curlik += loglike[newInd];
      if (indices[i] != newInd) {
        changed = true;
        changedCenter[indices[i]] = true;
        changedCenter[newInd] = true;
        indices[i] = newInd;
      }
      ++centers[indices[i]].second.first;
    }

    if (cfg::get().bayes_debug_output > 1) {
#     pragma omp critical
      {
        cout << "      total likelihood=" << curlik << " as compared to previous " << totalLikelihood << endl;
      }
    }
    improved = (curlik > totalLikelihood);
    if (improved)
      totalLikelihood = curlik;

    // M step: find new cluster centers
    for (unsigned j=0; j < l; ++j) {
      if (!changedCenter[j])
        continue; // nothing has changed
      centers[j].first = ConsensusWithMask(kmerinds, indices, j);
    }
  }

  // last M step
  for (unsigned j=0; j < l; ++j) {
    centers[j].first = ConsensusWithMask(kmerinds, indices, j);
  }

  if (cfg::get().bayes_debug_output > 1) {
#   pragma omp critical
    {
      std::cout << "    final centers:\n";
      for (size_t i=0; i < centers.size(); ++i) {
        std::cout << "    " << centers[i].first << "\n";
      }
    }
  }

  return ClusterBIC(kmerinds, centers, indices);
}

size_t KMerClustering::SubClusterSingle(const std::vector<unsigned> & block, vector< vector<unsigned> > & vec) {
  size_t newkmers = 0;

  if (cfg::get().bayes_debug_output > 0) {
#   pragma omp critical
    {
      std::cout << "  kmers:\n";
      for (size_t i = 0; i < block.size(); i++) {
        cout << data_[block[i]].kmer() << '\n';
      }
    }
  }

  size_t origBlockSize = block.size();
  if (origBlockSize == 0) return 0;

  // Ad-hoc max cluster limit: we start to consider only thous k-mers which
  // multiplicity differs from maximum multiplicity by 10x.
  unsigned maxcls = 0;
  unsigned cntthr = std::max(10u, data_[block[0]].count / 10);
  for (size_t i = 0; i < block.size(); ++i)
    maxcls += (data_[block[i]].count > cntthr);
  maxcls = std::max(1u, maxcls);

  if (cfg::get().bayes_debug_output > 0) {
    #pragma omp critical
    {
      std::cout << "\nClustering an interesting block. Maximum # of clusters estimated: " << maxcls << std::endl;
    }
  }

  std::vector<unsigned> indices(origBlockSize);
  double bestLikelihood = -std::numeric_limits<double>::infinity();
  std::vector<StringCount> bestCenters;
  std::vector<unsigned> bestIndices(block.size());

  unsigned max_l = cfg::get().bayes_hammer_mode ? 1 : origBlockSize;
  for (unsigned l = 1; l <= max_l; ++l) {
    std::vector<StringCount> centers(l);
    double curLikelihood = lMeansClustering(l, block, indices, centers);
    if (cfg::get().bayes_debug_output > 0) {
      #pragma omp critical
      {
        cout << "    indices: ";
        for (uint32_t i = 0; i < origBlockSize; i++) cout << indices[i] << " ";
        cout << "\n";
        cout << "  likelihood with " << l << " clusters is " << curLikelihood << endl;
      }
    }
    if (curLikelihood > bestLikelihood) {
      bestLikelihood = curLikelihood;
      bestCenters = centers; bestIndices = indices;
    } else if (l >= maxcls)
      break;
  }

  // find if centers are in clusters
  std::vector<unsigned> centersInCluster(bestCenters.size(), -1u);
  for (unsigned i = 0; i < origBlockSize; i++) {
    unsigned dist = hamdistKMer(data_[block[i]].kmer(), bestCenters[bestIndices[i]].first);
    if (dist == 0) {
      centersInCluster[bestIndices[i]] = i;
    }
  }

  bool cons_suspicion = false;
  for (size_t k=0; k<bestCenters.size(); ++k) if (centersInCluster[k] == -1u) cons_suspicion = true;
  if (cfg::get().bayes_debug_output > 0) {
#   pragma omp critical
    {
      std::cout << "Centers: \n";
      for (size_t k=0; k<bestCenters.size(); ++k) {
        std::cout << "  " << std::setw(4) << bestCenters[k].second.first << ": ";
        if (centersInCluster[k] != -1u) {
          const KMerStat &kms = data_[block[centersInCluster[k]]];
          std::cout << kms << " " << setw(8) << block[centersInCluster[k]] << "  ";
        } else {
          std::cout << bestCenters[k].first;
        }
        std::cout << '\n';
      }
      cout << "The entire block:" << endl;
      for (uint32_t i = 0; i < origBlockSize; i++) {
        const KMerStat &kms = data_[block[i]];
        cout << "  " << kms << " " << setw(8) << block[i] << "  ";
        for (uint32_t j=0; j<K; ++j) cout << setw(3) << (unsigned)getQual(kms, j) << " "; cout << "\n";
      }
      cout << endl;
    }
  }

  // it may happen that consensus string from one subcluster occurs in other subclusters
  // we need to check for that
  bool foundBadCenter = true;
  while (foundBadCenter) {
    foundBadCenter = false;
    for (size_t k=0; k<bestCenters.size(); ++k) {
      if (foundBadCenter) break; // restart if found one bad center
      if (bestCenters[k].second.first == 0) continue;
      if (centersInCluster[k] != -1u) continue;
      for (size_t s = 0; s<bestCenters.size(); ++s) {
        if (s == k || centersInCluster[s] == -1u) continue;
        unsigned dist = hamdistKMer(bestCenters[k].first, bestCenters[s].first);
        if (dist == 0) {
          // OK, that's the situation, cluster k should be added to cluster s
          for (uint32_t i = 0; i < origBlockSize; i++) {
            if (indices[i] == k) {
              indices[i] = s;
              bestCenters[s].second.first++;
            }
          }
          bestCenters[k].second.first = 0; // it will be skipped now
          foundBadCenter = true;
          break;
        }
      }
    }
  }

  if (cfg::get().bayes_debug_output > 0 && origBlockSize > 2) {
    #pragma omp critical
    {
      cout << "\nAfter the check we got centers: \n";
      for (size_t k=0; k<bestCenters.size(); ++k) {
        cout << "  " << bestCenters[k].first.str() << " " << bestCenters[k].second << " ";
        if (centersInCluster[k] != -1u) cout << block[centersInCluster[k]];
        cout << "\n";
      }
      cout << endl;
    }
  }

  for (size_t k=0; k<bestCenters.size(); ++k) {
    if (bestCenters[k].second.first == 0) {
      continue; // superfluous cluster with zero elements
    }
    std::vector<unsigned> v;
    if (bestCenters[k].second.first == 1) {
      for (size_t i = 0; i < origBlockSize; i++) {
        if (indices[i] == k) {
          v.push_back(block[i]);
          break;
        }
      }
    } else { // there are several kmers in this cluster
      for (size_t i = 0; i < origBlockSize; i++) {
        if (bestIndices[i] == k) {
          if (centersInCluster[k] == i) {
            v.insert(v.begin(), block[i]);
          } else {
            v.push_back(block[i]);
          }
        }
      }
      if (centersInCluster[k] == -1u) {
        size_t new_idx = 0;
        #pragma omp critical
        {
          KMer newkmer = bestCenters[k].first;

          KMerStat kms(0 /* cnt */, newkmer, 1.0 /* total quality */, NULL /*quality */);
          kms.status = KMerStat::GoodIter;
          new_idx = data_.push_back(kms);
          if (data_[newkmer].kmer() != newkmer)
            newkmers += 1;
        }
        v.insert(v.begin(), new_idx);
      }
    }
    vec.push_back(v);
  }

  return newkmers;
}

static void UpdateErrors(numeric::matrix<uint64_t> &m,
                         const KMerStat &kms, const KMerStat &center) {
  const KMer &kc = center.kmer();
  const KMer &k = kms.kmer();
  for (unsigned i = 0; i < K; ++i) {
    m(kc[i], k[i]) += 1;
  }
}

class KMerStatCountComparator {
  const KMerData &data_;
public:
  KMerStatCountComparator(const KMerData &data)
      : data_(data) {}
  bool operator()(unsigned a, unsigned b) {
    return data_[a].count > data_[b].count;
  }
};

void KMerClustering::process(std::vector<std::vector<unsigned> > classes) {
  size_t newkmers = 0;
  size_t gsingl = 0, tsingl = 0, tcsingl = 0, gcsingl = 0, tcls = 0, gcls = 0, tkmers = 0, tncls = 0;

  std::ofstream ofs, ofs_bad;
  if (cfg::get().bayes_write_solid_kmers)
    ofs.open(GetGoodKMersFname());
  if (cfg::get().bayes_write_bad_kmers)
    ofs_bad.open(GetBadKMersFname());

  // Underlying code expected classes to be sorted in count decreasing order.
#if 0
  // No lambdas with gcc 4.4x :(
  std::sort(classes.begin(), classes.end(),
            [](const std::vector<unsigned> &a, const std::vector<unsigned> &b) {
              return a.size() > b.size();
            });
  for (auto I = classes.begin(), E = classes.end(); I != E; ++I) {
    std::sort(I->begin(), I->end(),
              [=](unsigned a, unsigned b) {
                return data_[a].count > data_[b].count;
              });
  }
#else
  for (auto I = classes.begin(), E = classes.end(); I != E; ++I) {
    std::sort(I->begin(), I->end(),
              KMerStatCountComparator(data_));
  }
#endif

  std::vector<numeric::matrix<uint64_t> > errs(nthreads_, numeric::matrix<double>(4, 4, 0.0));
# pragma omp parallel for shared(classes, ofs, ofs_bad, errs) num_threads(nthreads_) schedule(dynamic) reduction(+:newkmers, gsingl, tsingl, tcsingl, gcsingl, tcls, gcls, tkmers, tncls)
  for (size_t i = 0; i < classes.size(); ++i) {
    auto cur_class = classes[i];

    // No need for clustering for singletons
    if (cur_class.size() == 1) {
      size_t idx = cur_class[0];
      KMerStat &singl = data_[idx];
      if ((1-singl.totalQual) > cfg::get().bayes_singleton_threshold) {
        singl.status = KMerStat::GoodIter;
        gsingl += 1;
        if (ofs.good()) {
#         pragma omp critical
          {
            ofs << " good singleton: " << idx << "\n  " << singl << '\n';
          }
        }
      } else {
        if (cfg::get().correct_use_threshold && (1-singl.totalQual) > cfg::get().correct_threshold)
          singl.status = KMerStat::GoodIterBad;
        else
          singl.status = KMerStat::Bad;
        if (ofs_bad.good()) {
#         pragma omp critical
          {
            ofs_bad << " bad singleton: " << idx << "\n  " << singl << '\n';
          }
        }
      }
      tsingl += 1;
    } else {
      std::vector<std::vector<unsigned> > blocksInPlace;

      if (cfg::get().bayes_debug_output) {
#       pragma omp critical
        {
          std::cout << "process_SIN with block idx= " << i << " size=" << cur_class.size() << std::endl;
        }
      }
      newkmers += SubClusterSingle(cur_class, blocksInPlace);

      tncls += 1;
      for (size_t m = 0; m < blocksInPlace.size(); ++m) {
        const std::vector<unsigned> &currentBlock = blocksInPlace[m];
        if (currentBlock.size() == 0)
          continue;

        size_t cidx = currentBlock[0];
        KMerStat &center = data_[cidx];
        double center_quality = 1 - center.totalQual;

        // Computing the overall quality of a cluster.
        double cluster_quality = 1;
        if (currentBlock.size() > 1) {
          for (size_t j = 1; j < currentBlock.size(); ++j) {
            cluster_quality *= data_[currentBlock[j]].totalQual;
          }
          cluster_quality = 1-cluster_quality;
        }

        if (currentBlock.size() == 1)
          tcsingl += 1;
        else
          tcls += 1;

        if ((center_quality > cfg::get().bayes_singleton_threshold &&
             cluster_quality > cfg::get().bayes_nonsingleton_threshold) ||
            cfg::get().bayes_hammer_mode) {
          center.status = KMerStat::GoodIter;

          if (currentBlock.size() == 1)
            gcsingl += 1;
          else
            gcls += 1;

          if (ofs.good()) {
#           pragma omp critical
            {
              ofs << " center of good cluster (" << currentBlock.size() << ", " << cluster_quality << ")" << "\n  "
                  << center << '\n';
            }
          }
        } else {
          if (cfg::get().correct_use_threshold && center_quality > cfg::get().correct_threshold)
            center.status = KMerStat::GoodIterBad;
          else
            center.status = KMerStat::Bad;
          if (ofs_bad.good()) {
#           pragma omp critical
            {
              ofs_bad << " center of bad cluster (" << currentBlock.size() << ", " << cluster_quality << ")" << "\n  "
                      << center << '\n';
            }
          }
        }

        tkmers += currentBlock.size();

        for (size_t j = 1; j < currentBlock.size(); ++j) {
          size_t eidx = currentBlock[j];
          KMerStat &kms = data_[eidx];

          kms.set_change(cidx);
          UpdateErrors(errs[omp_get_thread_num()], kms, center);
          if (ofs_bad.good()) {
#           pragma omp critical
            {
              ofs_bad << " part of cluster (" << currentBlock.size() << ", " << cluster_quality << ")" << "\n  "
                      << kms << '\n';
            }
          }
        }

      }
    }
  }

  for (unsigned i = 1; i < nthreads_; ++i)
    errs[0] += errs[i];
  numeric::matrix<uint64_t> rowsums = prod(errs[0], numeric::scalar_matrix<double>(4, 1, 1));
  numeric::matrix<double> err(4, 4);
  for (unsigned i = 0; i < 4; ++i)
    for (unsigned j = 0; j < 4; ++j)
      err(i, j) = 1.0 * errs[0](i, j) / rowsums(i, 0);

  INFO("Subclustering done. Total " << newkmers << " non-read kmers were generated.");
  INFO("Subclustering statistics:");
  INFO("  Total singleton hamming clusters: " << tsingl << ". Among them " << gsingl << " (" << 100.0 * gsingl / tsingl << "%) are good");
  INFO("  Total singleton subclusters: " << tcsingl << ". Among them " << gcsingl << " (" << 100.0 * gcsingl / tcsingl << "%) are good");
  INFO("  Total non-singleton subcluster centers: " << tcls << ". Among them " << gcls << " (" << 100.0 * gcls / tcls << "%) are good");
  INFO("  Average size of non-trivial subcluster: " << 1.0 * tkmers / tcls << " kmers");
  INFO("  Average number of sub-clusters per non-singleton cluster: " << 1.0 * (tcsingl + tcls) / tncls);
  INFO("  Total solid k-mers: " << gsingl + gcsingl + gcls);
  INFO("  Substitution probabilities: " << err);
}
