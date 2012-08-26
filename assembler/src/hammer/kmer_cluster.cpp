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

using std::max_element;
using std::min_element;

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
      res += - log(10) * getQual(x, i) / 10.0;
    } else {
      res += getProb(x, i, /* log */ true);
    }
  }

  return res;
}

KMer KMerClustering::find_consensus_with_mask(const vector<unsigned> & cl, const vector<unsigned> & mask, unsigned maskVal) {
  size_t blockSize = cl.size();

  // consensus of a single string is trivial
  if (blockSize == 1) return data_[cl[0]].kmer();

  std::string kmer(K, 'A');
  for (unsigned i = 0; i < K; i++) {
    int scores[4] = {0,0,0,0};
    for (uint32_t j = 0; j < blockSize; j++) {
      if (mask[j] == maskVal)
        scores[static_cast<uint8_t>(data_[cl[j]].kmer()[i])] += data_[cl[j]].count;
    }
    kmer[i] = nucl(std::max_element(scores, scores + 4) - scores);
  }

  return KMer(kmer);
}


KMer KMerClustering::find_consensus(const vector<unsigned> & cl) {
  size_t blockSize = cl.size();

  // consensus of a single string is trivial
  if (blockSize == 1) return data_[cl[0]].kmer();

  std::string kmer(K, 'A');
  for (size_t i = 0; i < K; i++) {
    int scores[4] = {0,0,0,0};
    for (size_t j = 0; j < blockSize; j++) {
      scores[static_cast<uint8_t>(data_[cl[j]].kmer()[i])]+=data_[cl[j]].count;
    }
    kmer[i] = nucl(std::max_element(scores, scores + 4) - scores);
  }
  return KMer(kmer);
}

double KMerClustering::trueClusterLogLikelihood(const vector<unsigned> & cl, const vector<StringCount> & centers, const vector<unsigned> & indices) {
  size_t blockSize = cl.size();
  size_t clusters = centers.size();
  if (blockSize == 0)
    return -std::numeric_limits<double>::infinity();
  assert(blockSize == indices.size());
  assert(centers.size() > 0);

  double loglik = 0;
  unsigned total = 0;
  for (size_t i=0; i<blockSize; ++i) {
    unsigned cnt = data_[cl[i]].count;
    loglik += cnt*logLikelihoodKMer(centers[indices[i]].first, data_[cl[i]]);
    total += cnt;
  }

  unsigned nparams = (clusters - 1) + clusters*K + 4*clusters*K;

  return loglik - nparams*log(total) / 2;
}

double KMerClustering::lMeansClustering(unsigned l, const std::vector<unsigned> &kmerinds,
                                        std::vector<unsigned> & indices, std::vector<StringCount> & centers) {
  centers.resize(l); // there are l centers

  // if l==1 then clustering is trivial
  if (l == 1) {
    centers[0].first = find_consensus(kmerinds);
    centers[0].second.first = kmerinds.size();
    centers[0].second.second = 1;
    for (size_t i = 0; i < kmerinds.size(); ++i)
      indices[i] = 0;
    return trueClusterLogLikelihood(kmerinds, centers, indices);
  }

  for (unsigned j = 0; j < l; ++j) {
    centers[j].first = data_[kmerinds[j]].kmer();
    centers[j].second.first = data_[kmerinds[j]].count;
    centers[j].second.second = data_[kmerinds[j]].totalQual;
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

  // auxiliary variables
  std::vector<int> dists(l);
  std::vector<double> loglike(l);
  std::vector<bool> changedCenter(l);

  // main loop
  bool changed = true, improved = true;
  double totalLikelihood = -std::numeric_limits<double>::infinity();

  std::vector<StringCount> bestCenters;
  std::vector<unsigned> bestIndices(kmerinds.size());

  while (changed && improved) {
    // fill everything with zeros
    changed = false;
    std::fill(changedCenter.begin(), changedCenter.end(), false);
    for (unsigned j=0; j < l; ++j)
      centers[j].second.first = 0;

    double cur_total = 0;

    // E step: find which clusters we belong to
    for (size_t i=0; i < kmerinds.size(); ++i) {
      const KMerStat &kms = data_[kmerinds[i]];
      const KMer &kmer = kms.kmer();
      for (uint32_t j=0; j < l; ++j) {
        dists[j] = hamdistKMer(kmer, centers[j].first);
        loglike[j] = logLikelihoodKMer(centers[j].first, kms);
      }

      if (cfg::get().bayes_debug_output > 1) {
#       pragma omp critical
        {
          cout << "      likelihoods for " << i << ": ";
          for (size_t j=0; j < l; ++j) {
            cout << loglike[j] << " ";
          }
          cout << endl;
        }
      }
      unsigned newInd = cfg::get().bayes_use_hamming_dist ?
                        (std::min_element(dists.begin(), dists.end()) - dists.begin()) :
                        (std::max_element(loglike.begin(), loglike.end()) - loglike.begin());
      cur_total += loglike[newInd];
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
        cout << "      total likelihood=" << cur_total << " as compared to previous " << totalLikelihood << endl;
      }
    }
    improved = (cur_total > totalLikelihood);
    if (improved)
      totalLikelihood = cur_total;

    // M step: find new cluster centers
    for (unsigned j=0; j < l; ++j) {
      if (!changedCenter[j])
        continue; // nothing has changed
      centers[j].first = find_consensus_with_mask(kmerinds, indices, j);
    }
  }

  // last M step
  for (unsigned j=0; j < l; ++j) {
    centers[j].first = find_consensus_with_mask(kmerinds, indices, j);
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

  return trueClusterLogLikelihood(kmerinds, centers, indices);
}

size_t KMerClustering::process_block_SIN(const std::vector<unsigned> & block, vector< vector<unsigned> > & vec) {
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

  unsigned maxcls = 0;
  for (size_t i = 0; i < block.size(); ++i) {
    maxcls += data_[block[i]].count > 10;
  }
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
    } else if (l > maxcls)
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
  for (size_t k=0; k<bestCenters.size(); ++k) if (centersInCluster[k] == -1) cons_suspicion = true;
  if (cfg::get().bayes_debug_output > 0) {
#   pragma omp critical
    {
      std::cout << "Centers: \n";
      for (size_t k=0; k<bestCenters.size(); ++k) {
        std::cout << "  " << bestCenters[k].second.first << ": ";
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
      if (centersInCluster[k] == -1) {
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

static void UpdateErrors(boost::numeric::ublas::matrix<uint64_t> &m,
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

  std::vector<boost::numeric::ublas::matrix<uint64_t> > errs(nthreads_, boost::numeric::ublas::matrix<double>(4, 4));
# pragma omp parallel for shared(classes, ofs, ofs_bad, errs) num_threads(nthreads_) schedule(dynamic) reduction(+:newkmers, gsingl, tsingl, tcsingl, gcsingl, tcls, gcls, tkmers, tncls)
  for (size_t i = 0; i < classes.size(); ++i) {
    //  {
    // size_t i = 8;
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
      newkmers += process_block_SIN(cur_class, blocksInPlace);

      tncls += 1;
      for (size_t m = 0; m < blocksInPlace.size(); ++m) {
        if (blocksInPlace[m].size() == 0)
          continue;

        if (blocksInPlace[m].size() == 1) {
          size_t idx = blocksInPlace[m][0];
          KMerStat &singl = data_[idx];
          if ((1-singl.totalQual) > cfg::get().bayes_singleton_threshold) {
            singl.status = KMerStat::GoodIter;
            gcsingl += 1;
            if (ofs.good()) {
#             pragma omp critical
              {
                ofs << " good cluster singleton: " << idx << "\n  " << singl << '\n';
              }
            }
          } else {
            if (cfg::get().correct_use_threshold && (1-singl.totalQual) > cfg::get().correct_threshold)
              singl.status = KMerStat::GoodIterBad;
            else
              singl.status = KMerStat::Bad;
            if (ofs_bad.good()) {
#             pragma omp critical
              {
                ofs_bad << " bad cluster singleton: " << idx << "\n  " << singl << '\n';
              }
            }
          }
          tcsingl += 1;
        } else {
          size_t cidx = blocksInPlace[m][0];
          KMerStat &center = data_[cidx];

          // we've got a nontrivial cluster; computing its overall quality
          double cluster_quality = 1;
          for (size_t j = 1; j < blocksInPlace[m].size(); ++j) {
            cluster_quality *= data_[blocksInPlace[m][j]].totalQual;
          }
          cluster_quality = 1-cluster_quality;

          // in regular hammer mode, all nonsingletons are good
          if (cfg::get().bayes_hammer_mode ||
              cluster_quality > cfg::get().bayes_nonsingleton_threshold) {
            center.status = KMerStat::GoodIter;
            gcls += 1;
            if (ofs.good()) {
#           pragma omp critical
              {
                ofs << " center of good cluster (" << blocksInPlace[m].size() - 1 << ", " << cluster_quality << ")" << "\n  "
                    << center << '\n';
              }
            }
          } else {
            if (cfg::get().correct_use_threshold && (1-center.totalQual) > cfg::get().correct_threshold)
              center.status = KMerStat::GoodIterBad;
            else
              center.status = KMerStat::Bad;
            if (ofs_bad.good()) {
#           pragma omp critical
              {
                ofs_bad << " center of bad cluster (" << blocksInPlace[m].size() - 1 << ", " << cluster_quality << ")" << "\n  "
                        << center << '\n';
              }
            }
          }
          tcls += 1;
          tkmers += (blocksInPlace[m].size() - 1);

          for (size_t j = 1; j < blocksInPlace[m].size(); ++j) {
            size_t eidx = blocksInPlace[m][j];
            KMerStat &kms = data_[eidx];

            kms.set_change(cidx);
            UpdateErrors(errs[omp_get_thread_num()], kms, center);
            if (ofs_bad.good()) {
#           pragma omp critical
              {
              ofs_bad << " part of cluster (" << blocksInPlace[m].size() - 1 << ", " << cluster_quality << ")" << "\n  "
                      << kms << '\n';
              }
            }
          }
        }
      }
    }
  }

  for (unsigned i = 1; i < nthreads_; ++i)
    errs[0] += errs[i];
  boost::numeric::ublas::matrix<uint64_t> rowsums = prod(errs[0], boost::numeric::ublas::scalar_matrix<double>(4, 1, 1));
  boost::numeric::ublas::matrix<double> err(4, 4);
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
