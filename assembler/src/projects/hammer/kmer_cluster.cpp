//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io/reads/ireadstream.hpp"
#include "utils/parallel/openmp_wrapper.h"

#include "hammer_tools.hpp"
#include "hamcluster.hpp"
#include "kmer_cluster.hpp"
#include "config_struct_hammer.hpp"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#include <iostream>
#include <fstream>
#include <algorithm>

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

static hammer::ExpandedSeq ConsensusWithMask(const std::vector<hammer::ExpandedKMer> &kmers,
                                             const std::vector<size_t> &mask, size_t maskVal) {
  size_t block_size = kmers.size();

  // consensus of a single string is trivial
  if (block_size == 1)
    return kmers[0].seq();

  uint64_t scores[4*K] = {0};
  for (size_t j = 0; j < block_size; ++j) {
    if (mask[j] != maskVal)
      continue;

    const ExpandedSeq &kmer = kmers[j].seq();

    for (unsigned i = 0; i < K; ++i)
      scores[4*i + kmer[i]] += kmers[j].count();
  }

  hammer::ExpandedSeq res;
  for (unsigned i = 0; i < K; ++i)
    res[i] = (char)(std::max_element(scores + 4*i, scores + 4*i + 4) - (scores + 4*i));

  return res;
}

static hammer::ExpandedSeq Consensus(const std::vector<hammer::ExpandedKMer> &kmers) {
  size_t block_size = kmers.size();

  // consensus of a single string is trivial
  if (block_size == 1)
    return kmers[0].seq();

  uint64_t scores[4*K] = {0};
  for (size_t j = 0; j < block_size; ++j) {
    const ExpandedSeq &kmer = kmers[j].seq();

    for (unsigned i = 0; i < K; ++i)
      scores[4*i + kmer[i]] += kmers[j].count();
  }

  hammer::ExpandedSeq res;
  for (unsigned i = 0; i < K; ++i)
    res[i] = (char)(std::max_element(scores + 4*i, scores + 4*i + 4) - (scores + 4*i));

  return res;
}

double KMerClustering::ClusterBIC(const std::vector<Center> &centers,
                                  const std::vector<size_t> &indices, const std::vector<hammer::ExpandedKMer> &kmers) const {
  size_t block_size = indices.size();
  size_t clusters = centers.size();
  if (block_size == 0)
    return -std::numeric_limits<double>::infinity();
  assert(centers.size() > 0);

  double loglik = 0;
  unsigned total = 0;
  for (size_t i = 0; i < block_size; ++i) {
    loglik += kmers[i].count()*kmers[i].logL(centers[indices[i]].center_);
    total += kmers[i].count();
  }

  size_t nparams = (clusters - 1) + clusters*K + 2*clusters*K;

  if (cfg::get().bayes_debug_output > 1) {
#   pragma omp critical
    {
        std::cout << "  logL: " << loglik << ", clusters: " << clusters << ", nparams: " << nparams << ", N: " << block_size << std::endl;
    }
  }
  
  return loglik - (double)nparams * log((double)total) / 2.0;
}


double KMerClustering::lMeansClustering(unsigned l, const std::vector<hammer::ExpandedKMer> &kmers,
                                        std::vector<size_t> &indices, std::vector<Center> &centers) {
  centers.resize(l); // there are l centers

  // if l==1 then clustering is trivial
  if (l == 1) {
    centers[0].center_ = Consensus(kmers);
    centers[0].count_ = kmers.size();
    for (size_t i = 0; i < kmers.size(); ++i)
      indices[i] = 0;
    return ClusterBIC(centers, indices, kmers);
  }

  // Provide the initial approximation.
  double totalLikelihood = 0.0;
  if (cfg::get().bayes_initial_refine) {
    // Refine the current approximation
    centers[l-1].center_ = kmers[l-1].seq();
    for (size_t i = 0; i < kmers.size(); ++i) {
      size_t cidx = indices[i];
      unsigned cdist = kmers[i].hamdist(centers[cidx].center_, K);
      unsigned mdist = kmers[i].hamdist(centers[l-1].center_, cdist);
      if (mdist < cdist) {
        indices[i] = l - 1;
        cidx = l - 1;
      }
      totalLikelihood += kmers[i].logL(centers[cidx].center_);
    }
  } else {
    // We assume that kmers are sorted wrt the count.
    for (size_t j = 0; j < l; ++j)
      centers[j].center_ = kmers[j].seq();

    for (size_t i = 0; i < kmers.size(); ++i) {
      unsigned mdist = K;
      unsigned cidx = 0;
      for (unsigned j = 0; j < l; ++j) {
        unsigned cdist = kmers[i].hamdist(centers[j].center_, mdist);
        if (cdist < mdist) {
          mdist = cdist;
          cidx = j;
        }
      }
      indices[i] = cidx;
      totalLikelihood += kmers[i].logL(centers[cidx].center_);
    }
  }

  if (cfg::get().bayes_debug_output > 1) {
#   pragma omp critical
    {
      std::cout << "    centers:\n";
      for (size_t i=0; i < centers.size(); ++i) {
        std::cout << "    " << centers[i].center_ << "\n";
      }
    }
  }

  // Main loop
  bool changed = true, improved = true;

  // auxiliary variables
  std::vector<size_t> dists(l);
  std::vector<double> loglike(l);
  std::vector<bool> changedCenter(l);

  while (changed && improved) {
    // fill everything with zeros
    changed = false;
    std::fill(changedCenter.begin(), changedCenter.end(), false);
    for (unsigned j = 0; j < l; ++j)
      centers[j].count_ = 0;

    double curlik = 0;

    // E step: find which clusters we belong to
    for (size_t i = 0; i < kmers.size(); ++i) {
      size_t newInd = 0;
      if (cfg::get().bayes_use_hamming_dist) {
        for (unsigned j = 0; j < l; ++j)
          dists[j] = kmers[i].hamdist(centers[j].center_);

        newInd = std::min_element(dists.begin(), dists.end()) - dists.begin();
      } else {
        for (unsigned j = 0; j < l; ++j)
          loglike[j] = kmers[i].logL(centers[j].center_);
        newInd = std::max_element(loglike.begin(), loglike.end()) - loglike.begin();
      }

      curlik += loglike[newInd];
      if (indices[i] != newInd) {
        changed = true;
        changedCenter[indices[i]] = true;
        changedCenter[newInd] = true;
        indices[i] = newInd;
      }
      ++centers[indices[i]].count_;
    }

    if (cfg::get().bayes_debug_output > 1) {
#     pragma omp critical
      {
        std::cout << "      total likelihood=" << curlik << " as compared to previous " << totalLikelihood << std::endl;
      }
    }
    improved = (curlik > totalLikelihood);
    if (improved)
      totalLikelihood = curlik;

    // M step: find new cluster centers
    for (unsigned j=0; j < l; ++j) {
      if (!changedCenter[j])
        continue; // nothing has changed

      centers[j].center_ = ConsensusWithMask(kmers, indices, j);
    }
  }

  // last M step
  for (unsigned j=0; j < l; ++j)
    centers[j].center_ = ConsensusWithMask(kmers, indices, j);

  if (cfg::get().bayes_debug_output > 1) {
#   pragma omp critical
    {
      std::cout << "    final centers:\n";
      for (size_t i=0; i < centers.size(); ++i) {
        std::cout << "    " << centers[i].center_ << "\n";
      }
    }
  }

  return ClusterBIC(centers, indices, kmers);
}


size_t KMerClustering::SubClusterSingle(const std::vector<size_t> & block, std::vector< std::vector<size_t> > & vec) {
  size_t newkmers = 0;

  if (cfg::get().bayes_debug_output > 0) {
#   pragma omp critical
    {
      std::cout << "  kmers:\n";
      for (size_t i = 0; i < block.size(); i++) {
        std::cout << data_.kmer(block[i]) << '\n';
      }
    }
  }

  size_t origBlockSize = block.size();
  if (origBlockSize == 0) return 0;

  // Ad-hoc max cluster limit: we start to consider only those k-mers which
  // multiplicity differs from maximum multiplicity by 10x.
  size_t maxcls = 0;
  size_t cntthr = std::max(10u, data_[block[0]].count() / 10);
  for (size_t i = 0; i < block.size(); ++i)
    maxcls += (data_[block[i]].count() > cntthr);
  // Another limit: we're interested in good centers only
  size_t maxgcnt = 0;
  for (size_t i = 0; i < block.size(); ++i) {
    float center_quality = 1 - data_[block[i]].total_qual;
    if ((center_quality > cfg::get().bayes_singleton_threshold) ||
        (cfg::get().correct_use_threshold && center_quality > cfg::get().correct_threshold))
      maxgcnt += 1;
  }
  
  maxcls = std::min(maxcls, maxgcnt) + 1;
  if (cfg::get().bayes_debug_output > 0) {
    #pragma omp critical
    {
      std::cout << "\nClustering an interesting block. Maximum # of clusters estimated: " << maxcls << std::endl;
    }
  }

  // Prepare the expanded k-mer structure
  std::vector<hammer::ExpandedKMer> kmers;
  for (size_t idx : block)
    kmers.emplace_back(data_.kmer(idx), data_[idx]);

  double bestLikelihood = -std::numeric_limits<double>::infinity();
  std::vector<Center> bestCenters;
  std::vector<size_t> indices(origBlockSize);
  std::vector<size_t> bestIndices(origBlockSize);

  unsigned max_l = cfg::get().bayes_hammer_mode ? 1 : (unsigned) origBlockSize;
  std::vector<Center> centers;
  for (unsigned l = 1; l <= max_l; ++l) {
    double curLikelihood = lMeansClustering(l, kmers, indices, centers);
    if (cfg::get().bayes_debug_output > 0) {
      #pragma omp critical
      {
        std::cout << "    indices: ";
        for (uint32_t i = 0; i < origBlockSize; i++) std::cout << indices[i] << " ";
        std::cout << "\n";
        std::cout << "  likelihood with " << l << " clusters is " << curLikelihood << std::endl;
      }
    }
    if (curLikelihood > bestLikelihood) {
      bestLikelihood = curLikelihood;
      bestCenters = centers; bestIndices = indices;
    } else if (l >= maxcls)
      break;
  }

  // find if centers are in clusters
  std::vector<size_t> centersInCluster(bestCenters.size(), -1u);
  for (unsigned i = 0; i < origBlockSize; i++) {
    unsigned dist = kmers[i].hamdist(bestCenters[bestIndices[i]].center_);
    if (dist == 0)
      centersInCluster[bestIndices[i]] = i;
  }

  if (cfg::get().bayes_debug_output > 0) {
#   pragma omp critical
    {
      std::cout << "Centers: \n";
      for (size_t k=0; k<bestCenters.size(); ++k) {
        std::cout << "  " << std::setw(4) << bestCenters[k].count_ << ": ";
        if (centersInCluster[k] != -1u) {
          const KMerStat &kms = data_[block[centersInCluster[k]]];
          std::cout << kms << " " << std::setw(8) << block[centersInCluster[k]] << "  ";
        } else {
          std::cout << bestCenters[k].center_;
        }
        std::cout << '\n';
      }
      std::cout << "The entire block:" << std::endl;
      for (uint32_t i = 0; i < origBlockSize; i++) {
        const KMerStat &kms = data_[block[i]];
        std::cout << "  " << kms << " " << std::setw(8) << block[i] << "  ";
        for (uint32_t j=0; j<K; ++j) std::cout << std::setw(3) << (unsigned)getQual(kms, j) << " "; std::cout << "\n";
      }
      std::cout << std::endl;
    }
  }

  // it may happen that consensus string from one subcluster occurs in other subclusters
  // we need to check for that
  bool foundBadCenter = true;
  while (foundBadCenter) {
    foundBadCenter = false;
    for (size_t k=0; k<bestCenters.size(); ++k) {
      if (foundBadCenter) break; // restart if found one bad center
      if (bestCenters[k].count_ == 0) continue;
      if (centersInCluster[k] != -1u) continue;
      for (size_t s = 0; s< bestCenters.size(); ++s) {
        if (s == k || centersInCluster[s] == -1u) continue;
        unsigned dist = hamdist(bestCenters[k].center_, bestCenters[s].center_);
        if (dist == 0) {
          // OK, that's the situation, cluster k should be added to cluster s
          for (uint32_t i = 0; i < origBlockSize; i++) {
            if (indices[i] == k) {
              indices[i] = (unsigned)s;
              bestCenters[s].count_++;
            }
          }
          bestCenters[k].count_ = 0; // it will be skipped now
          foundBadCenter = true;
          break;
        }
      }
    }
  }

  if (cfg::get().bayes_debug_output > 0 && origBlockSize > 2) {
    #pragma omp critical
    {
      std::cout << "\nAfter the check we got centers: \n";
      for (size_t k=0; k<bestCenters.size(); ++k) {
        std::cout << "  " << bestCenters[k].center_ << " (" << bestCenters[k].count_ << ")";
        if (centersInCluster[k] != -1u) std::cout << block[centersInCluster[k]];
        std::cout << "\n";
      }
      std::cout << std::endl;
    }
  }

  for (size_t k = 0; k < bestCenters.size(); ++k) {
    if (bestCenters[k].count_ == 0)
      continue; // superfluous cluster with zero elements

    std::vector<size_t> v;
    if (bestCenters[k].count_ == 1) {
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
        KMer newkmer(bestCenters[k].center_);
        size_t new_idx = data_.checking_seq_idx(newkmer);
        if (new_idx == -1ULL) {
          #pragma omp critical
          {
            KMerStat kms(0 /* cnt */, 1.0 /* total quality */, NULL /*quality */);
            kms.mark_good();
            new_idx = data_.push_back(newkmer, kms);
            newkmers += 1;
          }
        }
        v.insert(v.begin(), new_idx);
      }
    }
    vec.push_back(v);
  }

  return newkmers;
}

static void UpdateErrors(numeric::matrix<uint64_t> &m,
                         const KMer k, const KMer kc) {
  for (unsigned i = 0; i < K; ++i) {
    m(kc[i], k[i]) += 1;
  }
}

size_t KMerClustering::ProcessCluster(const std::vector<size_t> &cur_class,
                                      numeric::matrix<uint64_t> &errs,
                                      std::ofstream &ofs, std::ofstream &ofs_bad,
                                      size_t &gsingl, size_t &tsingl, size_t &tcsingl, size_t &gcsingl,
                                      size_t &tcls, size_t &gcls, size_t &tkmers, size_t &tncls) {
    size_t newkmers = 0;

    // No need for clustering for singletons
    if (cur_class.size() == 1) {
        size_t idx = cur_class[0];
        KMerStat &singl = data_[idx];
        if ((1-singl.total_qual) > cfg::get().bayes_singleton_threshold) {
            singl.mark_good();
            gsingl += 1;

            if (ofs.good()) {
#               pragma omp critical
                {
                    ofs << " good singleton: " << idx << "\n  " << singl << '\n';
                }
            }
        } else {
            if (cfg::get().correct_use_threshold && (1-singl.total_qual) > cfg::get().correct_threshold)
                singl.mark_good();
            else
                singl.mark_bad();

            if (ofs_bad.good()) {
#               pragma omp critical
                {
                    ofs_bad << " bad singleton: " << idx << "\n  " << singl << '\n';
                }
            }
        }
        tsingl += 1;
        return 0;
    }

    std::vector<std::vector<size_t> > blocksInPlace;
    if (cfg::get().bayes_debug_output) {
#       pragma omp critical
        {
          std::cout << "process_SIN with size=" << cur_class.size() << std::endl;
        }
      }
    newkmers += SubClusterSingle(cur_class, blocksInPlace);

    tncls += 1;
    for (size_t m = 0; m < blocksInPlace.size(); ++m) {
        const std::vector<size_t> &currentBlock = blocksInPlace[m];
        if (currentBlock.size() == 0)
            continue;

        size_t cidx = currentBlock[0];
        KMerStat &center = data_[cidx];
        KMer ckmer = data_.kmer(cidx);
        double center_quality = 1 - center.total_qual;

        // Computing the overall quality of a cluster.
        double cluster_quality = 1;
        if (currentBlock.size() > 1) {
            for (size_t j = 1; j < currentBlock.size(); ++j)
                cluster_quality *= data_[currentBlock[j]].total_qual;

            cluster_quality = 1-cluster_quality;
        }

        if (currentBlock.size() == 1)
            tcsingl += 1;
        else
            tcls += 1;

        if ((center_quality > cfg::get().bayes_singleton_threshold &&
             cluster_quality > cfg::get().bayes_nonsingleton_threshold) ||
            cfg::get().bayes_hammer_mode) {
          center.mark_good();

          if (currentBlock.size() == 1)
              gcsingl += 1;
          else
              gcls += 1;

          if (ofs.good()) {
#             pragma omp critical
              {
                  ofs << " center of good cluster (" << currentBlock.size() << ", " << cluster_quality << ")" << "\n  "
                  << center << '\n';
              }
          }
        } else {
            if (cfg::get().correct_use_threshold && center_quality > cfg::get().correct_threshold)
                center.mark_good();
            else
                center.mark_bad();
            if (ofs_bad.good()) {
#               pragma omp critical
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

            UpdateErrors(errs, data_.kmer(eidx), ckmer);

            if (ofs_bad.good()) {
#               pragma omp critical
                {
                    ofs_bad << " part of cluster (" << currentBlock.size() << ", " << cluster_quality << ")" << "\n  "
                            << kms << '\n';
                }
            }
        }
    }

    return newkmers;
}


class KMerStatCountComparator {
  const KMerData &data_;
public:
  KMerStatCountComparator(const KMerData &data)
      : data_(data) {}
  bool operator()(size_t a, size_t b) {
    return data_[a].count() > data_[b].count();
  }
};

void KMerClustering::process(const std::string &Prefix) {
  size_t newkmers = 0;
  size_t gsingl = 0, tsingl = 0, tcsingl = 0, gcsingl = 0, tcls = 0, gcls = 0, tkmers = 0, tncls = 0;

  std::ofstream ofs, ofs_bad;
  if (cfg::get().bayes_write_solid_kmers)
    ofs.open(GetGoodKMersFname());
  if (cfg::get().bayes_write_bad_kmers)
    ofs_bad.open(GetBadKMersFname());

  // Open and read index file
  MMappedRecordReader<size_t> findex(Prefix + ".idx",  /* unlink */ !debug_, -1ULL);

  std::vector<numeric::matrix<uint64_t> > errs(nthreads_, numeric::matrix<double>(4, 4, 0.0));

# pragma omp parallel for shared(ofs, ofs_bad, errs) num_threads(nthreads_) schedule(guided) reduction(+:newkmers, gsingl, tsingl, tcsingl, gcsingl, tcls, gcls, tkmers, tncls)
  for (size_t chunk = 0; chunk < nthreads_ * nthreads_; ++chunk) {
      size_t *current = findex.data() + findex.size() * chunk / nthreads_ / nthreads_;
      size_t *next = findex.data() + findex.size() * (chunk + 1)/ nthreads_ / nthreads_;
      std::ifstream is(Prefix, std::ios::in | std::ios::binary);

      // Calculate how much we need to seek
      size_t soff = 0;
      for (size_t *csz = findex.data(); csz != current; ++csz)
          soff += *csz;

      // Now see the stream and start processing
      is.seekg(soff * sizeof(size_t));

      for (; current != next; ++current) {
          std::vector<size_t> cluster(*current);
          VERIFY(is.good());
          is.read((char*)&cluster[0], *current * sizeof(cluster[0]));

          // Underlying code expected classes to be sorted in count decreasing order.
          std::sort(cluster.begin(), cluster.end(), KMerStatCountComparator(data_));

          newkmers += ProcessCluster(cluster,
                                     errs[omp_get_thread_num()],
                                     ofs, ofs_bad,
                                     gsingl, tsingl, tcsingl, gcsingl,
                                     tcls, gcls, tkmers, tncls);
      }
  }

  if (!debug_) {
      int res = unlink(Prefix.c_str());
      VERIFY_MSG(res == 0,
                 "unlink(2) failed. Reason: " << strerror(errno) << ". Error code: " << errno);
  }

  for (unsigned i = 1; i < nthreads_; ++i)
    errs[0] += errs[i];

  numeric::matrix<uint64_t> rowsums = prod(errs[0], numeric::scalar_matrix<double>(4, 1, 1));
  numeric::matrix<double> err(4, 4);
  for (unsigned i = 0; i < 4; ++i)
    for (unsigned j = 0; j < 4; ++j)
      err(i, j) = 1.0 * (double)errs[0](i, j) / (double)rowsums(i, 0);

  INFO("Subclustering done. Total " << newkmers << " non-read kmers were generated.");
  INFO("Subclustering statistics:");
  INFO("  Total singleton hamming clusters: " << tsingl << ". Among them " << gsingl << " (" << 100.0 * (double)gsingl / (double)tsingl << "%) are good");
  INFO("  Total singleton subclusters: " << tcsingl << ". Among them " << gcsingl << " (" << 100.0 * (double)gcsingl / (double)tcsingl << "%) are good");
  INFO("  Total non-singleton subcluster centers: " << tcls << ". Among them " << gcls << " (" << 100.0 * (double)gcls / (double)tcls << "%) are good");
  INFO("  Average size of non-trivial subcluster: " << 1.0 * (double)tkmers / (double)tcls << " kmers");
  INFO("  Average number of sub-clusters per non-singleton cluster: " << 1.0 * (double)(tcsingl + tcls) / (double)tncls);
  INFO("  Total solid k-mers: " << gsingl + gcsingl + gcls);
  INFO("  Substitution probabilities: " << err);
}
