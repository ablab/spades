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

struct entry_logger
{
	entry_logger(std::string const& name)
		: name_(name)
		, tid_ (syscall(SYS_gettid))
	{
		#pragma omp critical
		{
			std::cout << ">> entered " << name_ << "; thread " << tid_ << std::endl;
		}
	}

	~entry_logger()
	{
		#pragma omp critical
		{
			std::cout << ">> finished " << name_ << "; thread " << tid_ << std::endl;
		}
	}

private:
	string name_;
	pid_t  tid_;
};

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

static double logLikelihoodSingleton(const KMerStat & x) {
	double res = 0;
	for (unsigned i = 0; i < K; ++i) {
		res += getProb(x, i, /* log */ true);
	}
	return res;
}

KMer KMerClustering::find_consensus_with_mask(const vector<unsigned> & cl, const vector<int> & mask, int maskVal) {
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

/**
  * @return total log-likelihood of this particular clustering with real quality values
  */
double KMerClustering::trueClusterLogLikelihood(const vector<unsigned> & cl, const vector<StringCount> & centers, const vector<int> & indices) {
	size_t blockSize = cl.size();
	if (blockSize == 0) return 0.0;
	assert(blockSize == indices.size());
	assert(centers.size() > 0);

	// if there is only one center, there are no beta coefficients
	if (centers.size() == 1) {
		double logLikelihood = 0;
		for (size_t i=0; i<blockSize; ++i) {
			logLikelihood += logLikelihoodKMer(centers[indices[i]].first, data_[cl[i]]);
		}
		return ( lMultinomial(cl, data_) + logLikelihood );
	}

	// compute sufficient statistics
	vector<int> count(centers.size(), 0);		// how many kmers in cluster i
	vector<double> totalLogLikelihood(centers.size(), 0);	// total distance from kmers of cluster i to its center
	for (size_t i=0; i<blockSize; ++i) {
		count[indices[i]]+=data_[cl[i]].count;
		totalLogLikelihood[indices[i]] += logLikelihoodKMer(centers[indices[i]].first, data_[cl[i]]);
	}

	// sum up the likelihood
	double res = lBetaPlusOne(count);   // 1/B(count)
	res += lMultinomial(centers); 		// {sum(centers.count) \choose centers.count}
	for (size_t i=0; i<centers.size(); ++i) {
		res += lMultinomialWithMask(cl, data_, indices, i) + totalLogLikelihood[i];
	}
	res -= logSimplexVolume(centers.size());
	return res;
}

/**
  * @return total log-likelihood of this particular clustering with real quality values
  */
double KMerClustering::trueSingletonLogLikelihood(size_t kmind) {
	return logLikelihoodSingleton(data_[kmind]);
}


double KMerClustering::lMeansClustering(uint32_t l, const KMerHamDistMatrix &distances, const vector<unsigned> & kmerinds, vector<int> & indices, vector<StringCount> & centers) {
	centers.resize(l); // there are l centers

	// if l==1 then clustering is trivial
	if (l == 1) {
		centers[0].first = find_consensus(kmerinds);
		centers[0].second.first = kmerinds.size();
		centers[0].second.second = 1;
		for (size_t i=0; i < kmerinds.size(); ++i) indices[i] = 0;
		return trueClusterLogLikelihood(kmerinds, centers, indices);
	}

	int restartCount = 1;
	
	double bestLikelihood = -1000000.0;
	double curLikelihood = -1000000.0;
	vector<StringCount> bestCenters;
	vector<int> bestIndices(kmerinds.size());
	
	for (int restart = 0; restart < restartCount; ++restart) {

	// fill in initial approximations
	// idea: get centers at maximum distances from each other (with a greedy algorithm, of course)

	vector<uint32_t> init(l);
	for (uint32_t j=0; j < l; ++j) centers[j].second.first = 0;
	pair<uint32_t, uint32_t> fp(0, 1); // first pair: furthest, then count, then quality
	for (uint32_t i=0; i<kmerinds.size(); ++i) {
		for (uint32_t j=i+1; j<kmerinds.size(); ++j) {
			if (distances(i, j) > distances(fp.first, fp.second) ||
          (distances(i, j) == distances(fp.first, fp.second) &&
           ( data_[kmerinds[i]].count + data_[kmerinds[j]].count >
             data_[kmerinds[fp.first]].count + data_[kmerinds[fp.second]].count ||
             ( data_[kmerinds[i]].count + data_[kmerinds[j]].count ==
               data_[kmerinds[fp.first]].count + data_[kmerinds[fp.second]].count &&
               data_[kmerinds[i]].totalQual + data_[kmerinds[j]].totalQual >
               data_[kmerinds[fp.first]].totalQual + data_[kmerinds[fp.second]].totalQual)))) {
				fp.first = i; fp.second = j;
			}
		}
	}
	if (cfg::get().bayes_debug_output) {
		#pragma omp critical
		{
		cout << "    fp.first=" << fp.first << "  fp.second=" << fp.second << "\n";
		}
	}
	init[0] = fp.first; init[1] = fp.second;
	for (uint32_t j=2; j<l; ++j) {
		// to find each next center candidate, maximize total distance to previous centers
		int bestDist = -1; int ind = -1;
		for (uint32_t i=0; i < kmerinds.size(); ++i) {
			// check if it is a new k-mer
			bool newIndex = true;
			for (uint32_t k=0; k<j; ++k) if (i == init[k]) { newIndex = false; break; }
			if (!newIndex) continue;

			int curDist = 0;
			for (uint32_t k=0; k<j; ++k) curDist += distances(i, init[k]);

			if ( curDist > bestDist ) {
				ind = i; bestDist = curDist;
				continue;
			}
		}
		init[j] = ind;
	}

	for (uint32_t j=0; j<l; ++j) {
		centers[j].first = data_[kmerinds[init[j]]].kmer();
		centers[j].second.first = data_[kmerinds[init[j]]].count;
		centers[j].second.second = data_[kmerinds[init[j]]].totalQual;
	}


	if (cfg::get().bayes_debug_output) {
		#pragma omp critical
		{
		cout << "    centers:\n";
		for (size_t i=0; i < centers.size(); ++i) {
			cout << "    " << centers[i].first << "\n";
		}
		}
	}

	// TODO: make random restarts better!!! they don't really work now
	if (restart > 0) { // introduce random noise
		vector<bool> good(kmerinds.size(), false);
		for (size_t i=0; i < kmerinds.size(); ++i) {
			for (uint32_t j=0; j < l; ++j) {
				if (data_[kmerinds[i]].count >= centers[j].second.first ||
						( data_[kmerinds[i]].count == centers[j].second.first &&
						  data_[kmerinds[i]].totalQual > centers[j].second.second )) {
					good[i] = true; break;
				}
			}
		}
		for (uint32_t k=0; k < l/2 + 1; ++k) {
			int newNo = rand() % (kmerinds.size() - l - k);
			int indexOld = (rand() % (l-1)) + 1;
			int indexNew = kmerinds.size()-1;
			for (size_t i=0; i < kmerinds.size(); ++i) {
				if (!good[i]) --newNo;
				if (newNo < 0) { indexNew = i; good[indexNew] = true; break; }
			}
			centers[indexOld].first = data_[kmerinds[indexNew]].kmer();
		}
	}
	
	// auxiliary variables
	vector<int> dists(l);
	vector<double> loglike(l);
	vector<bool> changedCenter(l);

	// main loop
	bool changed = true;
	bool likelihood_improved = true;
	double total_likelihood = - std::numeric_limits<double>::infinity();
	while (changed && likelihood_improved) {
		// fill everything with zeros
		changed = false;
		fill(changedCenter.begin(), changedCenter.end(), false);
		for (uint32_t j=0; j < l; ++j) centers[j].second.first = 0;
		
		double cur_total = 0;
		// E step: find which clusters we belong to
		for (size_t i=0; i < kmerinds.size(); ++i) {
      const KMerStat &kms = data_[kmerinds[i]];
      const KMer &kmer = kms.kmer();
			for (uint32_t j=0; j < l; ++j) {
				dists[j] = hamdistKMer(kmer, centers[j].first);
				loglike[j] = logLikelihoodKMer(centers[j].first, kms);
			}
			if (cfg::get().bayes_debug_output) {
				#pragma omp critical
				{
				cout << "      likelihoods for " << i << ": ";
				for (size_t j=0; j < l; ++j) {
					cout << loglike[j] << " ";
				}
				cout << endl;
				}
			}
			int newInd = cfg::get().bayes_use_hamming_dist ?
				 (min_element(dists.begin(), dists.end()) - dists.begin()) :
				 (max_element(loglike.begin(), loglike.end()) - loglike.begin());
			cur_total += loglike[newInd];
			if (indices[i] != newInd) {
				changed = true;
				changedCenter[indices[i]] = true;
				changedCenter[newInd] = true;
				indices[i] = newInd;
			}
			++centers[indices[i]].second.first;
		}
		if (cfg::get().bayes_debug_output) {
			#pragma omp critical
			{
			cout << "      total likelihood=" << cur_total << " as compared to previous " << total_likelihood << endl;
			}
		}
		likelihood_improved = (cur_total > total_likelihood);
		if ( likelihood_improved ) total_likelihood = cur_total;
		// M step: find new cluster centers
		for (uint32_t j=0; j < l; ++j) {
			if (!changedCenter[j]) continue; // nothing has changed
			centers[j].first = find_consensus_with_mask(kmerinds, indices, j);
		}
	}

	// last M step
	for (uint32_t j=0; j < l; ++j) {
		centers[j].first = find_consensus_with_mask(kmerinds, indices, j);
	}
	
	curLikelihood = trueClusterLogLikelihood(kmerinds, centers, indices);

	if (restartCount > 1 && curLikelihood > bestLikelihood) {
		bestLikelihood = curLikelihood;
		bestCenters = centers; bestIndices = indices;
	}
	
	} // end restarts
	
	if (restartCount > 1) {
		centers = bestCenters; indices = bestIndices; curLikelihood = bestLikelihood;
	}
	return curLikelihood;
}

size_t KMerClustering::process_block_SIN(const std::vector<unsigned> & block, vector< vector<int> > & vec) {
  size_t newkmers = 0;
  
	if ( cfg::get().bayes_debug_output ) {
		#pragma omp critical
		{
		cout << "process_SIN with block size=" << block.size() << " total kmers=" << data_.size() << endl << "block:";
		for (size_t i = 0; i < block.size(); i++) {
			cout << " " << block[i];
		}
		cout << endl << "  kmers:\n";
		for (size_t i = 0; i < block.size(); i++) {
			cout << data_[block[i]].kmer().str() << endl;
		}
		}
	}

	size_t origBlockSize = block.size();
	if (origBlockSize == 0) return 0;
	
  KMerHamDistMatrix distances(origBlockSize, origBlockSize);
	string newkmer;
	string reason = "noreason";

	// Calculate distance matrix
	for (size_t i = 0; i < block.size(); i++) {
		distances(i, i) = 0;
    const KMerStat &kmsx = data_[block[i]];
    const KMer &kmerx = kmsx.kmer();
		for (size_t j = i + 1; j < block.size(); j++) {
      const KMerStat &kmsy = data_[block[j]];
      const KMer &kmery = kmsy.kmer();
			distances(i, j) = hamdistKMer(kmerx, kmery);
		}
	}

	if (cfg::get().bayes_debug_output) {
		#pragma omp critical
		{
			cout << "\nClustering an interesting block." << endl;
		}
	}

	vector<int> indices(origBlockSize);
	double bestLikelihood = - std::numeric_limits<double>::infinity();
	vector<StringCount> bestCenters;
	vector<int> bestIndices(block.size());

	uint32_t max_l = cfg::get().bayes_hammer_mode ? 1 : origBlockSize;

	for (uint32_t l = 1; l <= max_l; ++l) {
		vector<StringCount> centers(l);
		double curLikelihood = lMeansClustering(l, distances, block, indices, centers);
		if (cfg::get().bayes_debug_output) {
			#pragma omp critical
			{
				cout << "    indices: ";
				for (uint32_t i = 0; i < origBlockSize; i++) cout << indices[i] << " ";
				cout << "\n";
				cout << "  likelihood with " << l << " clusters is " << curLikelihood << endl;
			}
		}
		if (curLikelihood <= bestLikelihood && curLikelihood > -std::numeric_limits<double>::infinity()) {
			break;
		}
		bestLikelihood = curLikelihood;
		bestCenters = centers; bestIndices = indices;
	}


	// find if centers are in clusters
  std::vector<int> centersInCluster(bestCenters.size(), -1);
	for (uint32_t i = 0; i < origBlockSize; i++) {
		int dist = hamdistKMer(data_[block[i]].kmer(), bestCenters[bestIndices[i]].first);
		if (dist == 0) {
			centersInCluster[bestIndices[i]] = i;
		}
	}
	
	bool cons_suspicion = false;
	for (size_t k=0; k<bestCenters.size(); ++k) if (centersInCluster[k] == -1) cons_suspicion = true;
	if (cfg::get().bayes_debug_output) {
#   pragma omp critical
		{
      cout << "Centers: \n";
      for (size_t k=0; k<bestCenters.size(); ++k) {
        cout << "  " << bestCenters[k].first.str() << " " << bestCenters[k].second << " ";
        if (centersInCluster[k] >= 0) cout << block[centersInCluster[k]];
        cout << "\n";
      }
      cout << "The entire block:" << endl;
      for (uint32_t i = 0; i < origBlockSize; i++) {
        const KMerStat &kms = data_[block[i]];      
        cout << "  " << kms.kmer().str() << " (" << kms.count << ", " << kms.totalQual << ") " << block[i] << "\n";
        cout << "  "; for (uint32_t j=0; j<K; ++j) cout << setw(4) << (unsigned)getQual(kms, j) << " "; cout << "\n";
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
			if (centersInCluster[k] >= 0) continue;
			for (size_t s=0; s<bestCenters.size(); ++s) {
				if (s == k || centersInCluster[s] < 0) continue;
				int dist = hamdistKMer(bestCenters[k].first, bestCenters[s].first);
				if ( dist == 0 ) {
					// OK, that's the situation, cluster k should be added to cluster s
					for (uint32_t i = 0; i < origBlockSize; i++) {
						if ( indices[i] == (int)k ) {
							indices[i] = (int)s;
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

	if (cfg::get().bayes_debug_output && origBlockSize > 2) {
		#pragma omp critical
		{
			cout << "\nAfter the check we got centers: \n";
			for (size_t k=0; k<bestCenters.size(); ++k) {
				cout << "  " << bestCenters[k].first.str() << " " << bestCenters[k].second << " ";
				if (centersInCluster[k] >= 0) cout << block[centersInCluster[k]];
				cout << "\n";
			}
			cout << endl;
		}
	}

	for (size_t k=0; k<bestCenters.size(); ++k) {
		if (bestCenters[k].second.first == 0) {
			continue; // superfluous cluster with zero elements
		}
		vector<int> v;
		if (bestCenters[k].second.first == 1) {
			for (uint32_t i = 0; i < origBlockSize; i++) {
				if (indices[i] == (int)k) {
					v.push_back(block[i]);
					break;
				}
			}
		} else { // there are several kmers in this cluster
			for (uint32_t i = 0; i < origBlockSize; i++) {
				if (bestIndices[i] == (int)k) {
					if (centersInCluster[k] == (int)i) {
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

void KMerClustering::process(std::vector<std::vector<unsigned> > classes,
                             boost::shared_ptr<std::ofstream> ofs, boost::shared_ptr<std::ofstream> ofs_bad) {
  size_t newkmers = 0;
  size_t gsingl = 0, tsingl = 0, tcsingl = 0, gcsingl = 0, tcls = 0, gcls = 0;

  std::vector<boost::numeric::ublas::matrix<uint64_t> > errs(nthreads_, boost::numeric::ublas::matrix<double>(4, 4));
# pragma omp parallel for shared(classes, ofs, ofs_bad, errs) num_threads(nthreads_) schedule(dynamic) reduction(+:newkmers, gsingl, tsingl, tcsingl, gcsingl, tcls, gcls)
  for (size_t i = 0; i < classes.size(); ++i) {
    auto cur_class = classes[i];

    // No need for clustering for singletons
    if (cur_class.size() == 1) {
      size_t idx = cur_class[0];
      KMerStat &singl = data_[idx];
      if ((1-singl.totalQual) > cfg::get().bayes_singleton_threshold) {
        singl.status = KMerStat::GoodIter;
        gsingl += 1;
        if (std::ofstream *fs = ofs.get())
          (*fs) << " good singleton: " << idx << "\n  " << singl << '\n';
      } else {
        if (cfg::get().correct_use_threshold && (1-singl.totalQual) > cfg::get().correct_threshold)
          singl.status = KMerStat::GoodIterBad;
        else
          singl.status = KMerStat::Bad;
        if (std::ofstream *fs = ofs_bad.get())
          (*fs) << " bad singleton: " << idx << "\n  " << singl << '\n';
      }
      tsingl += 1;
    } else {
      std::vector<std::vector<int> > blocksInPlace;

      newkmers += process_block_SIN(cur_class, blocksInPlace);

      for (size_t m = 0; m < blocksInPlace.size(); ++m) {
        if (blocksInPlace[m].size() == 0)
          continue;

        if (blocksInPlace[m].size() == 1) {
          size_t idx = blocksInPlace[m][0];
          KMerStat &singl = data_[idx];
          if ((1-singl.totalQual) > cfg::get().bayes_singleton_threshold) {
            singl.status = KMerStat::GoodIter;
            gcsingl += 1;
            if (std::ofstream *fs = ofs.get())
              (*fs) << " good cluster singleton: " << idx << "\n  " << singl << '\n';
          } else {
            if (cfg::get().correct_use_threshold && (1-singl.totalQual) > cfg::get().correct_threshold)
              singl.status = KMerStat::GoodIterBad;
            else
              singl.status = KMerStat::Bad;
            if (std::ofstream *fs = ofs_bad.get())
              (*fs) << " bad cluster singleton: " << idx << "\n  " << singl << '\n';
          }
          tcsingl += 1;
        } else {
          size_t cidx = blocksInPlace[m][0];
          KMerStat &center = data_[cidx];

          // we've got a nontrivial cluster; computing its overall quality
          double cluster_quality = 1;
          for (size_t j=1; j < blocksInPlace[m].size(); ++j) {
            cluster_quality *= data_[blocksInPlace[m][j]].totalQual;
          }
          cluster_quality = 1-cluster_quality;

          // in regular hammer mode, all nonsingletons are good
          if (cfg::get().bayes_hammer_mode ||
              cluster_quality > cfg::get().bayes_nonsingleton_threshold) {
            center.status = KMerStat::GoodIter;
            gcls += 1;
            if (std::ofstream *fs = ofs.get())
              (*fs) << " center of good cluster (" << blocksInPlace[m].size() - 1 << ", " << cluster_quality << ")" << "\n  "
                    << center << '\n';
          } else {
            if (cfg::get().correct_use_threshold && (1-center.totalQual) > cfg::get().correct_threshold)
              center.status = KMerStat::GoodIterBad;
            else
              center.status = KMerStat::Bad;
            if (std::ofstream *fs = ofs_bad.get())
              (*fs) << " center of bad cluster (" << blocksInPlace[m].size() - 1 << ", " << cluster_quality << ")" << "\n  "
                    << center << '\n';
          }
          tcls += 1;

          for (size_t j = 1; j < blocksInPlace[m].size(); ++j) {
            size_t eidx = blocksInPlace[m][j];
            KMerStat &kms = data_[eidx];

            kms.set_change(cidx);
            UpdateErrors(errs[omp_get_thread_num()], kms, center);
            if (std::ofstream *fs = ofs_bad.get())
              (*fs) << " part of cluster (" << blocksInPlace[m].size() - 1 << ", " << cluster_quality << ")" << "\n  "
                    << kms << '\n';
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
  INFO("  Total solid k-mers: " << gsingl + gcsingl + gcls);
  INFO("  Substitution probabilities: " << err);
}

