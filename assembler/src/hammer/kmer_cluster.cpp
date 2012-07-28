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

double KMerClustering::logLikelihoodKMer(const string & center, const KMerCount & x) {
	double res = 0;
	for (unsigned i = 0; i < K; ++i) {
		if (center[i] != x.first[i]) {
			res += - log(10) * getQual(x, i) / 10.0;
		} else {
      res += getProb(x, i, /* log */ true);
    }
  }

	return res;
}

double KMerClustering::logLikelihoodSingleton(const KMerCount & x) {
	double res = 0;
	for (unsigned i = 0; i < K; ++i) {
		res += getProb(x, i, /* log */ true);
	}
	return res;
}

int KMerClustering::hamdistKMer(const PositionKMer & x, const PositionKMer & y, int tau) {
	int dist = 0;
	for (uint32_t i = 0; i < K; ++i) {
		if (x[i] != y[i]) {
			++dist; if (dist > tau) return dist;
		}
	}
	return dist;
}

int KMerClustering::hamdistKMer(const hint_t & x, const hint_t & y, int tau) {
	int dist = 0;
	for (uint32_t i = 0; i < K; ++i) {
		if (Globals::blob[x + i] != Globals::blob[y + i]) {
			++dist; if (dist > tau) return dist;
		}
	}
	return dist;
}

int KMerClustering::hamdistKMer(const string & x, const string & y, int tau) {
	int dist = 0;
	for (uint32_t i = 0; i < K; ++i) {
		if (x[i] != y[i]) {
			++dist; if (dist > tau) return dist;
		}
	}
	return dist;
}

int KMerClustering::hamdistKMer(const PositionKMer & x, const string & y, int tau) {
	int dist = 0;
	for (uint32_t i = 0; i < K; ++i) {
		if ( x[i] != y[i]) {
			++dist; if (dist > tau) return dist;
		}
	}
	return dist;
}

/**
  * SIN
  * find consensus with mask
  * @param mask is a vector of integers of the same size as the block
  * @param maskVal is the integer that we use
  */
string KMerClustering::find_consensus_with_mask(const vector<int> & cl, const vector<int> & mask, int maskVal) {

  //entry_logger el ("find_consensus_with_mask");

  size_t blockSize = cl.size();

  // consensus of a single string is trivial
  if (blockSize == 1) return k_[cl[0]].first.str();

  string c(K, 'A');
  for (uint32_t i = 0; i < K; i++) {
    int scores[4] = {0,0,0,0};
    for (uint32_t j = 0; j < blockSize; j++) {
      if (mask[j] == maskVal)
        scores[static_cast<uint8_t>(dignucl(k_[cl[j]].first[i]))] += k_[cl[j]].second.count;
    }
    c[i] = nucl(max_element(scores, scores + 4) - scores);
  }

  return c;
}


string KMerClustering::find_consensus(const vector<int> & cl) {
  size_t blockSize = cl.size();

  // consensus of a single string is trivial
  if (blockSize == 1) return k_[cl[0]].first.str();

  string c(K, 'A');
  for (size_t i = 0; i < K; i++) {
    int scores[4] = {0,0,0,0};
    for (size_t j = 0; j < blockSize; j++) {
      scores[static_cast<uint8_t>(dignucl(k_[cl[j]].first[i]))]+=k_[cl[j]].second.count;
    }
    c[i] = nucl(max_element(scores, scores + 4) - scores);
  }
  return c;
}

/**
  * @return total log-likelihood of this particular clustering with real quality values
  */
double KMerClustering::trueClusterLogLikelihood(const vector<int> & cl, const vector<StringCount> & centers, const vector<int> & indices) {
	size_t blockSize = cl.size();
	if (blockSize == 0) return 0.0;
	assert(blockSize == indices.size());
	assert(centers.size() > 0);

	// if there is only one center, there are no beta coefficients
	if (centers.size() == 1) {
		double logLikelihood = 0;
		for (size_t i=0; i<blockSize; ++i) {
			logLikelihood += logLikelihoodKMer(centers[indices[i]].first, k_[cl[i]]);
		}
		return ( lMultinomial(cl, k_) + logLikelihood );
	}

	// compute sufficient statistics
	vector<int> count(centers.size(), 0);		// how many kmers in cluster i
	vector<double> totalLogLikelihood(centers.size(), 0);	// total distance from kmers of cluster i to its center
	for (size_t i=0; i<blockSize; ++i) {
		count[indices[i]]+=k_[cl[i]].second.count;
		totalLogLikelihood[indices[i]] += logLikelihoodKMer(centers[indices[i]].first, k_[cl[i]]);
	}

	// sum up the likelihood
	double res = lBetaPlusOne(count);   // 1/B(count)
	res += lMultinomial(centers); 		// {sum(centers.count) \choose centers.count}
	for (size_t i=0; i<centers.size(); ++i) {
		res += lMultinomialWithMask(cl, k_, indices, i) + totalLogLikelihood[i];
	}
	res -= logSimplexVolume(centers.size());
	return res;
}

/**
  * @return total log-likelihood of this particular clustering with real quality values
  */
double KMerClustering::trueSingletonLogLikelihood(const hint_t & kmind) {
	return logLikelihoodSingleton(k_[kmind]);
}


double KMerClustering::lMeansClustering(uint32_t l, vector< vector<int> > & distances, const vector<int> & kmerinds, vector<int> & indices, vector<StringCount> & centers) {
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
			if (distances[i][j] > distances[fp.first][fp.second] ||
			  (distances[i][j] == distances[fp.first][fp.second] &&
					( k_[kmerinds[i]].second.count + k_[kmerinds[j]].second.count >
					  k_[kmerinds[fp.first]].second.count + k_[kmerinds[fp.second]].second.count ||
					  ( k_[kmerinds[i]].second.count + k_[kmerinds[j]].second.count ==
  					    k_[kmerinds[fp.first]].second.count + k_[kmerinds[fp.second]].second.count &&
  					    k_[kmerinds[i]].second.totalQual + k_[kmerinds[j]].second.totalQual >
  					    k_[kmerinds[fp.first]].second.totalQual + k_[kmerinds[fp.second]].second.totalQual)))) {
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
			for (uint32_t k=0; k<j; ++k) curDist += distances[i][init[k]];

			if ( curDist > bestDist ) {
				ind = i; bestDist = curDist;
				continue;
			}
		}
		init[j] = ind;
	}

	for (uint32_t j=0; j<l; ++j) {
		centers[j].first = k_[kmerinds[init[j]]].first.str();
		centers[j].second.first = k_[kmerinds[init[j]]].second.count;
		centers[j].second.second = k_[kmerinds[init[j]]].second.totalQual;
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
				if (k_[kmerinds[i]].second.count >= centers[j].second.first ||
						( k_[kmerinds[i]].second.count == centers[j].second.first &&
						  k_[kmerinds[i]].second.totalQual > centers[j].second.second )) {
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
			centers[indexOld].first = k_[kmerinds[indexNew]].first.str();
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
			for (uint32_t j=0; j < l; ++j) {
				dists[j] = hamdistKMer(k_[kmerinds[i]].first, centers[j].first);
				loglike[j] = logLikelihoodKMer(centers[j].first, k_[kmerinds[i]]);
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

void KMerClustering::process_block_SIN(const vector<int> & block, vector< vector<int> > & vec) {

	if ( cfg::get().bayes_debug_output ) {
		#pragma omp critical
		{
		cout << "process_SIN with block size=" << block.size() << " total kmers=" << k_.size() << endl << "block:";
		for (size_t i = 0; i < block.size(); i++) {
			cout << " " << block[i];
		}
		cout << endl << "  kmers:\n";
		for (size_t i = 0; i < block.size(); i++) {
			cout << k_[block[i]].first.str() << endl;
		}
		}
	}

	uint32_t origBlockSize = block.size();
	if (origBlockSize == 0) return;
	
	vector<double> multiCoef(origBlockSize,1000000);
	vector<int> distance(origBlockSize, 0);
	vector< vector<int> > distances(origBlockSize, distance);
	string newkmer;
	string reason = "noreason";

	//Calculate distance matrix
	for (size_t i = 0; i < block.size(); i++) {
		distances[i][i] = 0;
		for (size_t j = i + 1; j < block.size(); j++) {
			distances[i][j] = hamdistKMer(k_[block[i]].first, k_[block[j]].first);
			distances[j][i] = distances[i][j];
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
	vector<int> centersInCluster(bestCenters.size(), -1);
	for (uint32_t i = 0; i < origBlockSize; i++) {
		int dist = hamdistKMer(k_[block[i]].first, bestCenters[bestIndices[i]].first);
		if (dist == 0) {
			centersInCluster[bestIndices[i]] = i;
		}
	}
	
	bool cons_suspicion = false;
	for (size_t k=0; k<bestCenters.size(); ++k) if (centersInCluster[k] == -1) cons_suspicion = true;
	if (cfg::get().bayes_debug_output) {
		#pragma omp critical
		{
		cout << "Centers: \n";
		for (size_t k=0; k<bestCenters.size(); ++k) {
			cout << "  " << bestCenters[k].first.data() << " " << bestCenters[k].second << " ";
			if ( centersInCluster[k] >= 0 ) cout << k_[block[centersInCluster[k]]].first.start();
			cout << "\n";
		}
		cout << "The entire block:" << endl;
		for (uint32_t i = 0; i < origBlockSize; i++) {
			cout << "  " << k_[block[i]].first.str().data() << " " << k_[block[i]].first.start() << "\n";
			cout << "  " << k_[block[i]].first.strQual().data() << " " << k_[block[i]].second.totalQual << "\n";
			cout << "  "; for (uint32_t j=0; j<K; ++j) cout << (int)getQual(k_[block[i]], j) << " "; cout << "\n";
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
				cout << "  " << bestCenters[k].first.data() << " " << bestCenters[k].second << " ";
				if ( centersInCluster[k] >= 0 ) cout << k_[block[centersInCluster[k]]].first.start();
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
				#pragma omp critical
				{
					// change blob
					for (uint32_t j = 0; j < bestCenters[k].first.size(); ++j) {
						Globals::blob[Globals::blob_size + j] = bestCenters[k].first[j];
					}
					Globals::blob_size += bestCenters[k].first.size();

					// add position read
					PositionRead rs(Globals::blob_size - bestCenters[k].first.size(), bestCenters[k].first.size(), Globals::pr->size() - 1);
					Globals::pr->push_back(rs);

					PositionKMer pkm(Globals::pr->size() -1, 0);
					KMerStat kms(Globals::use_common_quality, 0, KMERSTAT_GOODITER, 1);
					k_.push_back(KMerCount(pkm, kms));
				}
				v.insert(v.begin(), k_.size() - 1);
			}
		}
		vec.push_back(v);
	}
}

void KMerClustering::process(boost::shared_ptr<std::ofstream> ofs, boost::shared_ptr<std::ofstream> ofs_bad) {
  TIMEDLN("Starting subclustering in " << nthreads_ << " threads.");
  TIMEDLN("Estimated: size=" << k_.size() << " mem=" << sizeof(KMerCount)*k_.size() << " clustering buffer size=" << cfg::get().hamming_class_buffer);

  std::string fname = HammerTools::getFilename(cfg::get().input_working_dir, Globals::iteration_no, "kmers.hamming");
  MMappedReader ifs(fname, /* unlink */ true);

	vector< vector< vector<int> > > blocks(nthreads_);

	string buf;

	size_t cur_class_num = 0;
  std::vector< std::vector<int> > curClasses;
  std::vector<int> cur_class;

  std::vector< std::vector< std::vector<int> > > blocksInPlace(cfg::get().hamming_class_buffer);

	while (ifs.good()) {
		curClasses.clear();
#if 0
    curClasses.shrink_to_fit();
#else
    std::vector<std::vector<int> >().swap(curClasses);
#endif

		size_t i_nontriv = 0;
		size_t cur_total_size = 0;
		size_t orig_class_num = cur_class_num;
		while (cur_total_size < (size_t)cfg::get().hamming_class_buffer && ifs.good()) {
			cur_class.clear();
#if 0
      cur_class.shrink_to_fit();
#else
      std::vector<int>().swap(cur_class);
#endif
      size_t classNum, sizeClass;
      ifs.read(&classNum, sizeof(classNum));
      ifs.read(&sizeClass, sizeof(sizeClass));
      cur_class.resize(sizeClass);
      ifs.read(&cur_class[0], sizeClass * sizeof(cur_class[0]));
			++cur_class_num;

			// processing singletons immediately
			if ( cur_class.size() == 1 ) {
				if ( (1-k_[cur_class[0]].second.totalQual) > cfg::get().bayes_singleton_threshold) {
					k_[cur_class[0]].second.changeto = KMERSTAT_GOODITER;
					if (std::ofstream *fs = ofs.get()) {
						(*fs) << k_[cur_class[0]].first.str() << "\n>" << k_[cur_class[0]].first.start()
							<< " good singleton " << "  ind=" << cur_class[0]
							<< "  cnt=" << k_[cur_class[0]].second.count
							<< "  tql=" << (1-k_[cur_class[0]].second.totalQual) << "\n";
					}
				} else {
					if (cfg::get().correct_use_threshold && (1-k_[cur_class[0]].second.totalQual) > cfg::get().correct_threshold)
						k_[cur_class[0]].second.changeto = KMERSTAT_GOODITER_BAD;
					else
						k_[cur_class[0]].second.changeto = KMERSTAT_BAD;
					if (std::ofstream *fs = ofs_bad.get()) {
						(*fs) << k_[cur_class[0]].first.str() << "\n>" << k_[cur_class[0]].first.start()
							<< " bad singleton "
							<< "  ind=" << cur_class[0]
							<< "  cnt=" << k_[cur_class[0]].second.count
							<< "  tql=" << (1-k_[cur_class[0]].second.totalQual) << "\n";
					}
				}
			} else {
				curClasses.push_back(cur_class);
				++i_nontriv;
				cur_total_size += cur_class.size();
			}
		}
		TIMEDLN("Processing " << i_nontriv << " nontrivial clusters from " << orig_class_num << " to " << cur_class_num << " in " << nthreads_ << " threads.");

		VERIFY(blocksInPlace.size() >= i_nontriv && curClasses.size() >= i_nontriv);

		#pragma omp parallel for shared(blocksInPlace, curClasses) num_threads(nthreads_)
		for (size_t i=0; i < i_nontriv; ++i) {
			blocksInPlace[i].clear();
#if 0
      blocksInPlace[i].shrink_to_fit();
#else
      std::vector< std::vector<int> >().swap(blocksInPlace[i]);
#endif
      process_block_SIN(curClasses[i], blocksInPlace[i]);
		}

		for (size_t n=0; n < i_nontriv; ++n) {
			for (uint32_t m = 0; m < blocksInPlace[n].size(); ++m) {
				if (blocksInPlace[n][m].size() == 0) continue;
				if (blocksInPlace[n][m].size() == 1) {
					if ( (1-k_[blocksInPlace[n][m][0]].second.totalQual) > cfg::get().bayes_singleton_threshold) {
						k_[blocksInPlace[n][m][0]].second.changeto = KMERSTAT_GOODITER;
						if (std::ofstream *fs = ofs.get()) {
							(*fs) << k_[blocksInPlace[n][m][0]].first.str() << "\n>" << k_[blocksInPlace[n][m][0]].first.start()
                    << " good singleton " << "  ind=" << blocksInPlace[n][m][0]
                    << "  cnt=" << k_[blocksInPlace[n][m][0]].second.count
                    << "  tql=" << (1-k_[blocksInPlace[n][m][0]].second.totalQual) << "\n";
						}
					} else {
						if (cfg::get().correct_use_threshold && (1-k_[blocksInPlace[n][m][0]].second.totalQual) > cfg::get().correct_threshold)
							k_[blocksInPlace[n][m][0]].second.changeto = KMERSTAT_GOODITER_BAD;
						else
							k_[blocksInPlace[n][m][0]].second.changeto = KMERSTAT_BAD;
						if (std::ofstream *fs = ofs_bad.get()) {
							(*fs) << k_[blocksInPlace[n][m][0]].first.str() << "\n>" << k_[blocksInPlace[n][m][0]].first.start()
                    << " bad singleton " << "  ind=" << blocksInPlace[n][m][0]
                    << "  cnt=" << k_[blocksInPlace[n][m][0]].second.count
                    << "  tql=" << (1-k_[blocksInPlace[n][m][0]].second.totalQual) << "\n";
						}
					}
				} else {
					// we've got a nontrivial cluster; computing its overall quality
					double cluster_quality = 1;
					for (uint32_t j=1; j < blocksInPlace[n][m].size(); ++j) {
						cluster_quality *= k_[blocksInPlace[n][m][j]].second.totalQual;
					}
					cluster_quality = 1-cluster_quality;

					// in regular hammer mode, all nonsingletons are good
					if (cfg::get().bayes_hammer_mode || cluster_quality > cfg::get().bayes_nonsingleton_threshold) {
						k_[blocksInPlace[n][m][0]].second.changeto = KMERSTAT_GOODITER;
						if (std::ofstream *fs = ofs.get()) {
							(*fs) << k_[blocksInPlace[n][m][0]].first.str() << "\n>" << k_[blocksInPlace[n][m][0]].first.start()
                    << " center clust=" << cluster_quality
                    << " ind=" << blocksInPlace[n][m][0]
                    << " cnt=" << k_[blocksInPlace[n][m][0]].second.count
                    << " tql=" << (1-k_[blocksInPlace[n][m][0]].second.totalQual) << "\n";
						}
					} else {
						if (cfg::get().correct_use_threshold && (1-k_[blocksInPlace[n][m][0]].second.totalQual) > cfg::get().correct_threshold)
							k_[blocksInPlace[n][m][0]].second.changeto = KMERSTAT_GOODITER_BAD;
						else
							k_[blocksInPlace[n][m][0]].second.changeto = KMERSTAT_BAD;
						if (std::ofstream *fs = ofs_bad.get()) {
							(*fs) << k_[blocksInPlace[n][m][0]].first.str() << "\n>" << k_[blocksInPlace[n][m][0]].first.start()
                    << " center of bad cluster clust=" << cluster_quality
                    << " ind=" << blocksInPlace[n][m][0]
                    << " cnt=" << k_[blocksInPlace[n][m][0]].second.count
                    << " tql=" << (1-k_[blocksInPlace[n][m][0]].second.totalQual) << "\n";
						}
					}
					for (uint32_t j=1; j < blocksInPlace[n][m].size(); ++j) {
						k_[blocksInPlace[n][m][j]].second.changeto = blocksInPlace[n][m][0];
						if (std::ofstream *fs = ofs_bad.get()) {
							(*fs) << k_[blocksInPlace[n][m][j]].first.str() << "\n>" << k_[blocksInPlace[n][m][j]].first.start()
                    << " part of cluster " << k_[blocksInPlace[n][m][0]].first.start() << " clust=" << cluster_quality
                    << " ind=" << blocksInPlace[n][m][j]
                    << " cnt=" << k_[blocksInPlace[n][m][j]].second.count
                    << " tql=" << (1-k_[blocksInPlace[n][m][j]].second.totalQual) << "\n";
						}
					}
				}
			}
		}
	}
}

