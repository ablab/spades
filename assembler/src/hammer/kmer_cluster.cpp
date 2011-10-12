/*
 * kmer_cluster.cpp
 *
 *  Created on: 16.07.2011
 *      Author: snikolenko
 */

#include <omp.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include "read/ireadstream.hpp"
#include "defs.hpp"
#include "mathfunctions.hpp"
#include "hammer_tools.hpp"
#include "kmer_cluster.hpp"

using std::max_element;
using std::min_element;

double KMerClustering::logLikelihoodKMer(const string & center, const KMerCount * x) {
	double res = 0;
	bool change = false;
	for (uint32_t i = 0; i < K; ++i) {
		if (center[i] != x->first[i]) {
			change = true;
			res += - log(10) * x->second.qual[i] / 10.0;
		} else {
			res += log( 1 - pow( 10, -x->second.qual[i] / 10.0 ) );
		}
	}
	if (change) {
		cout << "   logLikelihood btw\n";
		cout << "  " << center.data() << "\n  " <<  x->first.str().data() << "\n  " << x->first.strQual() << "\n  ";
		for (uint32_t i=0; i<K; ++i) cout << x->second.qual[i] << "  ";
		cout << endl;
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

void KMerClustering::processBlock(unionFindClass * uf, vector<hint_t> & block, int cur_subkmer) {
	uint32_t blockSize = block.size();

	if (blockSize < (uint32_t)Globals::blocksize_quadratic_threshold) {
		processBlockQuadratic(uf, block);
	} else {
		int nthreads_per_subkmer = max( (int)(nthreads_ / (tau_ + 1)), 1);
		SubKMerSorter subsubsorter( &block, &k_, nthreads_per_subkmer, tau_, cur_subkmer,
			SubKMerSorter::SorterTypeChequered, SubKMerSorter::SorterTypeStraight );
		subsubsorter.runSort();
		for (int sub_i = 0; sub_i < tau_+1; ++sub_i) {
			vector<hint_t> subblock;
			while ( subsubsorter.getNextBlock(sub_i, subblock) ) {
				processBlockQuadratic(uf, subblock);
			}
		}		
	}
}

void KMerClustering::processBlockQuadratic(unionFindClass * uf, vector<hint_t> & block) {
	uint32_t blockSize = block.size();

	for (uint32_t i = 0; i < blockSize; ++i) {
		uf->find_set(block[i]);
		for (uint32_t j = i + 1; j < blockSize; j++) {
			if (hamdistKMer(k_[block[i]]->first, k_[block[j]]->first, tau_ ) <= tau_) {
				uf->unionn(block[i], block[j]);
			}
		}
	}
	return;
}


void KMerClustering::clusterMerge(vector<unionFindClass *>uf, unionFindClass * ufMaster) {
	vector<string> row;
	vector<vector<int> > classes;
	for (uint32_t i = 0; i < uf.size(); i++) {
		classes.clear();
		uf[i]->get_classes(classes);
		delete uf[i];
		for (uint32_t j = 0; j < classes.size(); j++) {
			uint32_t first = classes[j][0];
			ufMaster->find_set(first);
			for (uint32_t k = 0; k < classes[j].size(); k++) {
				ufMaster->unionn(first, classes[j][k]);
			}
		}
	}
}

double KMerClustering::calcMultCoef(vector<int> & distances, const vector<int> & cl) {
	double prob = 0;
	double theta;
	for (size_t i = 0; i < cl.size(); i++) {
		theta = k_[cl[i]]->second.count * -((K - distances[i]) * log(1 - Globals::error_rate) + distances[i] * log(Globals::error_rate));
		prob = prob + theta;
	}
	return prob;
}


/**
  * SIN
  * find consensus with mask
  * @param mask is a vector of integers of the same size as the block
  * @param maskVal is the integer that we use
  */
string KMerClustering::find_consensus_with_mask(const vector<int> & cl, const vector<int> & mask, int maskVal) {
	size_t blockSize = cl.size();

	// consensus of a single string is trivial
	if (blockSize == 1) return k_[cl[0]]->first.str();

	string c(K, 'A');
	for (uint32_t i = 0; i < K; i++) {
		int scores[4] = {0,0,0,0};
		for (uint32_t j = 0; j < blockSize; j++) {
			if (mask[j] == maskVal)
                          scores[static_cast<uint8_t>(dignucl(k_[cl[j]]->first[i]))] += k_[cl[j]]->second.count;
		}
		c[i] = num2nt(max_element(scores, scores + 4) - scores);
	}
	return c;
}


string KMerClustering::find_consensus(const vector<int> & cl) {
	size_t blockSize = cl.size();

	// consensus of a single string is trivial
	if (blockSize == 1) return k_[cl[0]]->first.str();

	string c(K, 'A');
	for (size_t i = 0; i < K; i++) {
		int scores[4] = {0,0,0,0};
		for (size_t j = 0; j < blockSize; j++) {
                  scores[static_cast<uint8_t>(dignucl(k_[cl[j]]->first[i]))]+=k_[cl[j]]->second.count;
		}
		c[i] = num2nt(max_element(scores, scores + 4) - scores);
	}
	return c;
}


/**
  * @return total log-likelihood of this particular clustering
  */
double KMerClustering::clusterLogLikelihood(const vector<int> & cl, const vector<StringCount> & centers, const vector<int> & indices) {
	size_t blockSize = cl.size();
	if (blockSize == 0) return 0.0;
	assert(blockSize == indices.size());
	assert(centers.size() > 0);
	
	// if there is only one center, there are no beta coefficients
	if (centers.size() == 1) {
		int dist = 0;
		for (size_t i=0; i<blockSize; ++i) {
			dist += hamdistKMer(k_[cl[i]]->first, centers[indices[i]].first);
		}
		return ( lMultinomial(cl, k_) + log(Globals::error_rate) * dist + log(1-Globals::error_rate) * (K * blockSize - dist) );
	}
	
	// compute sufficient statistics
	vector<int> count(centers.size(), 0);		// how many kmers in cluster i
	vector<int> totaldist(centers.size(), 0);	// total distance from kmers of cluster i to its center
	for (size_t i=0; i<blockSize; ++i) {
		count[indices[i]]+=k_[cl[i]]->second.count;
		totaldist[indices[i]] += hamdistKMer(k_[cl[i]]->first, centers[indices[i]].first);
	}
	
	// sum up the likelihood
	double res = lBetaPlusOne(count);   // 1/B(count)
	res += lMultinomial(centers); 		// {sum(centers.count) \choose centers.count}
	for (size_t i=0; i<centers.size(); ++i) {
		res += lMultinomialWithMask(cl, k_, indices, i) + 
			   log(Globals::error_rate) * totaldist[i] + log(1-Globals::error_rate) * (K * count[i] - totaldist[i]);
	}
	return res;
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

	cout << "\nCluster of " << centers.size() << " subclusters\n";

	// compute sufficient statistics
	vector<int> count(centers.size(), 0);		// how many kmers in cluster i
	vector<double> totalLogLikelihood(centers.size(), 0);	// total distance from kmers of cluster i to its center
	for (size_t i=0; i<blockSize; ++i) {
		count[indices[i]]+=k_[cl[i]]->second.count;
		totalLogLikelihood[indices[i]] += logLikelihoodKMer(centers[indices[i]].first, k_[cl[i]]);
	}

	cout << "  counts: ";
	for (size_t i=0; i<centers.size(); ++i) {
		cout <<  count[i] << " ";
	}
	cout << "\n  loglikelihoods: ";
	for (size_t i=0; i<centers.size(); ++i) {
			cout <<  totalLogLikelihood[i] << " ";
	}
	cout << endl;

	// sum up the likelihood
	double res = lBetaPlusOne(count);   // 1/B(count)
	res += lMultinomial(centers); 		// {sum(centers.count) \choose centers.count}
	for (size_t i=0; i<centers.size(); ++i) {
		res += lMultinomialWithMask(cl, k_, indices, i) + totalLogLikelihood[i];
	}
	return res;
}


double KMerClustering::lMeansClustering(int l, vector< vector<int> > & distances, const vector<int> & kmerinds, vector<int> & indices, vector<StringCount> & centers) {
	centers.resize(l); // there are l centers

	// if l==1 then clustering is trivial
	if (l == 1) {
		centers[0].first = find_consensus(kmerinds);
		centers[0].second = kmerinds.size();
		for (size_t i=0; i < kmerinds.size(); ++i) indices[i] = 0;
		return ( Globals::use_true_likelihood ? trueClusterLogLikelihood(kmerinds, centers, indices) : clusterLogLikelihood(kmerinds, centers, indices) );
	}

	int restartCount = 1;
	
	double bestLikelihood = -1000000.0;
	double curLikelihood = -1000000.0;
	vector<StringCount> bestCenters;
	vector<int> bestIndices(kmerinds.size());
	
	for (int restart = 0; restart < restartCount; ++restart) {

	// fill in initial approximations
	for (int j=0; j < l; ++j) centers[j].second = 0;
	for (size_t i=0; i < kmerinds.size(); ++i) {
		for (int j=0; j < l; ++j) {
			if (k_[kmerinds[i]]->second.count > centers[j].second) {
				for (int s=j; s<l-1; ++s) {
					centers[s+1].first = centers[s].first;
					centers[s+1].second = centers[s].second;
				}
				centers[j].first = k_[kmerinds[i]]->first.str();
				centers[j].second = k_[kmerinds[i]]->second.count;
				break;
			}
		}
	}

	// TODO: make random restarts better!!! they don't really work now
	if (restart > 0) { // introduce random noise
		vector<bool> good(kmerinds.size(), false);
		for (size_t i=0; i < kmerinds.size(); ++i) {
			for (int j=0; j < l; ++j) {
				if (k_[kmerinds[i]]->second.count >= centers[j].second) {
					good[i] = true; break;
				}
			}
		}
		for (int k=0; k < l/2 + 1; ++k) {
			int newNo = rand() % (kmerinds.size() - l - k);
			int indexOld = (rand() % (l-1)) + 1;
			int indexNew = kmerinds.size()-1;
			for (size_t i=0; i < kmerinds.size(); ++i) {
				if (!good[i]) --newNo;
				if (newNo < 0) { indexNew = i; good[indexNew] = true; break; }
			}
			centers[indexOld].first = k_[kmerinds[indexNew]]->first.str();
		}
	}
	
	// auxiliary variables
	vector<int> dists(l);
	vector<bool> changedCenter(l);

	// main loop
	bool changed = true;
	while (changed) {
		// fill everything with zeros
		changed = false;
		fill(changedCenter.begin(), changedCenter.end(), false);
		for (int j=0; j < l; ++j) centers[j].second = 0;
		
		// E step: find which clusters we belong to
		for (size_t i=0; i < kmerinds.size(); ++i) {
			for (int j=0; j < l; ++j) {
				dists[j] = hamdistKMer(k_[kmerinds[i]]->first, centers[j].first);
			}
			int newInd = min_element(dists.begin(), dists.end()) - dists.begin();
			if (indices[i] != newInd) {
				changed = true;
				changedCenter[indices[i]] = true;
				changedCenter[newInd] = true;
				indices[i] = newInd;
			}
			++centers[indices[i]].second;
		}
		// M step: find new cluster centers
		for (int j=0; j < l; ++j) {
			if (!changedCenter[j]) continue; // nothing has changed
			centers[j].first = find_consensus_with_mask(kmerinds, indices, j);
		}
	}

	// last M step
	for (int j=0; j < l; ++j) {
		centers[j].first = find_consensus_with_mask(kmerinds, indices, j);
	}
	
	curLikelihood = ( Globals::use_true_likelihood ? trueClusterLogLikelihood(kmerinds, centers, indices) : clusterLogLikelihood(kmerinds, centers, indices) );

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
	int origBlockSize = block.size();
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
			distances[i][j] = hamdistKMer(k_[block[i]]->first, k_[block[j]]->first);
			distances[j][i] = distances[i][j];
		}
	}

	// Multinomial coefficients -- why? TODO: remove
	for (int i = 0; i < origBlockSize; i++) multiCoef[i] = calcMultCoef(distances[i], block);

	vector<int> indices(origBlockSize);
	double bestLikelihood = -1000000;
	vector<StringCount> bestCenters;
	vector<int> bestIndices(block.size());

	for (int l = 1; l <= origBlockSize; ++l) {
		vector<StringCount> centers(l);
		double curLikelihood = lMeansClustering(l, distances, block, indices, centers);
		if (curLikelihood <= bestLikelihood) {
			break;
		}
		bestLikelihood = curLikelihood;
		bestCenters = centers; bestIndices = indices;
	}


	// find if centers are in clusters
	vector<int> centersInCluster(bestCenters.size(), -1);
	for (int i = 0; i < origBlockSize; i++) {
		int dist = hamdistKMer(k_[block[i]]->first, bestCenters[bestIndices[i]].first);
		if (dist == 0) {
			centersInCluster[bestIndices[i]] = i;
		}
	}
	
	bool cons_suspicion = false;
	for (size_t k=0; k<bestCenters.size(); ++k) if (centersInCluster[k] == -1) cons_suspicion = true;
	/*if (cons_suspicion) {
		#pragma omp critical
		{
		cout << "\nCenters: \n";
		for (size_t k=0; k<bestCenters.size(); ++k) {
			cout << "  " << bestCenters[k].first.data() << " " << bestCenters[k].second << " ";
			if ( centersInCluster[k] >= 0 ) cout << k_[block[centersInCluster[k]]].first.start();
			cout << "\n";
		}
		cout << "The entire block:\n";
		for (int i = 0; i < origBlockSize; i++) {
			cout << "  " << k_[block[i]].first.str().data() << " " << k_[block[i]].first.start() << "\n";
		}
		cout << endl;
		}
	//}*/
	
	// it may happen that consensus string from one subcluster occurs in other subclusters
	// we need to check for that
	for (size_t k=0; k<bestCenters.size(); ++k) {
		if (bestCenters[k].second == 0) continue;
		if (centersInCluster[k] >= 0) continue;
		for (size_t s=0; s<bestCenters.size(); ++s) {
			if (s == k || centersInCluster[s] < 0) continue;
			int dist = hamdistKMer(bestCenters[k].first, bestCenters[s].first);
			if ( dist == 0 ) {
				// OK, that's the situation, cluster k should be added to cluster s
				for (int i = 0; i < origBlockSize; i++) {
					if ( indices[i] == (int)k ) {
						indices[i] = (int)s;
						bestCenters[s].second++;
					}
				}
				bestCenters[k].second = 0; // it will be skipped now
				break;
			}
		}
	}

	/*#pragma omp critical
	{
	cout << "\nAfter the check we got centers: \n";
	for (size_t k=0; k<bestCenters.size(); ++k) {
		cout << "  " << bestCenters[k].first.data() << " " << bestCenters[k].second << " ";
		if ( centersInCluster[k] >= 0 ) cout << k_[block[centersInCluster[k]]].first.start();
		cout << "\n";
	}
	cout << "\n";
	}*/

	for (size_t k=0; k<bestCenters.size(); ++k) {
		if (bestCenters[k].second == 0) {
			continue; // superfluous cluster with zero elements
		}
		vector<int> v;
		if (bestCenters[k].second == 1) {
			for (int i = 0; i < origBlockSize; i++) {
				if (indices[i] == (int)k) {
					v.push_back(block[i]);
					break;
				}
			}
		} else { // there are several kmers in this cluster
			for (int i = 0; i < origBlockSize; i++) {
				if (bestIndices[i] == (int)k) {
					if (centersInCluster[k] == i) {
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

					// add read
					Read r("Consensus", bestCenters[k].first, bestCenters[k].first);
					Globals::rv->push_back(r);

					// add position read
					PositionRead rs(Globals::blob_size - r.size(), r.size(), Globals::rv->size() - 1);
					Globals::pr->push_back(rs);

					PositionKMer pkm(Globals::pr->size()-1, 0);
					KMerStat kms( 0, KMERSTAT_GOODITER, 1 );
					k_.push_back( new KMerCount( pkm, kms ) );
				}
				v.insert(v.begin(), k_.size() - 1);

			}
		}
		vec.push_back(v);
	}
}

void KMerClustering::process(string dirprefix, SubKMerSorter * skmsorter, ofstream * ofs, ofstream * ofs_bad) {
	
	int effective_threads = min(nthreads_, tau_+1);
	vector<unionFindClass *> uf(tau_ + 1);
	
	#pragma omp parallel for shared(uf, skmsorter) num_threads(effective_threads)
	for (int i = 0; i < tau_ + 1; i++) {
		uf[i] = new unionFindClass(k_.size()); 

		vector<hint_t> block;
		while ( skmsorter->getNextBlock(i, block) ) {
			processBlock(uf[i], block, i);
		}

	}
	TIMEDLN("All split kmer threads finished. Starting merge.");
	
	unionFindClass * ufMaster;
	ufMaster = new unionFindClass(k_.size());
	clusterMerge(uf, ufMaster);
	TIMEDLN("Merging finished. Centering begins.");
	vector<vector<int> > classes;
	ufMaster->get_classes(classes);
	delete ufMaster;

	vector< vector< vector<int> > > blocks(nthreads_);

	vector< vector< vector<int> > > blocksInPlace(nthreads_);

	#pragma omp parallel for shared(blocksInPlace, classes) num_threads(nthreads_)
	for (size_t i=0; i < classes.size(); ++i) {
		int n = omp_get_thread_num();
		blocksInPlace[n].clear();
		process_block_SIN(classes[i], blocksInPlace[n]);
		for (uint32_t m = 0; m < blocksInPlace[n].size(); ++m) {
			if (blocksInPlace[n][m].size() == 0) continue;
			if (blocksInPlace[n][m].size() == 1) {
				if ( (1-k_[blocksInPlace[n][m][0]]->second.totalQual) > Globals::good_cluster_threshold) {
					k_[blocksInPlace[n][m][0]]->second.changeto = KMERSTAT_GOOD;
					if ( (1-k_[blocksInPlace[n][m][0]]->second.totalQual) > Globals::iterative_reconstruction_threshold) {
						k_[blocksInPlace[n][m][0]]->second.changeto = KMERSTAT_GOODITER;
					} else {
						k_[blocksInPlace[n][m][0]]->second.changeto = KMERSTAT_GOODITER_BAD;
					}
					#pragma omp critical
					{
					(*ofs) << k_[blocksInPlace[n][m][0]]->first.str() << "\n>" << k_[blocksInPlace[n][m][0]]->first.start()
					       << " good singleton "
					       << "  ind=" << blocksInPlace[n][m][0]
					       << "  cnt=" << k_[blocksInPlace[n][m][0]]->second.count 
					       << "  tql=" << (1-k_[blocksInPlace[n][m][0]]->second.totalQual) << "\n";
					}
				} else {
					k_[blocksInPlace[n][m][0]]->second.changeto = KMERSTAT_BAD;
					#pragma omp critical
					{
					(*ofs_bad) << k_[blocksInPlace[n][m][0]]->first.str() << "\n>" << k_[blocksInPlace[n][m][0]]->first.start()
						   << " bad singleton "
					       	   << "  ind=" << blocksInPlace[n][m][0]
						   << "  cnt=" << k_[blocksInPlace[n][m][0]]->second.count 
						   << "  tql=" << (1-k_[blocksInPlace[n][m][0]]->second.totalQual) << "\n";
					}
				}
			} else {
				// we've got a nontrivial cluster; computing its overall quality
				double cluster_quality = 1;
				for (uint32_t j=1; j < blocksInPlace[n][m].size(); ++j) {
					cluster_quality *= k_[blocksInPlace[n][m][j]]->second.totalQual;
				}
				cluster_quality = 1-cluster_quality;

				if ( cluster_quality > Globals::good_cluster_threshold ||
				     ( Globals::regular_threshold_for_correction && (cluster_quality > Globals::special_nonsingleton_threshold) ) ) {
					k_[blocksInPlace[n][m][0]]->second.changeto = KMERSTAT_GOOD;
					if ( cluster_quality > Globals::iterative_reconstruction_threshold ||
					    (Globals::regular_threshold_for_correction && (cluster_quality > Globals::special_nonsingleton_threshold) )) {
						k_[blocksInPlace[n][m][0]]->second.changeto = KMERSTAT_GOODITER;
					} else {
						k_[blocksInPlace[n][m][0]]->second.changeto = KMERSTAT_GOODITER_BAD;
					}
					#pragma omp critical
					{
					(*ofs) << k_[blocksInPlace[n][m][0]]->first.str() << "\n>" << k_[blocksInPlace[n][m][0]]->first.start()
					       << " center clust=" << cluster_quality
					       << " ind=" << blocksInPlace[n][m][0]
					       << " cnt=" << k_[blocksInPlace[n][m][0]]->second.count
					       << " tql=" << (1-k_[blocksInPlace[n][m][0]]->second.totalQual) << "\n";
					}
				} else {
					k_[blocksInPlace[n][m][0]]->second.changeto = KMERSTAT_BAD;
					#pragma omp critical
					{
					(*ofs_bad) << k_[blocksInPlace[n][m][0]]->first.str() << "\n>" << k_[blocksInPlace[n][m][0]]->first.start()
						   << " center of bad cluster clust=" << cluster_quality
					  	   << " ind=" << blocksInPlace[n][m][0]
						   << " cnt=" << k_[blocksInPlace[n][m][0]]->second.count
						   << " tql=" << (1-k_[blocksInPlace[n][m][0]]->second.totalQual) << "\n";
					}
				}
				for (uint32_t j=1; j < blocksInPlace[n][m].size(); ++j) {
					k_[blocksInPlace[n][m][j]]->second.changeto = blocksInPlace[n][m][0];
					#pragma omp critical
					{
					(*ofs_bad) << k_[blocksInPlace[n][m][j]]->first.str() << "\n>" << k_[blocksInPlace[n][m][j]]->first.start() 
						   << " part of cluster " << k_[blocksInPlace[n][m][0]]->first.start() << " clust=" << cluster_quality
					           << " ind=" << blocksInPlace[n][m][j]
						   << " cnt=" << k_[blocksInPlace[n][m][j]]->second.count
						   << " tql=" << (1-k_[blocksInPlace[n][m][j]]->second.totalQual) << "\n";
					}
				}
			}
		}
	}
	TIMEDLN("Centering finished.");	
}

