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
#include "read/ireadstream.hpp"
#include "defs.hpp"
#include "mathfunctions.hpp"
#include "hammer_tools.hpp"
#include "kmer_cluster.hpp"

using std::max_element;
using std::min_element;

int KMerClustering::hamdistKMer(const PositionKMer & x, const PositionKMer & y, int tau) {
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

void KMerClustering::processBlock(unionFindClass * uf, vector<uint64_t> & block) {
	for (uint32_t i = 0; i < block.size(); i++) {
		uf->find_set(block[i]);
		for (uint32_t j = i + 1; j < block.size(); j++) {
			//cout << "Comparing " << block[i].kmer.str() << "\n          " << block[j].kmer.str() << "\n";
			if (hamdistKMer(k_[block[i]].first, k_[block[j]].first, tau_ ) <= tau_) {
				//cout << "    ok!\n"; 
				uf->unionn(block[i], block[j]);
			}
		}
	}
	return;
}

void KMerClustering::clusterMerge(vector<unionFindClass *>uf, unionFindClass * ufMaster) {
	cout << "Merging union find files..." << endl;
	vector<string> row;
	vector<vector<int> > classes;
	for (uint32_t i = 0; i < uf.size(); i++) {
		classes.clear();
		uf[i]->get_classes(classes);
		delete uf[i];
		/*cout << "classes[" << i << "]: ";
		for (uint32_t j = 0; j < classes.size(); j++) {
			for (uint32_t k = 0; k < classes[j].size(); k++) cout << classes[j][k] << " ";
			cout << "|";
		}
		cout << endl;*/
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
		theta = k_[cl[i]].second.count * -((K - distances[i]) * log(1 - ERROR_RATE) + distances[i] * log(ERROR_RATE));
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
	string c;
	for (uint32_t i = 0; i < K; i++) {
		int scores[4] = {0,0,0,0};
		for (uint32_t j = 0; j < cl.size(); j++) {
			if (mask[j] == maskVal)
                          scores[static_cast<uint8_t>(dignucl(k_[cl[j]].first[i]))] += k_[cl[j]].second.count;
		}
		c.push_back(num2nt(max_element(scores, scores + 4) - scores));
	}
	return c;
}


string KMerClustering::find_consensus(const vector<int> & cl) {
	string c;
	for (size_t i = 0; i < K; i++) {
		int scores[4] = {0,0,0,0};
		for (size_t j = 0; j < cl.size(); j++) {
                  scores[static_cast<uint8_t>(dignucl(k_[cl[j]].first[i]))]+=k_[cl[j]].second.count;
		}
		c.push_back(num2nt(max_element(scores, scores + 4) - scores));
	}
	return c;
}


/**
  * @return total log-likelihood of this particular clustering
  */
double KMerClustering::clusterLogLikelihood(const vector<int> & cl, const vector<StringCount> & centers, const vector<int> & indices) {
	if (cl.size() == 0) return 0.0;
	assert(cl.size() == indices.size());
	assert(centers.size() > 0);
	
	// if there is only one center, there are no beta coefficients
	if (centers.size() == 1) {
		int dist = 0;
		for (size_t i=0; i<cl.size(); ++i) {
			dist += hamdistKMer(k_[cl[i]].first, centers[indices[i]].first);
		}
		return ( lMultinomial(cl, k_) + log(ERROR_RATE) * dist + log(1-ERROR_RATE) * (K * cl.size() - dist) );
	}
	
	// compute sufficient statistics
	vector<int> count(centers.size(), 0);		// how many kmers in cluster i
	vector<int> totaldist(centers.size(), 0);	// total distance from kmers of cluster i to its center
	for (size_t i=0; i<cl.size(); ++i) {
		count[indices[i]]+=k_[cl[i]].second.count;
		totaldist[indices[i]] += hamdistKMer(k_[cl[i]].first, centers[indices[i]].first);
	}
	
	// sum up the likelihood
	double res = lBetaPlusOne(count);   // 1/B(count)
	res += lMultinomial(centers); 		// {sum(centers.count) \choose centers.count}
	for (size_t i=0; i<centers.size(); ++i) {
		res += lMultinomialWithMask(cl, k_, indices, i) + 
			   log(ERROR_RATE) * totaldist[i] + log(1-ERROR_RATE) * (K * count[i] - totaldist[i]);
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
		return clusterLogLikelihood(kmerinds, centers, indices);
	}
	
	/*cout << l << "  " << kmerinds.size() << ": "; 
	for (size_t i=0; i < kmerinds.size(); ++i) cout << kmerinds[i] << " ";
	cout << endl;*/


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
			if (k_[kmerinds[i]].second.count > centers[j].second) {
				for (int s=j; s<l-1; ++s) {
					centers[s+1].first = centers[s].first;
					centers[s+1].second = centers[s].second;
				}
				centers[j].first = k_[kmerinds[i]].first.str();
				centers[j].second = k_[kmerinds[i]].second.count;
				break;
			}
		}
	}

	// TODO: make random restarts better!!! they don't really work now
	if (restart > 0) { // introduce random noise
		vector<bool> good(kmerinds.size(), false);
		for (size_t i=0; i < kmerinds.size(); ++i) {
			for (int j=0; j < l; ++j) {
				if (k_[kmerinds[i]].second.count >= centers[j].second) {
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
			centers[indexOld].first = k_[kmerinds[indexNew]].first.str();
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
				dists[j] = hamdistKMer(k_[kmerinds[i]].first, centers[j].first);
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
	
	curLikelihood = clusterLogLikelihood(kmerinds, centers, indices);

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
			distances[i][j] = hamdistKMer(k_[block[i]].first, k_[block[j]].first);
			distances[j][i] = distances[i][j];
		}
		//cout << k_[block[i]].first.str() << "\t" << k_[block[i]].second.count << endl;
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
		//cout << "   ...likelihood with " << l << " clusters = " << curLikelihood << endl;
		if (curLikelihood <= bestLikelihood) {
			break;
		}
		bestLikelihood = curLikelihood;
		bestCenters = centers; bestIndices = indices;
	}
	
	for (size_t k=0; k<bestCenters.size(); ++k) {
		if (bestCenters[k].second == 0) {
			continue; // superfluous cluster with zero elements
		}
		//cout << "subcluster with " << bestCenters[k].second << " elements" << endl;
		vector<int> v;
		if (bestCenters[k].second == 1) {
			for (int i = 0; i < origBlockSize; i++) {
				if (indices[i] == (int)k) {
					v.push_back(block[i]);
					break;
				}
			}
		} else { // there are several kmers in this cluster
			// cout << "center: " << bestCenters[k].first.str().data() << "\n";
			bool centerInCluster = false;
			for (int i = 0; i < origBlockSize; i++) {
				if (bestIndices[i] == (int)k) {
					int dist = hamdistKMer(k_[block[i]].first, bestCenters[k].first);
					// cout << "        " << k_[block[i]].first.str().data() << " " << dist << "\n";
					if (dist == 0) {
						// cout << "  found center\n";
						centerInCluster = true;
						v.insert(v.begin(), block[i]);
					} else {
						v.push_back(block[i]);
					}
				}
			}
			if (!centerInCluster) {
				// cout << "  pushing consensus\n";
				// cout << "Consensus: " << bestCenters[k].first;

				Read consensus_read("Consensus", bestCenters[k].first, bestCenters[k].first);
				
				#pragma omp critical
				{
					ReadStat rs; rs.read = consensus_read; rs.blobpos = PositionKMer::blob_size;
					for (uint32_t j = 0; j < rs.read.size(); ++j) {
						PositionKMer::blob[PositionKMer::blob_size + j] = rs.read.getSequenceString()[j];
					}
					PositionKMer::blob_size += rs.read.size();
					PositionKMer::rv->push_back(rs);

					PositionKMer pkm(PositionKMer::rv->size()-1, 0);
					KMerStat kms; kms.count = 0; kms.change = false; kms.good = true; kms.changeto = 0;
					k_.push_back( make_pair( pkm, kms ) );
					// cout << "  " << (PositionKMer::rv->size() - 1) << " " << (k_.size()-1) << " " << PositionKMer::blob_size << "     " << pkm.str() << endl;
				}
				v.insert(v.begin(), k_.size() - 1);

			}
		}
		vec.push_back(v);
	}
}

void KMerClustering::process(string dirprefix, const vector< vector<uint64_t> > & vs) {
	
	int effective_threads = min(nthreads_, tau_+1);
	vector<unionFindClass *> uf(tau_ + 1);

	#pragma omp parallel for shared(uf) num_threads(effective_threads)
	for (int i = 0; i < tau_ + 1; i++) {
		uf[i] = new unionFindClass(k_.size()); 

		cout << "Processing split kmers " << i << ", total " << vs[i].size() << "\n";
		
		string sbuf;
		uint64_t last = vs[i][0];
		vector<uint64_t> block;
		for (size_t j=0; j<vs[i].size(); ++j) {
			if (j % 10000000 == 0) cout << "Processed (" << i << ") " << j << ", ";
			if ( PositionKMer::equalSubKMers(last, vs[i][j], &k_, tau_, i) ) { //add to current reads
				block.push_back(vs[i][j]);
			} else {
				processBlock(uf[i], block);
				block.clear();
				block.push_back(vs[i][j]);
				last = vs[i][j];
			}
		}
		processBlock(uf[i], block);
		cout << "Finished(" << i << ") " << endl;
	}
	cout << "All threads finished.\n";flush(cout);
	
	unionFindClass * ufMaster;
	ufMaster = new unionFindClass(k_.size());
	clusterMerge(uf, ufMaster);
	cout << "Merging finished.\n"; flush(cout);
	
	
	cout << "Centering begins...\n"; flush(cout);
	vector<vector<int> > classes;
	ufMaster->get_classes(classes);
	delete ufMaster;
	vector< vector< vector<int> > > blocks(nthreads_);

	/*	cout << "classes: ";
		for (uint32_t j = 0; j < classes.size(); j++) {
			for (uint32_t k = 0; k < classes[j].size(); k++) cout << classes[j][k] << " " << k_[classes[j][k]].first.str() << endl;
			cout << "------" << endl;
		}
		cout << endl;*/

	
	#pragma omp parallel for shared(blocks, classes) num_threads(nthreads_)
	for (int n=0; n < nthreads_; ++n) {
		for (uint32_t i=n; i < classes.size(); i+=nthreads_) {
			process_block_SIN(classes[i], blocks[n]);
		}
	}
	
	cout << "Centering finished."  << endl;
	
	
	//ofstream outf; outf.open(dirprefix + "/reads.uf.corr");
	//size_t blockNum = 0;
	for (int n=0; n < nthreads_; ++n) {
		for (uint32_t i=0; i < blocks[n].size(); ++i) {
			if (blocks[n][i].size() == 0) continue;
			if (blocks[n][i].size() == 1) {
				//cout << "  kmer " << blocks[n][i][0] << " is good with count " << k_[blocks[n][i][0]].second.count << endl;
				if (k_[blocks[n][i][0]].second.count > GOOD_SINGLETON_THRESHOLD) {
					k_[blocks[n][i][0]].second.change = false;
					k_[blocks[n][i][0]].second.good = true;
				}
				//outf << blockNum << "\t" << k_[blocks[n][i][0]].first.str() << "\t" << k_[blocks[n][i][0]].second.count << "\t" << blocks[n][i].size() << "\t0.0\tgoodSingleton\t0" << endl;
			} else {
				//cout << "  kmer " << blocks[n][i][0] << " is good as a center with count " << k_[blocks[n][i][0]].second.count << endl;
				k_[blocks[n][i][0]].second.change = false;
				k_[blocks[n][i][0]].second.good = true;
				//outf << blockNum << "\t" << k_[blocks[n][i][0]].first.str() << "\t" << k_[blocks[n][i][0]].second.count << "\t" << blocks[n][i].size() << "\t0.0\tcenter\t0" << endl;
				for (uint32_t j=1; j < blocks[n][i].size(); ++j) {
					//cout << "  kmer " << blocks[n][i][j] << " should be changed to " << blocks[n][i][0] << endl;
					//cout << "       " << k_[blocks[n][i][j]].first.str() << "\n       " << k_[blocks[n][i][0]].first.str() << endl;

					//outf << blockNum << "\t" << k_[blocks[n][i][j]].first.str() << "\t" << k_[blocks[n][i][j]].second.count << "\t" << blocks[n][i].size() << "\t0.0\tchange\t0" << endl;
					k_[blocks[n][i][j]].second.good = false;
					k_[blocks[n][i][j]].second.change = true;
					k_[blocks[n][i][j]].second.changeto = blocks[n][i][0];
				}
			}
			//++blockNum; outf << endl;
		}
	}
	//outf.close();
}

