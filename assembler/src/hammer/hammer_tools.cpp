/*
 * hammer_tools.cpp
 *
 *  Created on: 08.07.2011
 *      Author: snikolenko
 */

#include <omp.h>
#include <iostream>
#include <fstream>
#include "read/ireadstream.hpp"
#include "hammer/kmer_functions.hpp"
#include "defs.hpp"
#include "hammerread.hpp"
#include "mathfunctions.hpp"
#include "hammer_tools.hpp"

string encode3toabyte (const string & s)  {
	string retval;
	char c = 48;
	int weight = 16;
	size_t i;
	for (i = 0; i < s.length(); i += 1) {
		if (i % 3 == 0) {
			c= 48;
			weight = 16;
		}
		c += weight * nt2num(s[i]);
		weight /= 4;
		if (i % 3 == 2) retval += c;
	}
	if (i % 3 != 0) retval += c;
	return retval;
}

void join_maps(KMerStatMap & v1, const KMerStatMap & v2) {
	KMerStatMap::iterator itf;
	for (KMerStatMap::const_iterator it = v2.begin(); it != v2.end(); ++it) {
		itf = v1.find(it->first);
		if (itf != v1.end()) {
			itf->second.count = itf->second.count + it->second.count;
			itf->second.freq  = itf->second.freq  + it->second.freq;
		} else {
			pair<KMer, KMerStat> p;
			p.first = it->first;
			p.second.count = it->second.count;
			p.second.freq = it->second.freq;
			v1.insert(p);
		}
	}
}


size_t ReadStatMapContainer::size() {
	size_t res = 0;
	for (size_t i=0; i < v_.size(); ++i) {
		res += v_.size();
	}
	return res;
}
void ReadStatMapContainer::init() {
	i_.clear();
	for (size_t i=0; i<v_.size(); ++i) {
		i_.push_back(v_[i].begin());
	}
}

pair<KMer, KMerStat> ReadStatMapContainer::next() {
	pair<KMer, KMerStat> p;
	KMerStatMap::const_iterator imin = cur_min();
	if (imin == v_[0].end()) {
		p.second.count = MAX_INT_64;
		return p;
	}
	p.first = imin->first; p.second.count = 0; p.second.freq = 0;
	for (size_t i=0; i<v_.size(); ++i) {
		if (i_[i]->first == p.first) {
			p.second.count += i_[i]->second.count;
			p.second.freq  += i_[i]->second.freq;
			++i_[i];
		}
	}
	return p;
}

const KMerStatMap::const_iterator & ReadStatMapContainer::cur_min() {
	int min=0;
	for (size_t i=1; i<v_.size(); ++i) {
		if (i_[i] != v_[i].end()) {
			if (i_[min] == v_[min].end() || KMerLess(i_[i]->first, i_[min]->first)) {
				min = i;
			}
		}
	}
	return i_[min];
}

void ClusterProcessBlock(unionFindClass * uf, StringKMerVector & block, int tau) {
	for (int i = 0; i < block.size(); i++) {
		uf->find_set(block[i].count);
		for (int j = i + 1; j < block.size(); j++) {
			if (hamdist(block[i].kmer.str(), block[j].kmer.str(), SAME_STRAND, tau ) <= tau) {
				uf->unionn(block[i].count, block[j].count);
			}
		}
	}
	return;
}

void DoPreprocessing(int tau, int qvoffset, string readsFilename, int nthreads, vector<KMerStatMap> & vv) {
	for (int i=0; i<nthreads; ++i) {
		KMerStatMap v; v.clear();
		vv.push_back(v);
	}

	cout << "Starting preproc.\n";
	ireadstream ifs(readsFilename.data(), qvoffset);
	ofstream ofs;
	Read r;
	size_t tmpc = 0;
	size_t cur_maps = 0;
	vector<Read> rv;
	while (!ifs.eof()) {
		// reading a batch of reads
		for (int thr = 0; thr < READ_BATCH_SIZE; ++thr) {
			ifs >> r; 
			if (TrimBadQuality(&r) >= K) {
				rv.push_back(r);
			}
			if (ifs.eof()) break;
		}

		// trim the reads for bad quality and process only the ones with at least K "reasonable" elements
		// we now do it in parallel

		++tmpc;
		cout << "Batch " << tmpc << " read.\n"; flush(cout);
		#pragma omp parallel for shared(rv, vv, ofs) num_threads(nthreads)
		for(int i=0; i<rv.size(); ++i) {
			AddKMers(rv[i], &vv[omp_get_thread_num() + cur_maps * nthreads]);
			AddKMers(!(rv[i]), &vv[omp_get_thread_num() + cur_maps * nthreads]);
		}
		cout << "Batch " << tmpc << " added.\n"; flush(cout);
		rv.clear();
	}
	ifs.close();
	cout << "All k-mers added to maps.\n"; flush(cout);

	for (int i=0; i<vv.size(); ++i) {
		cout << "size(" << i << ")=" << vv[i].size() << "\n"; flush(cout);
	}
}

void DoSplitAndSort(int tau, int nthreads, ReadStatMapContainer & rsmc, vector<StringKMerVector> * vs, vector<KMerCount> * kmers) {
	int effective_threads = min(nthreads, tau+1);	
	size_t counter = 0;
	kmers->clear();
	cout << "hello from dosplitandsort\n";

	for (KMerCount p = rsmc.next(); p.second.count < MAX_INT_64; p = rsmc.next()) {
		kmers->push_back(p);

		#pragma omp parallel for shared(vs, p, counter, tau) num_threads(effective_threads)
		for (int j=0; j<tau+1; ++j) {
			string sub = "";
			for (int i = j; i < K; i += tau+1) {
				sub += p.first[i];
			}
			StringKMer skm; skm.sub = sub; skm.count = counter; skm.kmer = p.first;
			vs->at(j).push_back(skm);
		}
		++counter;
		if (counter % 10000000 == 0) cout << "Split " << counter << " kmers.\n"; flush(cout); 
	}
	cout << "Auxiliary vectors loaded.\n"; flush(cout);

	#pragma omp parallel for shared(vs) num_threads(effective_threads)
	for (int j=0; j<tau+1; ++j) {
		sort(vs->at(j).begin(), vs->at(j).end(), SKgreater);
	}
	cout << "Auxiliary vectors sorted.\n"; flush(cout);		
}

void ClusterMerge(vector<unionFindClass *>uf, const vector<KMerCount> & kmers, unionFindClass * ufMaster) {
	cout << "Merging union find files...\n"; flush(cout);
	vector<string> row;
	vector<vector<int> > classes;
	for (int i = 0; i < uf.size(); i++) {
		classes.clear();
		uf[i]->get_classes(classes);
		delete uf[i];
		for (int j = 0; j < classes.size(); j++) {
			int first = classes[j][0];
			ufMaster->find_set(first);
			for (int k = 0; k < classes[j].size(); k++) {
				ufMaster->unionn(first, classes[j][k]);
			}
		}
	}
}

void DoClustering(int tau, int nthreads, string dirprefix, const vector<StringKMerVector> & vs, const vector<KMerCount> & kmers) {
	
	int effective_threads = min(nthreads, tau+1);	
	vector<unionFindClass *> uf(effective_threads);

	for (int i = 0; i < effective_threads; i++) {
		uf[i] = new unionFindClass(kmers.size()); 

		cout << "Processing split kmers " << i << ", total " << vs[i].size() << "\n";
		
		string sbuf;
		StringKMer last;
		StringKMerVector block;
		for (size_t j=0; j<vs[i].size(); ++j) {
			if (j % 10000000 == 0) cerr << "Processed (" << i << ") " << add_commas(j) << ", ";
			if (last.sub == vs[i][j].sub) { //add to current reads
				block.push_back(vs[i][j]);
			} else {
				ClusterProcessBlock(uf[i], block, tau);
				block.clear();
				block.push_back(vs[i][j]);
				last = vs[i][j];
			}
		}
		ClusterProcessBlock(uf[i], block, tau);
		cout << "Finished(" << i << ") " << endl;
	}
	cout << "All threads finished.\n";flush(cout);
	
	unionFindClass * ufMaster;
	ufMaster = new unionFindClass(kmers.size());
	ClusterMerge(uf, kmers, ufMaster);
	cout << "Merging finished.\n"; flush(cout);
	
	
	cout << "Centering begins...\n"; flush(cout);
	vector<vector<int> > classes;
	ufMaster->get_classes(classes);
	delete ufMaster;
	vector< vector< vector<HammerRead> > > blocks(nthreads);
	
	#pragma omp parallel for shared(blocks, classes, kmers) num_threads(nthreads)
	for (int n=0; n < nthreads; ++n) {
		for (int i=n; i < classes.size(); i+=nthreads) {
			vector<HammerRead> sinBlock;
			for (int j = 0; j < classes[i].size(); j++) {
				HammerRead hr;
				hr.seq = kmers[classes[i][j]].first.str();
				hr.count = kmers[classes[i][j]].second.count;
				sinBlock.push_back(hr);
			}
			process_block_SIN(sinBlock, blocks[n]);
		}
	}
	
	cout << "Centering finished. Outputting clusters...\n"; flush(cout);
	
	ofstream outf; open_file(outf, dirprefix + "/reads.uf.corr");
	size_t blockNum = 0;
	for (int n=0; n < nthreads; ++n) {
		for (int i=0; i < blocks[n].size(); ++i) {
			if (blocks[n][i].size() == 0) continue;
			if (blocks[n][i].size() == 1) {
				outf << blockNum << "\t" << blocks[n][i][0].seq << "\t" << blocks[n][i][0].count << "\t" << blocks[n][i][0].freq << "\t" << blocks[n][i].size() << "\t0.0\tgoodSingleton\t0" << endl;
			} else {
				outf << blockNum << "\t" << blocks[n][i][0].seq << "\t" << blocks[n][i][0].count << "\t" << blocks[n][i][0].freq << "\t" << blocks[n][i].size() << "\t0.0\tcenter\t0" << endl;
				for (int j=1; j < blocks[n][i].size(); ++j) {
					outf << blockNum << "\t" << blocks[n][i][j].seq << "\t" << blocks[n][i][j].count << "\t" << blocks[n][i][j].freq << "\t" << blocks[n][i].size() << "\t0.0\tchange\t0" << endl;
				}
			}
			++blockNum; outf << endl;
		}
	}
	outf.close();
}


double calcMultCoef(vector<int> & distances, vector<HammerRead> & kmers) {
	int kmersize = kmers[0].seq.size();
	double prob = 0;
	double theta;
	for (size_t i = 0; i < kmers.size(); i++) {
		theta = kmers[i].count * -((kmersize - distances[i]) * log(1 - ERROR_RATE) + distances[i] * log(ERROR_RATE));
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
string find_consensus_with_mask(vector<HammerRead> & block, const vector<int> & mask, int maskVal) {
	string c;
	for (size_t i = 0; i < block[0].seq.length(); i++) {
		int scores[4] = {0,0,0,0};
		for (size_t j = 0; j < block.size(); j++) {
			if (mask[j] == maskVal)
				scores[nt2num(block[j].seq.at(i))] += block[j].count;
		}
		c.push_back(num2nt(argmax(scores, 4)));
	}
	return c;
}


string find_consensus(vector<HammerRead> & block) {
	string c;
	for (size_t i = 0; i < block[0].seq.length(); i++) {
		int scores[4] = {0,0,0,0};
		for (size_t j = 0; j < block.size(); j++) {
			scores[nt2num(block[j].seq.at(i))]+=block[j].count;
		}
		c.push_back(num2nt(argmax(scores, 4)));
	}
	return c;
}



/**
  * @return total log-likelihood of this particular clustering
  */
double clusterLogLikelihood(const vector<HammerRead> & block, const vector<HammerRead> & centers, const vector<int> & indices) {
	if (block.size() == 0) return 0.0;
	assert(block.size() == indices.size());
	assert(centers.size() > 0);
	
	// if there is only one center, there are no beta coefficients
	if (centers.size() == 1) {
		int dist = 0;
		for (size_t i=0; i<block.size(); ++i) {
			dist += hamdist(block[i].seq, centers[indices[i]].seq, SAME_STRAND);
		}
		return ( lMultinomial(block) + log(ERROR_RATE) * dist + log(1-ERROR_RATE) * (block[0].seq.size() * block.size() - dist) );
	}
	
	// compute sufficient statistics
	vector<int> count(centers.size(), 0);		// how many kmers in cluster i
	vector<int> totaldist(centers.size(), 0);	// total distance from kmers of cluster i to its center
	for (size_t i=0; i<block.size(); ++i) {
		count[indices[i]]+=block[i].count;
		totaldist[indices[i]] += hamdist(block[i].seq, centers[indices[i]].seq, SAME_STRAND);
	}
	
	// sum up the likelihood
	double res = lBetaPlusOne(count);   // 1/B(count)
	res += lMultinomial(centers); 		// {sum(centers.count) \choose centers.count}
	for (size_t i=0; i<centers.size(); ++i) {
		res += lMultinomialWithMask(block, indices, i) + 
			   log(ERROR_RATE) * totaldist[i] + log(1-ERROR_RATE) * (block[0].seq.size() * count[i] - totaldist[i]);
	}
	return res;
}


double lMeansClustering(int l, vector< vector<int> > & distances, vector<HammerRead> & kmers, vector<int> & indices, vector<HammerRead> & centers) {
	centers.resize(l); // there are l centers
	// if l==1 then clustering is trivial
	if (l == 1) {
		centers[0].seq = find_consensus(kmers);
		centers[0].count = kmers.size();
		for (size_t i=0; i < kmers.size(); ++i) indices[i] = 0;
		return clusterLogLikelihood(kmers, centers, indices);
	}
	
	int restartCount = 1;
	
	double bestLikelihood = -1000000.0;
	double curLikelihood = -1000000.0;
	vector<HammerRead> bestCenters;
	vector<int> bestIndices(kmers.size());
	
	for (int restart = 0; restart < restartCount; ++restart) {

	// fill in initial approximations
	for (int j=0; j < l; ++j) centers[j].count = 0;
	for (size_t i=0; i < kmers.size(); ++i) {
		for (int j=0; j < l; ++j) {
			if (kmers[i].count > centers[j].count) {
				for (int s=j; s<l-1; ++s) {
					centers[s+1].seq = centers[s].seq;
					centers[s+1].count = centers[s].count;
				}
				centers[j].seq = kmers[i].seq;
				centers[j].count = kmers[i].count;
				break;
			}
		}
	}

	// TODO: make random restarts better!!! they don't really work now
	if (restart > 0) { // introduce random noise
		vector<bool> good(kmers.size(), false);
		for (size_t i=0; i < kmers.size(); ++i) {
			for (int j=0; j < l; ++j) {
				if (kmers[i].count >= centers[j].count) {
					good[i] = true; break;
				}
			}
		}
		for (int k=0; k < l/2 + 1; ++k) {
			int newNo = rand() % (kmers.size() - l - k);
			int indexOld = (rand() % (l-1)) + 1;
			int indexNew = kmers.size()-1;
			for (size_t i=0; i < kmers.size(); ++i) {
				if (!good[i]) --newNo;
				if (newNo < 0) { indexNew = i; good[indexNew] = true; break; }
			}
			centers[indexOld].seq = kmers[indexNew].seq;
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
		for (int j=0; j < l; ++j) centers[j].count = 0;
		
		// E step: find which clusters we belong to
		for (size_t i=0; i < kmers.size(); ++i) {
			for (int j=0; j < l; ++j)
				dists[j] = hamdist(kmers[i].seq, centers[j].seq, ANY_STRAND);
			int newInd = argmin(dists);
			if (indices[i] != newInd) {
				changed = true;
				changedCenter[indices[i]] = true;
				changedCenter[newInd] = true;
				indices[i] = newInd;
			}
			++centers[indices[i]].count;
		}
		
		// M step: find new cluster centers
		for (int j=0; j < l; ++j) {
			if (!changedCenter[j]) continue; // nothing has changed
			centers[j].seq = find_consensus_with_mask(kmers, indices, j);
		}
	}

	// last M step
	for (int j=0; j < l; ++j) {
		centers[j].seq = find_consensus_with_mask(kmers, indices, j);
	}
	
	curLikelihood = clusterLogLikelihood(kmers, centers, indices);
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

void process_block_SIN(vector<HammerRead> & block, vector< vector<HammerRead> > & vec) {
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
			distances[i][j] = hamdist(block[i].seq, block[j].seq, SAME_STRAND);
			distances[j][i] = distances[i][j];
		}
		//cout << block[i].seq << "\t" << block[i].count << "\n";
	}

	/*for (size_t i = 0; i < block.size(); i++) {
		for (size_t j = 0; j < block.size(); j++) {
			cout << distances[i][j] << " ";
		}
		cout << "\n";
	}*/
	
	// Multinomial coefficients -- why? TODO: remove
	for (int i = 0; i < origBlockSize; i++) multiCoef[i] = calcMultCoef(distances[i], block);

	vector<int> indices(origBlockSize);
	double bestLikelihood = -1000000;
	vector<HammerRead> bestCenters;
	vector<int> bestIndices(block.size());

	for (int l = 1; l <= origBlockSize; ++l) {
		vector<HammerRead> centers(l);
		double curLikelihood = lMeansClustering(l, distances, block, indices, centers);
		//cout << "   ...likelihood with " << l << " clusters = " << curLikelihood << "\n";
		if (curLikelihood <= bestLikelihood) {
			break;
		}
		bestLikelihood = curLikelihood;
		bestCenters = centers; bestIndices = indices;
	}
	
	for (size_t k=0; k<bestCenters.size(); ++k) {
		if (bestCenters[k].count == 0) {
			continue; // superfluous cluster with zero elements
		}
		//cout << "subcluster with " << bestCenters[k].count << " elements\n";
		vector<HammerRead> v;
		if (bestCenters[k].count == 1) {
			for (int i = 0; i < origBlockSize; i++) {
				if (indices[i] == (int)k) {
					v.push_back(block[i]);
					break;
				}
			}
		} else { // there are several kmers in this cluster
			//cout << "center: " << bestCenters[k].seq << "\n";
			bool centerInCluster = false;
			for (int i = 0; i < origBlockSize; i++) {
				if (bestIndices[i] == (int)k) {
					int dist = hamdist(block[i].seq, bestCenters[k].seq, SAME_STRAND);
					//cout << "        " << block[i].seq << " " << dist << "\n";
					if (dist == 0) {
						//cout << "  found center\n";
						centerInCluster = true;
						v.insert(v.begin(), block[i]);
					} else {
						v.push_back(block[i]);
					}
				}
			}
			if (!centerInCluster) {
				//cout << "  pushing consensus\n";
				HammerRead consensus; consensus.seq = bestCenters[k].seq; consensus.count = 0; consensus.freq = 0.0;
				v.insert(v.begin(), consensus);
			}
		}
		vec.push_back(v);
	}
}
