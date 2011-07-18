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

KMerCount ReadStatMapContainer::next() {
	KMerCount p;
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
			p.second.pos.insert(p.second.pos.end(), i_[i]->second.pos.begin(), i_[i]->second.pos.end());
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

void DoPreprocessing(int tau, int qvoffset, string readsFilename, int nthreads, vector<KMerStatMap> * vv, vector<ReadStat> * rv) {
	vv->clear();
	for (int i=0; i<nthreads; ++i) {
		KMerStatMap v; v.clear();
		vv->push_back(v);
	}

	cout << "Starting preproc.\n";
	#pragma omp parallel for shared(rv, vv) num_threads(nthreads)
	for(uint64_t i=0; i< (rv->size()); ++i) {
		if (TrimBadQuality(&(rv->at(i).read)) < K) continue;  /// maybe we should get the quality into AddKMers
		AddKMers<K, KMerStatMap>(rv->at(i).read, i, &(vv->at(omp_get_thread_num())));
		AddKMers<K, KMerStatMap>(!(rv->at(i).read), i, &(vv->at(omp_get_thread_num())));
	}
	
	cout << "All k-mers added to maps." << endl;

	for (uint32_t i=0; i<vv->size(); ++i) {
		cout << "size(" << i << ")=" << vv->at(i).size() << endl;
		#pragma omp critical
		{
		for (KMerStatMap::iterator it = vv->at(i).begin(); it != vv->at(i).end(); ++it) {
			cout << it->first.str() << " ";
			for (uint32_t j=0; j<it->second.pos.size(); ++j)
				cout << it->second.pos[j].first << ":" << it->second.pos[j].second << " ";
			cout << endl;
		}
		}
	}
}

void DoSplitAndSort(int tau, int nthreads, ReadStatMapContainer & rsmc, vector<StringKMerVector> * vs, vector<KMerCount> * kmers, vector<ReadStat> * rv) {
	int effective_threads = min(nthreads, tau+1);	
	size_t kmerno = 0;
	kmers->clear();
	cout << "Starting split and sort...\n";

	for (KMerCount p = rsmc.next(); p.second.count < MAX_INT_64; p = rsmc.next()) {
		kmers->push_back(p);

		#pragma omp critical
		{
		cout << kmerno << "   " << p.first.str() << " ";
		for (uint32_t j=0; j<p.second.pos.size(); ++j) {
			rv->at(p.second.pos[j].first).kmers.insert( make_pair(p.second.pos[j].second, kmerno) );
			cout << p.second.pos[j].first << ":" << rv->at(p.second.pos[j].first).kmers[ p.second.pos[j].second ] << " ";
		}
		cout << endl;
		}
		
		#pragma omp parallel for shared(vs, p, kmerno, tau) num_threads(effective_threads)
		for (int j=0; j<tau+1; ++j) {
			string sub = "";
			for (int i = j; i < K; i += tau+1) {
				sub += p.first[i];
			}
			StringKMer skm; skm.sub = sub; skm.count = kmerno; skm.kmer = p.first;
			vs->at(j).push_back(skm);
		}
		++kmerno;
		if (kmerno % 10000000 == 0) cout << "Split " << kmerno << " kmers." << endl; 
	}
	cout << "Auxiliary vectors loaded." << endl;

	#pragma omp parallel for shared(vs) num_threads(effective_threads)
	for (int j=0; j<tau+1; ++j) {
		sort(vs->at(j).begin(), vs->at(j).end(), SKgreater);
	}
	cout << "Auxiliary vectors sorted." << endl;
}

void CorrectRead(const map<KMer, KMer, KMer::less2> & changes, const unordered_set<KMer, KMer::hash> & good, Read * r, bool print_debug) {
	string seq = r->getSequenceString();
	if (print_debug) cout << "Correcting " << r->getName() << "\n" << seq.data() << "\n";
	// create auxiliary structures for consensus
	vector<int> vA(r->size(), 0), vC(r->size(), 0), vG(r->size(), 0), vT(r->size(), 0);
	vector< vector<int> > v;  // A=0, C=1, G=2, T=3
	v.push_back(vA); v.push_back(vC); v.push_back(vG); v.push_back(vT);
			
	int32_t pos = -1;
	bool changedRead = false;
	KMer kmer;
	map<KMer, KMer, KMer::less2>::const_iterator it;
	unordered_set<KMer, KMer::hash>::const_iterator it_single;
	while ( (pos = NextValidKmer<K>(*r, pos, kmer)) >= 0 ) {
		it_single = good.find(kmer);
		if (it_single != good.end()) { //it's a good singleton
			for (size_t j=0; j<K; ++j) {
				v[kmer[j]][(pos-1)+j]++;
				//cout << nucl(kmer[j]);
			}
			//for (size_t j=0; j<i; ++j) cout << " ";
			//cout << kmer.str().data() << "\n";
		} else {
			it = changes.find(kmer);				
			if (it != changes.end()) { // we've found that we have to change this kmer
				for (size_t j=0; j<K; ++j) {
					v[it->second[j]][(pos-1)+j]++;
				}
				// pretty print the k-mer
				if (!changedRead) {
					changedRead = true;
					if (print_debug) cout << "\n " << r->getName() << "\n" << r->getSequenceString().data() << "\n";
				}
				if (print_debug) {
					for (size_t j=0; j<(pos-1); ++j) cout << " ";
					cout << it->second.str().data() << "\n";
				}
			}		
		}
	}
					
	// at this point the array v contains votes for consensus
		
	// find max consensus element
	for (size_t j=0; j<r->size(); ++j) {
		char cmax = seq[j]; int nummax = 0;
		for (size_t k=0; k<4; ++k) {
			if (v[k][j] > nummax) {
				cmax = nucl(k); nummax = v[k][j];
			}
		}
		seq[j] = cmax;
	}
	
	// print consensus array
	if (print_debug && changedRead) {
		for (size_t i=0; i<4; ++i) {
			for (size_t j=0; j<r->size(); ++j) {
				if (v[i][j] > 9) cout << "*"; else cout << v[i][j];
			}
			cout << "\n";
		}
		cout << seq.data() << "\n";
	}
	
	// change the read
	r->setSequence(seq.data());
}
