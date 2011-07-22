/*
 * hammer_tools.cpp
 *
 *  Created on: 08.07.2011
 *      Author: snikolenko
 */

#include <omp.h>
#include <iostream>
#include <fstream>
#include <boost/bind.hpp>
#include "read/ireadstream.hpp"
#include "defs.hpp"
#include "mathfunctions.hpp"
#include "valid_kmer_generator.hpp"
#include "position_kmer.hpp"

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
		} else {
			PositionKMer pkm = it->first;
			KMerStat kms; kms.count = it->second.count;
			v1.insert( make_pair( pkm, kms ) );
		}
	}
}


size_t ReadStatMapContainer::size() {
	size_t res = 0;
	for (size_t i=0; i < v_.size(); ++i) {
		res += v_[i].size();
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
	KMerStatMap::const_iterator imin = cur_min();
	if (imin == v_[0].end() ) {
		//cout << "!end" << endl;
		KMerCount p( make_pair(PositionKMer(0,0), KMerStat()) );
		p.second.count = MAX_INT_64;
		//cout << "!!end" << endl;
		return p;
	}
	KMerCount p = *imin;
	p.first = imin->first; p.second.count = 0;
	for (size_t i=0; i<v_.size(); ++i) {
		if (i_[i]->first == p.first) {
			p.second.count += i_[i]->second.count;
			p.second.pos.insert(p.second.pos.end(), i_[i]->second.pos.begin(), i_[i]->second.pos.end());
			++i_[i];
		}
	}
	return p;
}

const KMerStatMap::const_iterator & ReadStatMapContainer::cur_min() {
	int min=0; //cout << (i_[min] == v_[min].end()) << endl;
	for (size_t i=1; i<v_.size(); ++i) {
		if (i_[i] != v_[i].end()) {
			if (i_[min] == v_[min].end() || i_[i]->first < i_[min]->first) {
				min = i;
			}
		}
	}
	//cout << min << "  " << (i_[min] == v_[min].end()) << endl;
	return i_[min];
}

/**
 * add k-mers from read to map
 */
template<uint32_t kK, typename KMerStatMap>
void AddKMers(const Read &r, uint64_t readno, KMerStatMap *v) {
	ValidKMerGenerator<K> gen(r);
	while (gen.HasMore()) {
		PositionKMer pkm(readno, gen.pos() - 1);
		++(*v)[pkm].count;
		(*v)[pkm].pos.push_back( make_pair(readno, gen.pos() - 1) );
		gen.Next();
	}
}

void DoPreprocessing(int tau, int qvoffset, string readsFilename, int nthreads, vector<KMerStatMap> * vv, vector<ReadStat> * rv) {
	vv->clear();
	//for (int i=0; i<nthreads; ++i) {
		KMerStatMap v; v.clear();
		vv->push_back(v);
	//}

	cout << "Starting preproc.\n";

	// TODO: think about a parallelization -- for some reason, the previous version started producing segfaults

	for(int i=0; i < rv->size(); ++i) {
		AddKMers<K, KMerStatMap>(rv->at(i).read, i, &(vv->at(0)));
	}
	
	cout << "All k-mers added to maps." << endl;
	//for (int i=0; i<nthreads; ++i) {
		cout << "map " << 0 << " size = " << vv->at(0).size() << endl;
	//}

}

void DoSplitAndSort(int tau, int nthreads, ReadStatMapContainer & rsmc, vector< vector<uint64_t> > * vs, vector<KMerCount> * kmers, vector<ReadStat> * rv) {
	int effective_threads = min(nthreads, tau+1);	
	uint64_t kmerno = 0;
	kmers->clear();
	cout << "Starting split and sort..." << endl;

	for (KMerCount p = rsmc.next(); p.second.count < MAX_INT_64; p = rsmc.next()) {
		kmers->push_back(p);
		for (uint32_t j=0; j<p.second.pos.size(); ++j) {
			rv->at(p.second.pos[j].first).kmers.insert( make_pair(p.second.pos[j].second, kmerno) );
		}
		++kmerno;
	}
	cout << "Auxiliary vectors loaded. In total, we have " << kmerno << " kmers." << endl;

	#pragma omp parallel for shared(vs, kmers, tau) num_threads(effective_threads)
	for (int j=0; j<tau+1; ++j) {
		vs->at(j).resize( kmers->size() );
		for (size_t m = 0; m < kmers->size(); ++m) vs->at(j)[m] = m;

		sort(vs->at(j).begin(), vs->at(j).end(), boost::bind(PositionKMer::compareSubKMers, _1, _2, kmers, tau, j));
		//cout << "SubVector " << j << ":" << endl;
		//for (size_t m = 0; m < vs->at(j).size(); ++m ) cout << "  " << vs->at(j)[m] << "  " << kmers->at(vs->at(j)[m]).first.strSub(tau, j) << endl;
	}
	cout << "Auxiliary vectors sorted." << endl;
}



bool CorrectRead(const vector<KMerCount> & km, ReadStat * r, ReadStat * r_rev, ofstream * ofs) {
	string seq = r->read.getSequenceString();
//	cout << "Correcting " << r->read.getName() << "\n" << seq.data() << endl;

	// create auxiliary structures for consensus
	vector<int> vA(r->read.size(), 0), vC(r->read.size(), 0), vG(r->read.size(), 0), vT(r->read.size(), 0);
	vector< vector<int> > v;  // A=0, C=1, G=2, T=3
	v.push_back(vA); v.push_back(vC); v.push_back(vG); v.push_back(vT);

	bool changedRead = false;
	for (map<uint32_t, uint64_t>::const_iterator it = r->kmers.begin(); it != r->kmers.end(); ++it) {
		const PositionKMer & kmer = km[it->second].first;
		const uint32_t pos = it->first;
		const KMerStat & stat = km[it->second].second;

		if (stat.good == true) {
			for (size_t j=0; j<K; ++j) {
				v[dignucl(kmer[j])][pos+j]++;
			}

		} else {
			if (stat.change == true) {
				const PositionKMer & newkmer = km[ stat.changeto ].first;
				
//				for (size_t j=0; j<pos; ++j) cout << " ";
//				cout << newkmer.str().data() << "\n";
				
				for (size_t j=0; j<K; ++j) {
					v[dignucl(newkmer[j])][pos+j]++;
				}
				// pretty print the k-mer
				if (!changedRead) {
					changedRead = true;
					if (ofs != NULL) *ofs << "\n " << r->read.getName() << "\n" << r->read.getSequenceString().data() << "\n";
				}
				if (ofs != NULL) {
					for (size_t j=0; j<pos; ++j) *ofs << " ";
					*ofs << newkmer.str().data() << "\n";
				}
			}
		}
	}

	for (map<uint32_t, uint64_t>::const_iterator it = r_rev->kmers.begin(); it != r_rev->kmers.end(); ++it) {
		const PositionKMer & kmer = km[it->second].first;
		const uint32_t pos = it->first;
		const KMerStat & stat = km[it->second].second;
		if (stat.good == true) {
			for (size_t j=0; j<K; ++j) {
				v[complement(dignucl(kmer[j]))][r->read.size()-pos-j-1]++;
			}
			/*if (print_debug) {
				for (size_t j=0; j<r->read.size()-pos-K; ++j) cout << " ";
				cout << (!kmer).str().data() << "\n";
			}*/
		} else {
			if (stat.change == true) {
				const PositionKMer & newkmer = km[ stat.changeto ].first;

//				for (size_t j=0; j<r->read.size()-pos-K; ++j) cout << " ";
//				for (size_t j=0; j<K; ++j) cout << nucl_complement( newkmer[K-j-1] );
//				cout << endl;

				for (size_t j=0; j<K; ++j) {
					v[complement(dignucl(newkmer[j]))][r->read.size()-pos-j-1]++;
				}
				if (!changedRead) {
					changedRead = true;
					if (ofs != NULL) *ofs << "\n " << r->read.getName() << "\n" << r->read.getSequenceString().data() << "\n";
				}
				if (ofs != NULL) {
					for (size_t j=0; j<r->read.size()-pos-K; ++j) *ofs << " ";
					for (size_t j=0; j<K; ++j) *ofs << nucl_complement(newkmer[K-j-1]);
					*ofs << "\n";
				}
			}
		}
	}

	// at this point the array v contains votes for consensus

	bool res = false; // has anything really changed?
	// find max consensus element
	for (size_t j=0; j<r->read.size(); ++j) {
		char cmax = seq[j]; int nummax = 0;
		for (size_t k=0; k<4; ++k) {
			if (v[k][j] > nummax) {
				cmax = nucl(k); nummax = v[k][j];
			}
		}
		if (seq[j] != cmax) res = true;
		seq[j] = cmax;
	}
	
	// print consensus array
	if (ofs != NULL && changedRead) {
		for (size_t i=0; i<4; ++i) {
			for (size_t j=0; j<r->read.size(); ++j) {
				*ofs << (char)((int)'0' + v[i][j]);
			}
			*ofs << "\n";
		}
		*ofs << seq.data() << "\n";
	}
	
	// change the read
	r->read.setSequence(seq.data());
	return res;
}

