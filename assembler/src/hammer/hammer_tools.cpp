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
#include "defs.hpp"
#include "hammerread.hpp"
#include "mathfunctions.hpp"
#include "valid_kmer_generator.hpp"
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
			pair<KMer, KMerStat> p;
			p.first = it->first;
			p.second.count = it->second.count;
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

/**
 * add k-mers from read to map
 */
template<uint32_t kK, typename KMerStatMap>
void AddKMers(const Read &r, uint64_t readno, KMerStatMap *v) {
	ValidKMerGenerator<K> gen(r);
	while (gen.HasMore()) {
		++(*v)[gen.kmer()].count;
		(*v)[gen.kmer()].pos.push_back( make_pair(readno, gen.pos() - 1) );
		gen.Next();
	}
}

void DoPreprocessing(int tau, int qvoffset, string readsFilename, int nthreads, vector<KMerStatMap> * vv, vector<ReadStat> * rv) {
	vv->clear();
	for (int i=0; i<nthreads; ++i) {
		KMerStatMap v; v.clear();
		vv->push_back(v);
	}

	cout << "Starting preproc.\n";
	#pragma omp parallel for shared(rv, vv) num_threads(nthreads)
	for(int64_t i=0; i< (rv->size()); ++i) {
		// we have to enumerate reads from 1 rather than 0 because we need negative numbers for reverse complements
		AddKMers<K, KMerStatMap>(rv->at(i).read, i+1, &(vv->at(omp_get_thread_num())));
		AddKMers<K, KMerStatMap>(!(rv->at(i).read), -(i+1), &(vv->at(omp_get_thread_num())));
	}
	
	cout << "All k-mers added to maps." << endl;

	/*for (uint32_t i=0; i<vv->size(); ++i) {
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
	}*/
}

void DoSplitAndSort(int tau, int nthreads, ReadStatMapContainer & rsmc, vector<StringKMerVector> * vs, vector<KMerCount> * kmers, vector<ReadStat> * rv) {
	int effective_threads = min(nthreads, tau+1);	
	size_t kmerno = 0;
	kmers->clear();
	cout << "Starting split and sort...\n";

	for (KMerCount p = rsmc.next(); p.second.count < MAX_INT_64; p = rsmc.next()) {
		kmers->push_back(p);

		for (uint32_t j=0; j<p.second.pos.size(); ++j) {
			if (p.second.pos[j].first > 0) {
				rv->at(p.second.pos[j].first-1).kmers.insert( make_pair(p.second.pos[j].second, kmerno) );
				// cout << p.second.pos[j].first << ":" << p.second.pos[j].second << ":" << rv->at(p.second.pos[j].first-1).kmers[ p.second.pos[j].second ] << " ";
			} else {
				rv->at(-p.second.pos[j].first-1).kmers_rev.insert( make_pair(p.second.pos[j].second, kmerno) );
				// cout << p.second.pos[j].first << ":" << p.second.pos[j].second << ":" << rv->at(-p.second.pos[j].first-1).kmers_rev[ p.second.pos[j].second ] << " ";
			}
		}
		//cout << endl;
		// }
		
		//cout << kmerno << "   " << p.first.str() << " ";

		#pragma omp parallel for shared(vs, p, kmerno, tau) num_threads(effective_threads)
		for (int j=0; j<tau+1; ++j) {
			string sub = "";
			for (int i = j; i < K; i += tau+1) {
				sub += p.first[i];
			}
			StringKMer skm; skm.sub = sub; skm.count = kmerno; skm.kmer = p.first;
			vs->at(j).push_back(skm);
			//for (int m=0; m< vs->at(j)[vs->at(j).size() - 1].sub.size(); ++m) cout << nucl(vs->at(j)[vs->at(j).size() - 1].sub[m]);
			//cout << "  ";
		}
		//cout << endl;
		++kmerno;
		if (kmerno % 10000000 == 0) cout << "Split " << kmerno << " kmers." << endl;
	}
	cout << "Auxiliary vectors loaded." << endl;

	#pragma omp parallel for shared(vs) num_threads(effective_threads)
	for (int j=0; j<tau+1; ++j) {
		sort(vs->at(j).begin(), vs->at(j).end(), SKgreater);
		/*cout << "vs[" << j << "]:" << endl;
		for (int l=0; l<vs->at(j).size(); ++l) {
			cout << "    ";
			for (int m=0; m< vs->at(j)[l].sub.size(); ++m) cout << nucl(vs->at(j)[l].sub[m]);
			cout << "  " << vs->at(j)[l].count << endl;
		}*/
	}
	cout << "Auxiliary vectors sorted." << endl;
}



bool CorrectRead(const map<KMer, KMer, KMer::less2> & changes, const unordered_set<KMer, KMer::hash> & good, const vector<KMerCount> & km, ReadStat * r, ofstream * ofs) {
	string seq = r->read.getSequenceString();
	//if (print_debug) cout << "Correcting " << r->read.getName() << "\n" << seq.data() << "\n";
	// create auxiliary structures for consensus
	vector<int> vA(r->read.size(), 0), vC(r->read.size(), 0), vG(r->read.size(), 0), vT(r->read.size(), 0);
	vector< vector<int> > v;  // A=0, C=1, G=2, T=3
	v.push_back(vA); v.push_back(vC); v.push_back(vG); v.push_back(vT);

	bool changedRead = false;
	for (map<uint32_t, uint64_t>::const_iterator it = r->kmers.begin(); it != r->kmers.end(); ++it) {
		const KMer & kmer = km[it->second].first;
		const uint32_t pos = it->first;
		const KMerStat & stat = km[it->second].second;
		if (stat.good == true) {
			for (size_t j=0; j<K; ++j) {
				v[kmer[j]][pos+j]++;
			}
			/*if (print_debug) {
				for (size_t j=0; j<pos; ++j) cout << " ";
				cout << kmer.str().data() << "\n";
			}*/
		} else {
			if (stat.change == true) {
				const KMer & newkmer = km[ stat.changeto ].first;
				for (size_t j=0; j<K; ++j) {
					v[newkmer[j]][pos+j]++;
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

	for (map<uint32_t, uint64_t>::const_iterator it = r->kmers_rev.begin(); it != r->kmers_rev.end(); ++it) {
		const KMer & kmer = km[it->second].first;
		const uint32_t pos = it->first;
		const KMerStat & stat = km[it->second].second;
		if (stat.good == true) {
			for (size_t j=0; j<K; ++j) {
				v[complement(kmer[j])][r->read.size()-pos-j-1]++;
			}
			/*if (print_debug) {
				for (size_t j=0; j<r->read.size()-pos-K; ++j) cout << " ";
				cout << (!kmer).str().data() << "\n";
			}*/
		} else {
			if (stat.change == true) {
				const KMer & newkmer = km[ stat.changeto ].first;
				for (size_t j=0; j<K; ++j) {
					v[complement(newkmer[j])][r->read.size()-pos-j-1]++;
				}
				if (!changedRead) {
					changedRead = true;
					if (ofs != NULL) *ofs << "\n " << r->read.getName() << "\n" << r->read.getSequenceString().data() << "\n";
				}
				if (ofs != NULL) {
					for (size_t j=0; j<r->read.size()-pos-K; ++j) *ofs << " ";
					*ofs << (!newkmer).str().data() << "\n";
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

