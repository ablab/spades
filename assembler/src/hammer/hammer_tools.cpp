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
			KMerStat kms( it->second.count, KMERSTAT_GOOD );
			v1.insert( make_pair( pkm, kms ) );
		}
	}
}

/**
 * add k-mers from read to map
 */
template<uint32_t kK, typename KMerStatMap>
void AddKMers(const PositionRead &r, uint64_t readno, KMerStatMap *v) {
	string s = r.getSequenceString();
	ValidKMerGenerator<K> gen(r, s);
	while (gen.HasMore()) {
		PositionKMer pkm(readno, gen.pos() - 1);
		++(*v)[pkm].count;
		(*v)[pkm].pos.push_back( make_pair(readno, gen.pos() - 1) );
		gen.Next();
	}
}

/**
 * add k-mers from read to vector
 */
void AddKMerNos(const PositionRead &r, uint64_t readno, vector<KMerNo> *v) {
	string s = r.getSequenceString();
	ValidKMerGenerator<K> gen(r, s);
	while (gen.HasMore()) {
		v->push_back( PositionKMer::pr->at(readno).start() + gen.pos() - 1 );
		gen.Next();
	}
}

void DoPreprocessing(int tau, int qvoffset, string readsFilename, int nthreads, vector<KMerNo> * vv) {
	vv->clear();
	cout << "Starting preproc. " << PositionKMer::pr->size() << " reads.\n";

	// TODO: think about a parallelization -- for some reason, the previous version started producing segfaults

	for(size_t i=0; i < PositionKMer::pr->size(); ++i) {
		AddKMerNos(PositionKMer::pr->at(i), i, vv);
		if ( i % 1000000 == 0 ) cout << "Processed " << i << " reads." << endl;
	}	
	cout << "All k-mers added to maps." << endl;
}

void DoSplitAndSort(int tau, int nthreads, const vector<KMerNo> & vv, vector< vector<uint64_t> > * vs, vector<KMerCount> * kmers) {
	int effective_threads = min(nthreads, tau+1);
	uint64_t kmerno = 0;
	kmers->clear();
	cout << "Starting split and sort..." << endl;

	KMerNo curKMer = vv[0];
	KMerCount curKMerCount = make_pair( PositionKMer(vv[0].index), KMerStat(0, KMERSTAT_GOOD) );
	for (uint64_t i = 0; i < vv.size(); ++i) {
		if ( !curKMer.equal(vv[i]) ) {
			kmers->push_back(curKMerCount);
			curKMer = vv[i];
			curKMerCount = make_pair( PositionKMer(vv[i].index), KMerStat(0, KMERSTAT_GOOD) );
			++kmerno;
		}
		curKMerCount.second.count++;
		uint64_t readno = PositionKMer::readNoFromBlobPos( vv[i].index );
		PositionKMer::pr->at(readno).kmers().insert( make_pair( vv[i].index - PositionKMer::pr->at(readno).start(), kmerno ) );
	}
	kmers->push_back(curKMerCount);

	cout << "Auxiliary vectors loaded. In total, we have " << kmers->size() << " kmers." << endl;

	#pragma omp parallel for shared(vs, kmers, tau) num_threads(effective_threads)
	for (int j=0; j<tau+1; ++j) {
		vs->at(j).resize( kmers->size() );
		for (size_t m = 0; m < kmers->size(); ++m) vs->at(j)[m] = m;

		cout << j << " " << vs->at(j).size() << " iter=" << (size_t)(&vs->at(j)) << endl;

		sort(vs->at(j).begin(), vs->at(j).end(), boost::bind(PositionKMer::compareSubKMers, _1, _2, kmers, tau, j));
		cout << "Sorted auxiliary vector " << j << endl;
	}
	cout << "Auxiliary vectors sorted." << endl;
}


bool CorrectRead(const vector<KMerCount> & km, uint64_t readno, ofstream * ofs) {
	uint64_t readno_rev = PositionKMer::revNo + readno;
	const Read & r = PositionKMer::rv->at(readno);
	string seq = r.getSequenceString();
	const uint32_t read_size = PositionKMer::rv->at(readno).size();
	PositionRead & pr = PositionKMer::pr->at(readno);
	PositionRead & pr_rev = PositionKMer::pr->at(readno_rev);

	// create auxiliary structures for consensus
	vector<int> vA(read_size, 0), vC(read_size, 0), vG(read_size, 0), vT(read_size, 0);
	vector< vector<int> > v;  // A=0, C=1, G=2, T=3
	v.push_back(vA); v.push_back(vC); v.push_back(vG); v.push_back(vT);

	bool changedRead = false;
	for (map<uint32_t, uint64_t>::const_iterator it = pr.kmers().begin(); it != pr.kmers().end(); ++it) {
		const PositionKMer & kmer = km[it->second].first;
		const uint32_t pos = it->first;
		const KMerStat & stat = km[it->second].second;

		if (stat.changeto == KMERSTAT_GOOD) {
			for (size_t j=0; j<K; ++j) {
				v[dignucl(kmer[j])][pos+j]++;
			}

		} else {
			if (stat.changeto < KMERSTAT_CHANGE) {
				const PositionKMer & newkmer = km[ stat.changeto ].first;
				
				for (size_t j=0; j<K; ++j) {
					v[dignucl(newkmer[j])][pos+j]++;
				}
				// pretty print the k-mer
				if (!changedRead) {
					changedRead = true;
					if (ofs != NULL) *ofs << "\n " << r.getName() << "\n" << seq.data() << "\n";
				}
				if (ofs != NULL) {
					for (size_t j=0; j<pos; ++j) *ofs << " ";
					*ofs << newkmer.str().data() << "\n";
				}
			}
		}
	}

	for (map<uint32_t, uint64_t>::const_iterator it = pr_rev.kmers().begin(); it != pr_rev.kmers().end(); ++it) {
		const PositionKMer & kmer = km[it->second].first;
		const uint32_t pos = it->first;
		const KMerStat & stat = km[it->second].second;
		if (stat.changeto == KMERSTAT_GOOD) {
			for (size_t j=0; j<K; ++j) {
				v[complement(dignucl(kmer[j]))][read_size-pos-j-1]++;
			}
			/*if (print_debug) {
				for (size_t j=0; j<r->read.size()-pos-K; ++j) cout << " ";
				cout << (!kmer).str().data() << "\n";
			}*/
		} else {
			if (stat.changeto < KMERSTAT_CHANGE) {
				const PositionKMer & newkmer = km[ stat.changeto ].first;

				for (size_t j=0; j<K; ++j) {
					v[complement(dignucl(newkmer[j]))][read_size-pos-j-1]++;
				}
				if (!changedRead) {
					changedRead = true;
					if (ofs != NULL) *ofs << "\n " << r.getName() << "\n" << r.getSequenceString().data() << "\n";
				}
				if (ofs != NULL) {
					for (size_t j=0; j<read_size-pos-K; ++j) *ofs << " ";
					for (size_t j=0; j<K; ++j) *ofs << nucl_complement(newkmer[K-j-1]);
					*ofs << "\n";
				}
			}
		}
	}

	// at this point the array v contains votes for consensus

	bool res = false; // has anything really changed?
	// find max consensus element
	for (size_t j=0; j<read_size; ++j) {
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
			for (size_t j=0; j<read_size; ++j) {
				*ofs << (char)((int)'0' + v[i][j]);
			}
			*ofs << "\n";
		}
		*ofs << seq.data() << "\n";
	}
	
	// change the read
	PositionKMer::rv->at(readno).setSequence(seq.data());
	return res;
}

