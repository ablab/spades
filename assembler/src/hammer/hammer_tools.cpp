/*
 * hammer_tools.cpp
 *
 *  Created on: 08.07.2011
 *      Author: snikolenko
 */

#include <omp.h>
#include <iostream>
#include <fstream>
#include <queue>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <time.h>
#include <iomanip>

#include "read/ireadstream.hpp"
#include "defs.hpp"
#include "mathfunctions.hpp"
#include "valid_kmer_generator.hpp"
#include "position_kmer.hpp"
#include "globals.hpp"
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
			KMerStat kms( it->second.count, KMERSTAT_GOOD, it->second.totalQual );
			v1.insert( make_pair( pkm, kms ) );
		}
	}
}

/**
 * add k-mers from read to map
 */
template<uint32_t kK, typename KMerStatMap>
void AddKMers(const PositionRead &r, hint_t readno, KMerStatMap *v) {
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
void AddKMerNos(const PositionRead &r, hint_t readno, vector<KMerNo> *v) {
	string s = r.getSequenceString();
	ValidKMerGenerator<K> gen(r, s);
	while (gen.HasMore()) {
		v->push_back( PositionKMer::pr->at(readno).start() + gen.pos() - 1 );
		gen.Next();
	}
}

void DoPreprocessing(int tau, string readsFilename, int nthreads, vector<KMerNo> * vv) {
	vv->clear();

	vector< vector<KMerNo> > vtmp;
	for(int n=0; n < nthreads; ++n) {
		vector<KMerNo> v_cur;
		vtmp.push_back(v_cur);
	}

	#pragma omp parallel for shared(vtmp) num_threads(nthreads)
	for(size_t i=0; i < PositionKMer::pr->size(); ++i) {
		if (PositionKMer::pr->at(i).bad()) continue;
		AddKMerNos(PositionKMer::pr->at(i), i, &vtmp[omp_get_thread_num()]);
	}
	
	for(int n=0; n < nthreads; ++n) {
		vv->insert(vv->end(), vtmp[n].begin(), vtmp[n].end());
		vtmp[n].clear();
	}
}

void DoSplitAndSort(int tau, int nthreads, vector< vector<hint_t> > * vs, vector<KMerCount> * kmers, vector<SubKMerPQ> * vskpq) {
	int effective_threads = min(nthreads, tau+1);
	int subkmer_nthreads = max ( (tau + 1) * ( (int)(nthreads / (tau + 1)) ), tau+1 );
	int effective_subkmer_threads = min(subkmer_nthreads, nthreads);

	#pragma omp parallel for shared(vs, kmers, tau) num_threads(effective_threads)
	for (int j=0; j<tau+1; ++j) {
		vs->at(j).resize( kmers->size() );
		for (size_t m = 0; m < kmers->size(); ++m) vs->at(j)[m] = m;
	}

	for (int j=0; j < tau+1; ++j) {
		SubKMerCompType sort_routine = boost::bind(SubKMerPQElement::compareSubKMerPQElements, _1, _2, kmers, tau, PositionKMer::subKMerPositions->at(j), PositionKMer::subKMerPositions->at(j+1));
		SubKMerPQ skpq( &(vs->at(j)), max( (int)(nthreads / (tau + 1)), 1), sort_routine );
		vskpq->push_back(skpq);
	}

	// we divide each of (tau+1) subkmer vectors into nthreads/(tau+1) subvectors
	// as a result, we have subkmer_nthreads threads for sorting
	#pragma omp parallel for shared(vs, vskpq) num_threads( effective_subkmer_threads )
	for (int j=0; j < subkmer_nthreads; ++j) {
		// for each j, we sort subvector (j/(tau+1)) of the vector of subkmers at offset (j%(tau+1))
		boost::function< bool (const hint_t & kmer1, const hint_t & kmer2)  > sub_sort = boost::bind(PositionKMer::compareSubKMers, _1, _2, kmers, tau, PositionKMer::subKMerPositions->at(j % (tau+1)), PositionKMer::subKMerPositions->at((j % (tau+1))+1));
		(*vskpq)[ (j % (tau+1)) ].doSort( j / (tau+1), sub_sort );
	}
}


bool CorrectRead(const vector<KMerCount> & km, hint_t readno, ofstream * ofs) {
	hint_t readno_rev = PositionKMer::revNo + readno;
	const Read & r = PositionKMer::rv->at(readno);
	string seq = r.getSequenceString();
	const uint32_t read_size = PositionKMer::rv->at(readno).size();
	PositionRead & pr = PositionKMer::pr->at(readno);
	PositionRead & pr_rev = PositionKMer::pr->at(readno_rev);

	// create auxiliary structures for consensus
	vector<int> vA(read_size, 0), vC(read_size, 0), vG(read_size, 0), vT(read_size, 0);
	vector< vector<int> > v;  // A=0, C=1, G=2, T=3
	v.push_back(vA); v.push_back(vC); v.push_back(vG); v.push_back(vT);
	PositionKMer::rv_bad->at(readno) = true;

	bool changedRead = false;
	pair<uint32_t, hint_t> it = make_pair( -1, -1 );
	while ( pr.nextKMer( &it ) ) {
		const PositionKMer & kmer = km[it.second].first;
		const uint32_t pos = it.first;
		const KMerStat & stat = km[it.second].second;

		if (stat.changeto == KMERSTAT_GOOD) {
			PositionKMer::rv_bad->at(readno) = false;
			for (size_t j=0; j<K; ++j) {
				v[dignucl(kmer[j])][pos+j]++;
			}
		} else {
			if (stat.changeto < KMERSTAT_CHANGE) {
				PositionKMer::rv_bad->at(readno) = false;
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

	it = make_pair( -1, -1 );
	while ( pr_rev.nextKMer( &it ) ) {
		const PositionKMer & kmer = km[it.second].first;
		const uint32_t pos = it.first;
		const KMerStat & stat = km[it.second].second;

		if (stat.changeto == KMERSTAT_GOOD) {
			PositionKMer::rv_bad->at(readno) = false;
			for (size_t j=0; j<K; ++j) {
				v[complement(dignucl(kmer[j]))][read_size-pos-j-1]++;
			}
		} else {
			if (stat.changeto < KMERSTAT_CHANGE) {
				PositionKMer::rv_bad->at(readno) = false;
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

	/*if (changedRead) {
		cout << "\n" << PositionKMer::rv->at(readno).getSequenceString() << endl;
		for (size_t i=0; i<4; ++i) {
			for (size_t j=0; j<read_size; ++j) {
				cout << (char)((int)'0' + v[i][j]);
			}
			cout << endl;
		}
		cout << seq.data() << endl;
	}*/
	
	PositionKMer::rv->at(readno).setSequence(seq.data());
	return res;
}

struct PriorityQueueElement {
	KMerNo kmerno;
	int n;
	PriorityQueueElement( KMerNo km, int l) : kmerno(km), n(l) { }
};

bool operator < (const PriorityQueueElement & l, const PriorityQueueElement & r) {
	return KMerNo::greater(l.kmerno, r.kmerno);
}

bool operator == (const PriorityQueueElement & l, const PriorityQueueElement & r) {
	return l.kmerno.equal(r.kmerno);
}

void print_time() {
	time_t rawtime;
	tm * ptm;
	time ( &rawtime );
	ptm = gmtime( &rawtime );
	cout << setfill('0') << "[ " << setw(2) << ptm->tm_hour << ":" << setw(2) << ptm->tm_min << ":" << setw(2) << ptm->tm_sec << " ] ";
}

void ParallelSortKMerNos(vector<KMerNo> * v, vector<KMerCount> * kmers, int nthreads) {

	ofstream ofs;

	// find boundaries of the pieces
	vector< size_t > boundaries(nthreads + 1);
	size_t sub_size = (size_t)(v->size() / nthreads);
	for (int j=0; j<nthreads; ++j) {
		boundaries[j] = j * sub_size;
	}
	boundaries[nthreads] = v->size();


	#pragma omp parallel for shared(v, boundaries) num_threads(nthreads)
	for (int j = 0; j < nthreads; ++j) {
		sort(v->begin() + boundaries[j], v->begin() + boundaries[j+1], KMerNo::less);
	}
	TIMEDLN("Subvectors sorted.");

	std::priority_queue< PriorityQueueElement, vector<PriorityQueueElement> > pq;
	vector< vector<KMerNo>::iterator > it(nthreads);
	vector< vector<KMerNo>::iterator > it_end(nthreads);
	for (int j=0; j<nthreads; ++j) {
		it[j] = v->begin() + boundaries[j];
		it_end[j] = v->begin() + boundaries[j+1];
		pq.push( PriorityQueueElement(*(it[j]), j) );
	}

	PriorityQueueElement cur_min = pq.top();
	KMerCount curKMerCount = make_pair( PositionKMer(cur_min.kmerno.index), KMerStat(0, KMERSTAT_GOOD, 1) );
	
	hint_t kmerno = 0;
	double curErrorProb = 1;

	while(pq.size()) {
		const PriorityQueueElement & pqel = pq.top();
		if ( !(cur_min == pqel) ) {
			cur_min = pqel;
			curKMerCount.second.totalQual = 1-curErrorProb;
			kmers->push_back(curKMerCount);
			curKMerCount = make_pair( PositionKMer(cur_min.kmerno.index), KMerStat(0, KMERSTAT_GOOD, 1 ) );
			curErrorProb = 1;
			++kmerno;
		}
		curKMerCount.second.count++;
		curErrorProb *= (1 - PositionKMer::getKMerQuality(pqel.kmerno.index, Globals::qvoffset) );
		PositionKMer::blobkmers[ pqel.kmerno.index ] = kmerno;

		int nn = pqel.n;
		vector<KMerNo>::iterator & it_cur = it[nn];
		pq.pop();
		++it_cur;
		if ( it_cur != it_end[nn] ) pq.push( PriorityQueueElement(*it_cur, nn) );
	}
	curKMerCount.second.totalQual = 1-curErrorProb;
	kmers->push_back(curKMerCount);
}

void outputReads(bool paired, const char * fname, const char * fname_bad, const char * fname_right, const char * fname_right_bad,
			      const char * fname_left_unpaired, const char * fname_right_unpaired) {
	ofstream outf(fname); ofstream outf_bad(fname_bad); 
	if (paired) {
		cout << "outputting paired results" << endl;
		ofstream outf_right(fname_right);
		ofstream outf_right_bad(fname_right_bad);
		ofstream outf_left_unpaired(fname_left_unpaired);
		ofstream outf_right_unpaired(fname_right_unpaired);

		for (hint_t i = 0; i < PositionKMer::lastLeftNo; ++i) {
			if (PositionKMer::rv_bad->at(i)) {
				cout << "  read " << i << " -- left bad " << endl;
				PositionKMer::pr->at(i).print(outf_bad);
			} else {
				if ( i + PositionKMer::lastLeftNo < PositionKMer::revNo && PositionKMer::rv_bad->at(i + PositionKMer::lastLeftNo) ) {
					cout << "  read " << i << " -- left unpaired " << endl;
					PositionKMer::pr->at(i).print(outf_left_unpaired);
				} else {
					cout << "  read " << i << " -- left good " << endl;
					PositionKMer::pr->at(i).print(outf);
				}
			}
		}

		for (hint_t i = PositionKMer::lastLeftNo; i < PositionKMer::revNo; ++i) {
			if (PositionKMer::rv_bad->at(i)) {
					cout << "  read " << i << " -- right bad " << endl;
				PositionKMer::pr->at(i).print(outf_right_bad);
			} else {
				if ( PositionKMer::rv_bad->at(i - PositionKMer::lastLeftNo) ) {
					cout << "  read " << i << " -- right unpaired " << endl;
					PositionKMer::pr->at(i).print(outf_right_unpaired);
				} else {
					cout << "  read " << i << " -- right good " << endl;
					PositionKMer::pr->at(i).print(outf_right);
				}
			}
		}
		outf_right.close(); outf_right_bad.close(); outf_left_unpaired.close(); outf_right_unpaired.close();
	} else {
		for (hint_t i = 0; i < PositionKMer::revNo; ++i) {
			if (PositionKMer::rv_bad->at(i)) {
				PositionKMer::pr->at(i).print(outf_bad);
			} else {
				PositionKMer::pr->at(i).print(outf);
			}
		}
	}

	outf.close(); outf_bad.close();
}

