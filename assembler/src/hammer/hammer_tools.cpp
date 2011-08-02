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

void DoPreprocessing(int tau, int qvoffset, string readsFilename, int nthreads, vector<KMerNo> * vv) {
	vv->clear();

	// TODO: think about a parallelization -- for some reason, the previous version started producing segfaults

	vector< vector<KMerNo> > vtmp;
	for(size_t n=0; n < nthreads; ++n) {
		vector<KMerNo> v_cur;
		vtmp.push_back(v_cur);
	}

	#pragma omp parallel for shared(vtmp) num_threads(nthreads)
	for(size_t i=0; i < PositionKMer::pr->size(); ++i) {
		AddKMerNos(PositionKMer::pr->at(i), i, &vtmp[omp_get_thread_num()]);
		//if ( i % 1000000 == 0 ) cout << "Read no. " << i << " processed by thread " << omp_get_thread_num() << "." << endl;
	}
	//TIMEDLN("All k-mers added to vectors.");
	
	for(size_t n=0; n < nthreads; ++n) {
		vv->insert(vv->end(), vtmp[n].begin(), vtmp[n].end());
		vtmp[n].clear();
	}
	//cout << "Vectors of k-mers joined. Got a vector of size " << vv->size() << "." << endl;
}

void DoSplitAndSort(int tau, int nthreads, const vector<KMerNo> & vv, vector< vector<hint_t> > * vs, vector<KMerCount> * kmers, vector<SubKMerPQ> * vskpq) {
	int effective_threads = min(nthreads, tau+1);
	int subkmer_nthreads = max ( (tau + 1) * ( (int)(nthreads / (tau + 1)) ), tau+1 );
	int effective_subkmer_threads = min(subkmer_nthreads, nthreads);

	//cout << "Starting split and sort..." << endl;

	#pragma omp parallel for shared(vs, kmers, tau) num_threads(effective_threads)
	for (int j=0; j<tau+1; ++j) {
		vs->at(j).resize( kmers->size() );
		for (size_t m = 0; m < kmers->size(); ++m) vs->at(j)[m] = m;
	}

	for (size_t j=0; j < tau+1; ++j) {
		subkmer_comp_type sort_routine = boost::bind(SubKMerPQElement::compareSubKMerPQElements, _1, _2, kmers, tau, PositionKMer::subKMerPositions->at(j), PositionKMer::subKMerPositions->at(j+1));
		SubKMerPQ skpq( &(vs->at(j)), max( (int)(nthreads / (tau + 1)), 1), sort_routine );
		vskpq->push_back(skpq);
	}

	// we divide each of (tau+1) subkmer vectors into nthreads/(tau+1) subvectors
	// as a result, we have subkmer_nthreads threads for sorting
	#pragma omp parallel for shared(vs, vskpq) num_threads( effective_subkmer_threads )
	for (size_t j=0; j < subkmer_nthreads; ++j) {
		// for each j, we sort subvector (j/(tau+1)) of the vector of subkmers at offset (j%(tau+1))
		boost::function< bool (const hint_t & kmer1, const hint_t & kmer2)  > sub_sort = boost::bind(PositionKMer::compareSubKMers, _1, _2, kmers, tau, PositionKMer::subKMerPositions->at(j % (tau+1)), PositionKMer::subKMerPositions->at((j % (tau+1))+1));
		(*vskpq)[ (j % (tau+1)) ].doSort( j / (tau+1), sub_sort );
	}
	//cout << "Auxiliary subvectors sorted and sent to priority queues." << endl;
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

	bool changedRead = false;
	pair<uint32_t, hint_t> it = make_pair( -1, -1 );
	while ( pr.nextKMer( &it ) ) {
	//for (map<uint32_t, hint_t>::const_iterator it = pr.kmers().begin(); it != pr.kmers().end(); ++it) {
		const PositionKMer & kmer = km[it.second].first;
		const uint32_t pos = it.first;
		const KMerStat & stat = km[it.second].second;

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

	it = make_pair( -1, -1 );
	while ( pr_rev.nextKMer( &it ) ) {
//	for (map<uint32_t, hint_t>::const_iterator it = pr_rev.kmers().begin(); it != pr_rev.kmers().end(); ++it) {
		const PositionKMer & kmer = km[it.second].first;
		const uint32_t pos = it.first;
		const KMerStat & stat = km[it.second].second;
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
	
	// fill newblob
	//for (size_t j=0; j < pr.size(); ++j) {
	//	newblob[ pr.start() + j ] = seq[j];
	//}
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

	// find boundaries of the pieces
	vector< size_t > boundaries(nthreads + 1);
	size_t sub_size = (size_t)(v->size() / nthreads);
	for (size_t j=0; j<nthreads; ++j) {
		boundaries[j] = j * sub_size;
	}
	boundaries[nthreads] = v->size();
	//cout << "  thread boundaries: "; for (size_t j=0; j<nthreads+1; ++j) cout << boundaries[j] << " "; cout << endl;

	#pragma omp parallel for shared(v, boundaries) num_threads(nthreads)
	for (int j = 0; j < nthreads; ++j) {
		sort(v->begin() + boundaries[j], v->begin() + boundaries[j+1], KMerNo::less);
	}
	TIMEDLN("Subvectors sorted.");

	std::priority_queue< PriorityQueueElement, vector<PriorityQueueElement> > pq;
	vector< vector<KMerNo>::iterator > it(nthreads);
	vector< vector<KMerNo>::iterator > it_end(nthreads);
	for (size_t j=0; j<nthreads; ++j) {
		it[j] = v->begin() + boundaries[j];
		it_end[j] = v->begin() + boundaries[j+1];
		pq.push( PriorityQueueElement(*(it[j]), j) );
	}

	PriorityQueueElement cur_min = pq.top();
	KMerCount curKMerCount = make_pair( PositionKMer(cur_min.kmerno.index), KMerStat(0, KMERSTAT_GOOD) );
	
	hint_t kmerno = 0;

	while(pq.size()) {
		PriorityQueueElement pqel = pq.top(); pq.pop();
		// cout << "cur_min = " << cur_min.kmerno.index << "   popped " << pqel.kmerno.index << endl;
		if ( !(cur_min == pqel) ) {
			//cout << "push " << curKMerCount.first.str() << " cnt=" << curKMerCount.second.count << endl;
			cur_min = pqel;
			kmers->push_back(curKMerCount);
			curKMerCount = make_pair( PositionKMer(cur_min.kmerno.index), KMerStat(0, KMERSTAT_GOOD) );
			++kmerno;			
		}
		curKMerCount.second.count++;
		PositionKMer::blobkmers[ cur_min.kmerno.index ] = kmerno;

		++it[pqel.n];
		if ( it[pqel.n] != it_end[pqel.n] ) pq.push( PriorityQueueElement(*(it[pqel.n]), pqel.n) );
	}
	//cout << "push " << curKMerCount.first.str() << " cnt=" << curKMerCount.second.count << endl;
	kmers->push_back(curKMerCount);

	//cout << "Subvectors merged. In total, we have " << kmers->size() << " kmers." << endl;
}


