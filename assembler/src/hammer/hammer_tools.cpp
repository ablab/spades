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
			itf->second.totalQual = itf->second.totalQual * it->second.totalQual;
		} else {
			PositionKMer pkm = it->first;
			KMerStat kms( it->second.count, KMERSTAT_GOODITER, it->second.totalQual );
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
		v->push_back( KMerNo(Globals::pr->at(readno).start() + gen.pos() - 1, 1-gen.correct_probability()) );
		gen.Next();
	}
}

void DoPreprocessing(int tau, string readsFilename, int nthreads, vector<KMerCount*> * kmers, KMerNoHashMap * km) {
	vector< vector<KMerNo> > vtmp;
	for(int n=0; n < nthreads; ++n) {
		vector<KMerNo> v_cur;
		vtmp.push_back(v_cur);
	}

	#pragma omp parallel for shared(vtmp) num_threads(nthreads)
	for(size_t i=0; i < Globals::pr->size(); ++i) {
		if (Globals::pr->at(i).bad()) continue;
		AddKMerNos(Globals::pr->at(i), i, &vtmp[omp_get_thread_num()]);
	}
	
	TIMEDLN("Made KMerNo vector.");
	// fill in the KMerNoHashMap
	for(int n=0; n < nthreads; ++n) {
		for ( vector<KMerNo>::const_iterator it = vtmp[n].begin(); it != vtmp[n].end(); ++it ) {
			KMerNoHashMap::iterator it_hash = km->find( *it );
			if ( it_hash == km->end() ) {
				km->insert( make_pair( *it, new KMerCount( PositionKMer(it->index), KMerStat(1, KMERSTAT_GOODITER, it->errprob) ) ) );
				if (Globals::use_true_likelihood) {
					for (uint32_t j=0; j<K; ++j) {
						Globals::totalquality[it->index + j] = Globals::blobquality[it->index + j] - (char)Globals::qvoffset;
					}
				}
			} else {
				it_hash->second->second.count++;
				it_hash->second->second.totalQual *= it->errprob;
				if (Globals::use_true_likelihood) {
					for (uint32_t j=0; j<K; ++j) {
						Globals::totalquality[it_hash->second->first.start() + j] = (char)min(255,
								(int)Globals::totalquality[it_hash->second->first.start() + j] + (int)Globals::blobquality[it->index + j] - Globals::qvoffset);
					}
				}
			}
		}
	}
	TIMEDLN("Made KMerNo hash map.");
	kmers->clear();
	for ( KMerNoHashMap::const_iterator it_hash = km->begin(); it_hash != km->end(); ++it_hash ) {
		kmers->push_back(it_hash->second);
	}
}

void DoSplitAndSort(int tau, int nthreads, vector< vector<hint_t> > * vs, vector<KMerCount*> * kmers, vector<SubKMerPQ> * vskpq) {
	int effective_threads = min(nthreads, tau+1);
	int subkmer_nthreads = max ( (tau + 1) * ( (int)(nthreads / (tau + 1)) ), tau+1 );
	int effective_subkmer_threads = min(subkmer_nthreads, nthreads);

	#pragma omp parallel for shared(vs, kmers, tau) num_threads(effective_threads)
	for (int j=0; j<tau+1; ++j) {
		vs->at(j).resize( kmers->size() );
		for (size_t m = 0; m < kmers->size(); ++m) vs->at(j)[m] = m;
	}

	for (int j=0; j < tau+1; ++j) {
		SubKMerCompType sort_routine = boost::bind(SubKMerPQElement::compareSubKMerPQElements, _1, _2, kmers, tau, Globals::subKMerPositions->at(j), Globals::subKMerPositions->at(j+1));
		SubKMerPQ skpq( &(vs->at(j)), max( (int)(nthreads / (tau + 1)), 1), sort_routine );
		vskpq->push_back(skpq);
	}

	// we divide each of (tau+1) subkmer vectors into nthreads/(tau+1) subvectors
	// as a result, we have subkmer_nthreads threads for sorting
	#pragma omp parallel for shared(vs, vskpq) num_threads( effective_subkmer_threads )
	for (int j=0; j < subkmer_nthreads; ++j) {
		// for each j, we sort subvector (j/(tau+1)) of the vector of subkmers at offset (j%(tau+1))
		boost::function< bool (const hint_t & kmer1, const hint_t & kmer2)  > sub_sort = boost::bind(PositionKMer::compareSubKMers, _1, _2, kmers, tau, Globals::subKMerPositions->at(j % (tau+1)), Globals::subKMerPositions->at((j % (tau+1))+1));
		(*vskpq)[ (j % (tau+1)) ].doSort( j / (tau+1), sub_sort );
	}
}


size_t CorrectRead(const KMerNoHashMap & hm, const vector<KMerCount*> & km, hint_t readno, ofstream * ofs) {
	hint_t readno_rev = Globals::revNo + readno;
	const Read & r = Globals::rv->at(readno);
	string seq = r.getSequenceString();
	const uint32_t read_size = Globals::rv->at(readno).size();
	PositionRead & pr = Globals::pr->at(readno);
	PositionRead & pr_rev = Globals::pr->at(readno_rev);

	// create auxiliary structures for consensus
	vector<int> vA(read_size, 0), vC(read_size, 0), vG(read_size, 0), vT(read_size, 0);
	vector< vector<int> > v;  // A=0, C=1, G=2, T=3
	v.push_back(vA); v.push_back(vC); v.push_back(vG); v.push_back(vT);
	Globals::rv_bad->at(readno) = true;

	// getting the leftmost and rightmost positions of a solid kmer
	int left = read_size; int right = -1;

	bool changedRead = false;
	pair<int, KMerCount *> it = make_pair( -1, (KMerCount*)NULL );
	while ( (it = pr.nextKMer(it.first)).first > -1 ) {
		const PositionKMer & kmer = it.second->first;
		const uint32_t pos = it.first;
		const KMerStat & stat = it.second->second;

		if (  stat.isGoodForIterative() || 
			// if regular_threshold_for_correction = true, we use a (more relaxed) threshold isGood() for solid k-mers
		      ( (!Globals::use_iterative_reconstruction || Globals::regular_threshold_for_correction) && stat.isGood() ) ) {
			Globals::rv_bad->at(readno) = false;
			for (size_t j=0; j<K; ++j) {
				v[dignucl(kmer[j])][pos+j]++;
			}
			if ((int)pos < left) left = pos; if ((int)pos > right) right = pos;
		} else {
			// if discard_only_singletons = true, we always use centers of clusters that do not coincide with the current center
			if (stat.change() && ( Globals::discard_only_singletons || km[stat.changeto]->second.isGoodForIterative() || 
			((!Globals::use_iterative_reconstruction || Globals::regular_threshold_for_correction) && km[stat.changeto]->second.isGood()) ) ) {
				Globals::rv_bad->at(readno) = false;
				if ((int)pos < left) left = pos; if ((int)pos > right) right = pos;
				const PositionKMer & newkmer = km[ stat.changeto ]->first;
				
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

	it = make_pair( -1, (KMerCount*)NULL );
	while ( (it = pr_rev.nextKMer(it.first)).first > -1 ) {
		const PositionKMer & kmer = it.second->first;
		const uint32_t pos = it.first;
		const KMerStat & stat = it.second->second;

		if (stat.isGoodForIterative() || ( !Globals::use_iterative_reconstruction && stat.isGood() ) ) {
			Globals::rv_bad->at(readno) = false;
			for (size_t j=0; j<K; ++j) {
				v[complement(dignucl(kmer[j]))][read_size-pos-j-1]++;
			}
			if ((int)(read_size-K-pos) < left) left = (int)(read_size-K-pos); if ((int)(read_size-K-pos) > right) right = (int)(read_size-K-pos);
		} else {
			if (stat.change() && (km[stat.changeto]->second.isGoodForIterative() || 
						( !Globals::use_iterative_reconstruction && km[stat.changeto]->second.isGood() ) ) ) {
				Globals::rv_bad->at(readno) = false;
				if ((int)(read_size-K-pos) < left) left = (int)(read_size-K-pos); if ((int)(read_size-K-pos) > right) right = (int)(read_size-K-pos);
				const PositionKMer & newkmer = km[ stat.changeto ]->first;

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

	size_t res = 0; // how many nucleotides have really changed?
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
		cout << "\n" << Globals::rv->at(readno).getSequenceString() << endl;
		for (size_t i=0; i<4; ++i) {
			for (size_t j=0; j<read_size; ++j) {
				cout << (char)((int)'0' + v[i][j]);
			}
			cout << endl;
		}
		cout << seq.data() << endl;
	}*/
	
	Globals::rv->at(readno).setSequence(seq.data());
	if (ofs != NULL && changedRead) {
		*ofs << "Final result:  size=" << Globals::rv->at(readno).size() << "\n" << Globals::rv->at(readno).getSequenceString().c_str() << "\n" << Globals::rv->at(readno).getPhredQualityString(Globals::qvoffset).c_str() << endl;
	}
	if (Globals::trim_left_right) {
		Globals::rv->at(readno).trimLeftRight(left, right+K-1);
		if (ofs != NULL && changedRead) {
			*ofs << "Trimming to [ " << left << ", " << right+K-1 << "]" << endl;
			*ofs << "Trimmed: " << Globals::rv->at(readno).getSequenceString().c_str() << endl;
		}
	}
	return res;
}

void print_time() {
	time_t rawtime;
	tm * ptm;
	time ( &rawtime );
	ptm = gmtime( &rawtime );
	cout << setfill('0') << "[ " << setw(2) << ptm->tm_hour << ":" << setw(2) << ptm->tm_min << ":" << setw(2) << ptm->tm_sec << " ] ";
}

struct PriorityQueueElement {
        KMerNo kmerno;
        int n;
        PriorityQueueElement( KMerNo km, int l) : kmerno(km), n(l) { }
};

bool operator < (const PriorityQueueElement & l, const PriorityQueueElement & r) {
        return l.kmerno.greater(r.kmerno);
}

bool operator == (const PriorityQueueElement & l, const PriorityQueueElement & r) {
        return l.kmerno.equal(r.kmerno);
}

void outputReads(bool paired, const char * fname, const char * fname_bad, const char * fname_right, const char * fname_right_bad,
			      const char * fname_left_unpaired, const char * fname_right_unpaired) {
	ofstream outf(fname); ofstream outf_bad(fname_bad); 
	if (paired) {
		ofstream outf_right(fname_right);
		ofstream outf_right_bad(fname_right_bad);
		ofstream outf_left_unpaired(fname_left_unpaired);
		ofstream outf_right_unpaired(fname_right_unpaired);

		for (hint_t i = 0; i < Globals::lastLeftNo; ++i) {
			if (Globals::rv_bad->at(i)) {
				Globals::pr->at(i).print(outf_bad);
			} else {
				if ( i + Globals::lastLeftNo < Globals::revNo && Globals::rv_bad->at(i + Globals::lastLeftNo) ) {
					Globals::pr->at(i).print(outf_left_unpaired);
				} else {
					Globals::pr->at(i).print(outf);
				}
			}
		}

		for (hint_t i = Globals::lastLeftNo; i < Globals::revNo; ++i) {
			if (Globals::rv_bad->at(i)) {
				Globals::pr->at(i).print(outf_right_bad);
			} else {
				if ( Globals::rv_bad->at(i - Globals::lastLeftNo) ) {
					Globals::pr->at(i).print(outf_right_unpaired);
				} else {
					Globals::pr->at(i).print(outf_right);
				}
			}
		}
		outf_right.close(); outf_right_bad.close(); outf_left_unpaired.close(); outf_right_unpaired.close();
	} else {
		for (hint_t i = 0; i < Globals::revNo; ++i) {
			if (Globals::rv_bad->at(i)) {
				Globals::pr->at(i).print(outf_bad);
			} else {
				Globals::pr->at(i).print(outf);
			}
		}
	}

	outf.close(); outf_bad.close();
}

size_t IterativeReconstructionStep(int nthreads, const vector<KMerCount*> & kmers, ostream * ofs) {
	size_t res = 0;
	// cycle over the reads, looking for reads completely covered by solid k-mers and adding new solid k-mers on the fly
	#pragma omp parallel for shared(res, ofs) num_threads(nthreads)
	for (hint_t readno = 0; readno < Globals::revNo; ++readno) {
		if ( ofs != NULL ) {
			#pragma omp critical
			{
			(*ofs) << "Read: " << Globals::rv->at(readno).getSequenceString().c_str() << "\n";
			}
		}
		const uint32_t read_size = Globals::rv->at(readno).size();
		vector< bool > covered_by_solid( read_size, false );

		const PositionRead & pr = Globals::pr->at(readno);
		pair<int, KMerCount *> it = make_pair( -1, (KMerCount*)NULL );
		while ( (it = pr.nextKMer(it.first)).first > -1 ) {
			if ( it.second->second.isGoodForIterative() ) {
				for ( size_t j = it.first; j < it.first + K; ++j )
					covered_by_solid[j] = true;
			}
		}
		bool isGood = true;
		for ( size_t j = 0; j < read_size; ++j ) {
			if ( !covered_by_solid[j] ) { isGood = false; break; }
		}
		if ( !isGood ) continue;

		const PositionRead & pr_rev = Globals::pr->at(Globals::revNo + readno);
		std::fill( covered_by_solid.begin(), covered_by_solid.end(), false );
		it = make_pair( -1, (KMerCount*)NULL );
		while ( (it = pr_rev.nextKMer(it.first)).first > -1 ) {
			if ( it.second->second.isGoodForIterative() ) {
				for ( size_t j = read_size-it.first-K; j < read_size-it.first; ++j )
					covered_by_solid[j] = true;
			}
		}
		isGood = true;
		for ( size_t j = 0; j < read_size; ++j ) {
			if ( !covered_by_solid[j] ) { isGood = false; break; }
		}
		if ( !isGood ) continue;

		// ok, now we're sure that everything is covered
		// let's mark all k-mers as solid
		it = make_pair( -1, (KMerCount*)NULL );
		while ( (it = pr.nextKMer(it.first)).first > -1 ) {
			if ( !it.second->second.isGoodForIterative() && !it.second->second.isMarkedGoodForIterative() ) {
				#pragma omp critical
				{
				++res;
				if ( ofs != NULL ) (*ofs) << "    make solid: " << it.second->first.str().c_str() << "\n";
				if (Globals::reconstruction_in_full_iterations) it.second->second.markGoodForIterative();
				else it.second->second.makeGoodForIterative();
				}
			}
		}
		it = make_pair( -1, (KMerCount*)NULL );
		while ( (it = pr_rev.nextKMer(it.first)).first > -1 ) {
			if ( !it.second->second.isGoodForIterative() && !it.second->second.isMarkedGoodForIterative() ) {
				#pragma omp critical
				{
				++res;
				if ( ofs != NULL ) (*ofs) << "    make solid: " << it.second->first.str().c_str() << "\n";
				if (Globals::reconstruction_in_full_iterations) it.second->second.markGoodForIterative();
				else it.second->second.makeGoodForIterative();
				}
			}				
		}
	}

	if (Globals::reconstruction_in_full_iterations) {
		// ok, so now we've marked everything and simply need to go through marked k-mers actually making them good
		for ( hint_t i=0; i < kmers.size(); ++i ) {
			if ( kmers[i]->second.isMarkedGoodForIterative() ) {
				kmers[i]->second.makeGoodForIterative();
			}
		}
	}

	return res;
}


