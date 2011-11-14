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
#include <boost/format.hpp>
#include <time.h>
#include <sys/resource.h>
#include <iomanip>

#include "read/ireadstream.hpp"
#include "defs.hpp"
#include "mathfunctions.hpp"
#include "valid_kmer_generator.hpp"
#include "position_kmer.hpp"
#include "globals.hpp"
#include "config_struct_hammer.hpp"
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
				KMerCount * kmc = new KMerCount( PositionKMer(it->index), KMerStat(1, KMERSTAT_GOODITER, it->errprob) );
				if (Globals::use_true_likelihood) {
					for (uint32_t j=0; j<K; ++j) {
						kmc->second.qual.set(j, Globals::blobquality[it->index + j] - (char)Globals::qvoffset);
					}
				}
				km->insert( make_pair( *it, kmc ) );
			} else {
				it_hash->second->second.count++;
				it_hash->second->second.totalQual *= it->errprob;
				if (Globals::use_true_likelihood) {
					for (uint32_t j=0; j<K; ++j) {
						const int cur_qual = (int)it_hash->second->second.qual[j];
						if (cur_qual + (int)Globals::blobquality[it->index + j] - (int)Globals::qvoffset < MAX_SHORT) {
							it_hash->second->second.qual.set(j, cur_qual + (short)Globals::blobquality[it->index + j] - (short)Globals::qvoffset);
						} else {
							it_hash->second->second.qual.set(j, MAX_SHORT);
						}
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

bool internalCorrectReadProcedure( const Read & r, const hint_t readno, const string & seq, const vector<KMerCount*> & km,
		const PositionKMer & kmer, const uint32_t pos, const KMerStat & stat, vector< vector<int> > & v,
		int & left, int & right, bool & isGood, ofstream * ofs ) {
	bool res = false;
	//if (ofs != NULL)
	//	*ofs << "\n " << r.getName() << "\n" << seq.data() << "\n";
	if (  stat.isGoodForIterative() ||
			// if regular_threshold_for_correction = true, we use a (more relaxed) threshold isGood() for solid k-mers
			((!Globals::use_iterative_reconstruction
					|| Globals::regular_threshold_for_correction)
					&& stat.isGood())) {
		isGood = true;
		if (ofs != NULL) *ofs << "\t\t\tsolid";
		for (size_t j = 0; j < K; ++j) {
			v[dignucl(kmer[j])][pos + j]++;
		}
		if ((int) pos < left)
			left = pos;
		if ((int) pos > right)
			right = pos;
		/*if (ofs != NULL) {
				for (size_t j = 0; j < pos; ++j)
					*ofs << " ";
				*ofs << kmer.str().data() << "\n";
		}*/
	} else {
		// if discard_only_singletons = true, we always use centers of clusters that do not coincide with the current center
		if (stat.change() && (Globals::discard_only_singletons
				|| km[stat.changeto]->second.isGoodForIterative()
				|| ((!Globals::use_iterative_reconstruction
						|| Globals::regular_threshold_for_correction)
						&& km[stat.changeto]->second.isGood()))) {
			if (ofs != NULL) *ofs << "\tchange to\n";
			//cout << "  kmer " << kmer.str() << " wants to change to " << km[stat.changeto]->first.str() << endl;
			isGood = true;
			if ((int) pos < left)
				left = pos;
			if ((int) pos > right)
				right = pos;
			const PositionKMer & newkmer = km[stat.changeto]->first;

			for (size_t j = 0; j < K; ++j) {
				v[dignucl(newkmer[j])][pos + j]++;
			}
			// pretty print the k-mer
			res = true;
			if (ofs != NULL) {
				for (size_t j = 0; j < pos; ++j)
					*ofs << " ";
				*ofs << newkmer.str().data();
			}
		}
	}
	return res;
}

size_t CorrectRead(const KMerNoHashMap & hm, const vector<KMerCount*> & km, hint_t readno, Read & r, bool & isGood, ofstream * ofs) {
	hint_t readno_rev = Globals::revNo + readno;
	string seq (Globals::blob        + Globals::pr->at(readno).start(), Globals::pr->at(readno).size());
	string qual(Globals::blobquality + Globals::pr->at(readno).start(), Globals::pr->at(readno).size());
	const uint32_t read_size = Globals::pr->at(readno).size();
	PositionRead & pr = Globals::pr->at(readno);
	PositionRead & pr_rev = Globals::pr->at(readno_rev);

	/*if (ofs != NULL) {
	*ofs << "    readno=" << readno << "  size=" << read_size << " kmersize=" << km.size() << " pr=" << pr.start() << endl;
	*ofs << "    " << seq << endl;
	*ofs << "    " << qual << endl;

	*ofs << "    readno_rev=" << readno_rev << "  size=" << read_size << " kmersize=" << km.size() << " pr=" << pr.start() << endl;
	*ofs << "    " << string(Globals::blob        + pr_rev.start(), pr_rev.size()) << endl;
	*ofs << "    " << string(Globals::blobquality + pr_rev.start(), pr_rev.size()) << endl;
	}*/

	// create auxiliary structures for consensus
	vector<int> vA(read_size, 0), vC(read_size, 0), vG(read_size, 0), vT(read_size, 0);
	vector< vector<int> > v;  // A=0, C=1, G=2, T=3
	v.push_back(vA); v.push_back(vC); v.push_back(vG); v.push_back(vT);
	isGood = false;

	// getting the leftmost and rightmost positions of a solid kmer
	int left = read_size; int right = -1;

	if (ofs != NULL)
		*ofs << "\n " << r.getName() << "\n" << seq.data() << "\n";

	// cout << seq << endl;
	bool changedRead = false;
	if (Globals::conserve_memory) {
		pair<int, hint_t> it = make_pair( -1, BLOBKMER_UNDEFINED );
		while ( (it = pr.nextKMerNo(it.first)).first > -1 ) {
			const PositionKMer & kmer = km[it.second]->first;
			const uint32_t pos = it.first;
			const KMerStat & stat = km[it.second]->second;

			if (ofs != NULL) {
				for (uint32_t i=0; i<pos; ++i) *ofs << " ";
				*ofs << kmer.str();
			}

			changedRead = changedRead || internalCorrectReadProcedure( r, readno, seq, km, kmer, pos, stat, v, left, right, isGood, ofs );

			if (ofs != NULL) *ofs << "\n";
		}
	} else {
		pair<int, KMerCount *> it = make_pair( -1, (KMerCount*)NULL );
		while ( (it = pr.nextKMer(it.first)).first > -1 ) {
			const PositionKMer & kmer = it.second->first;
			const uint32_t pos = it.first;
			const KMerStat & stat = it.second->second;

			changedRead = changedRead || internalCorrectReadProcedure( r, readno, seq, km, kmer, pos, stat, v, left, right, isGood, ofs );
		}
	}


	if (ofs != NULL)
		*ofs << "\n " << r.getName() << " revcomp\n" << string(Globals::blob        + pr_rev.start(), pr_rev.size()) << "\n";

	int left_rev = read_size; int right_rev = -1;

	if (Globals::conserve_memory) {
		pair<int, hint_t> it = make_pair( -1, BLOBKMER_UNDEFINED );
		while ( (it = pr_rev.nextKMerNo(it.first)).first > -1 ) {
			const PositionKMer & kmer = km[it.second]->first;
			const uint32_t pos = it.first;
			const KMerStat & stat = km[it.second]->second;

			if (ofs != NULL) {
				for (uint32_t i=0; i<pos; ++i) *ofs << " ";
				*ofs << kmer.str();
			}

			changedRead = changedRead || internalCorrectReadProcedure( r, readno, seq, km, kmer, pos, stat, v, left_rev, right_rev, isGood, ofs );

			if (ofs != NULL) *ofs << "\n";
		}
	} else {
		pair<int, KMerCount *> it = make_pair( -1, (KMerCount*)NULL );
		while ( (it = pr_rev.nextKMer(it.first)).first > -1 ) {
			const PositionKMer & kmer = it.second->first;
			const uint32_t pos = it.first;
			const KMerStat & stat = it.second->second;

			changedRead = changedRead || internalCorrectReadProcedure( r, readno, seq, km, kmer, pos, stat, v, left_rev, right_rev, isGood, ofs );
		}
	}

	if (ofs != NULL) {
		*ofs << "    direct=[" << left << ", " << right << "]\treverse=[" << left_rev << ", " << right_rev << "]";
	}

	if ( left <= right && left_rev <= right_rev ) {
		left = std::min(left, (int)read_size - left_rev - (int)K);
		right = std::max(right, (int)read_size - right_rev - (int)K);
	} else if ( left > right && left_rev <= right_rev ) {
		left = (int)read_size - left_rev - (int)K;
		right = (int)read_size - right_rev - (int)K;
	}
	if (ofs != NULL) {
		*ofs << "\tresult=[" << left << ", " << right << "]" << endl;
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
		if (seq[j] != cmax) ++res;
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

	r.setSequence(seq.data());
	if (Globals::trim_left_right) {
		r.trimLeftRight(left, right+K-1);
		if ( left > 0 || right + K -1 < read_size ) changedRead = true;
		if (ofs != NULL && changedRead) {
			*ofs << "Trimming to [ " << left << ", " << right+K-1 << "]" << endl;
			// *ofs << "Trimmed: " << r.getSequenceString().c_str() << endl;
		}
	}

	if (ofs != NULL && changedRead) {
		*ofs << "Final result:  size=" << r.size() << "\n" << r.getSequenceString() << "\n" << r.getPhredQualityString(Globals::qvoffset) << endl;
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

void print_mem_usage() {
	static size_t pid = getpid();
	std::string str = (boost::format("pmap -d %d | grep writeable/private") % pid).str();
	std::cout << "==== MEM USAGE ==== " << std::endl;
	system(str.c_str());
}

void print_stats() {
	rusage ru;
	getrusage(RUSAGE_SELF, &ru);
	std::cout << "[";
	tm * ptm = gmtime( &ru.ru_utime.tv_sec );
	std::cout << " elapsed = " << setw(2) << setfill('0')  << ptm->tm_hour << ":" << setw(2) << setfill('0') << ptm->tm_min << ":" << setw(2) << setfill('0') << ptm->tm_sec;
	if ( ru.ru_maxrss < 1024 * 1024 ) {
		std::cout << " mem = " << setw(5) << setfill(' ') << (ru.ru_maxrss / 1024) << "M ";
	} else {
		std::cout << " mem = " << setw(5) << setprecision(1) << fixed << setfill(' ') << (ru.ru_maxrss / (1024.0 * 1024.0) ) << "G ";
	}
	std::cout << "] ";
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
		const uint32_t read_size = Globals::pr->at(readno).size();
		string seq = Globals::conserve_memory
				? string(Globals::blob + Globals::pr->at(readno).start(), read_size)
				: Globals::rv->at(readno).getSequenceString();
		if ( ofs != NULL ) {
			#pragma omp critical
			{
			(*ofs) << "Read: " << seq.c_str() << "\n";
			}
		}
		vector< bool > covered_by_solid( read_size, false );

		const PositionRead & pr = Globals::pr->at(readno);
		if (!Globals::conserve_memory) {
			pair<int, KMerCount *> it = make_pair( -1, (KMerCount*)NULL );
			while ( (it = pr.nextKMer(it.first)).first > -1 ) {
				if ( it.second->second.isGoodForIterative() ) {
					for ( size_t j = it.first; j < it.first + K; ++j )
						covered_by_solid[j] = true;
				}
			}
		} else {
			pair<int, hint_t > it = make_pair( -1, BLOBKMER_UNDEFINED );
			while ( (it = pr.nextKMerNo(it.first)).first > -1 ) {
				if ( kmers[it.second]->second.isGoodForIterative() ) {
					for ( size_t j = it.first; j < it.first + K; ++j )
						covered_by_solid[j] = true;
				}
			}
		}
		bool isGood = true;
		for ( size_t j = 0; j < read_size; ++j ) {
			if ( !covered_by_solid[j] ) { isGood = false; break; }
		}
		if ( !isGood ) continue;

		const PositionRead & pr_rev = Globals::pr->at(Globals::revNo + readno);
		string seq_rev = Globals::conserve_memory ?
						string(Globals::blob	+ Globals::pr->at(Globals::revNo + readno).start(),	read_size)
						: Globals::rv->at(Globals::revNo + readno).getSequenceString();
		std::fill(covered_by_solid.begin(), covered_by_solid.end(), false);
		if (!Globals::conserve_memory) {
			pair<int, KMerCount *> it = make_pair(-1, (KMerCount*) NULL );
			while ((it = pr_rev.nextKMer(it.first)).first > -1) {
				if (it.second->second.isGoodForIterative()) {
					for (size_t j = it.first; j < it.first + K; ++j) covered_by_solid[j] = true;
				}
			}
		} else {
			pair<int, hint_t> it = make_pair(-1, BLOBKMER_UNDEFINED );
			while ((it = pr_rev.nextKMerNo(it.first)).first > -1) {
				if (kmers[it.second]->second.isGoodForIterative()) {
					for (size_t j = it.first; j < it.first + K; ++j) covered_by_solid[j] = true;
				}
			}
		}
		isGood = true;
		for (size_t j = 0; j < read_size; ++j) {
			if (!covered_by_solid[j]) {
				isGood = false;
				break;
			}
		}
		if (!isGood) continue;

		// ok, now we're sure that everything is covered
		// let's mark all k-mers as solid
		if (!Globals::conserve_memory) {
			pair<int, KMerCount *> it = make_pair( -1, (KMerCount*)NULL );
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
		} else {
			pair<int, hint_t > it = make_pair( -1, BLOBKMER_UNDEFINED );
			while ( (it = pr.nextKMerNo(it.first)).first > -1 ) {
				if ( !kmers[it.second]->second.isGoodForIterative() && !kmers[it.second]->second.isMarkedGoodForIterative() ) {
					#pragma omp critical
					{
					++res;
					if ( ofs != NULL ) (*ofs) << "    make solid: " << kmers[it.second]->first.str().c_str() << "\n";
					if (Globals::reconstruction_in_full_iterations) kmers[it.second]->second.markGoodForIterative();
					else kmers[it.second]->second.makeGoodForIterative();
					}
				}
			}
			it = make_pair( -1, BLOBKMER_UNDEFINED );
			while ( (it = pr_rev.nextKMerNo(it.first)).first > -1 ) {
				if ( !kmers[it.second]->second.isGoodForIterative() && !kmers[it.second]->second.isMarkedGoodForIterative() ) {
					#pragma omp critical
					{
					++res;
					if ( ofs != NULL ) (*ofs) << "    make solid: " << kmers[it.second]->first.str().c_str() << "\n";
					if (Globals::reconstruction_in_full_iterations) kmers[it.second]->second.markGoodForIterative();
					else kmers[it.second]->second.makeGoodForIterative();
					}
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

void SplitToFiles(string dirprefix, int iter_count) {
	vector<ofstream *> ofiles(Globals::num_of_tmp_files);
	for (int i = 0; i < Globals::num_of_tmp_files; ++i) {
		ofiles[i] = new ofstream( getFilename( dirprefix, iter_count, "tmp.kmers", i ) );
	}
	Seq<K>::hash hash_function;
	for(size_t i=0; i < Globals::pr->size(); ++i) {
		string s(Globals::blob        + Globals::pr->at(i).start(), Globals::pr->at(i).size());
		string q(Globals::blobquality + Globals::pr->at(i).start(), Globals::pr->at(i).size());
		for ( size_t j=0; j < Globals::pr->at(i).size(); ++j) q[j] = (char)(q[j] - Globals::qvoffset);
		ValidKMerGenerator<K> gen(s, q);
		while (gen.HasMore()) {
			ofstream &cur_file = *ofiles[hash_function(gen.kmer()) % Globals::num_of_tmp_files];
			hint_t cur_pos = Globals::pr->at(i).start() + gen.pos() - 1;
			double correct_probability = 1 - gen.correct_probability();
			cur_file << cur_pos << "\t" << correct_probability << "\n";
			gen.Next();
		}
	}
	for (int i = 0; i < Globals::num_of_tmp_files; ++i) {
		ofiles[i]->close();
		delete ofiles[i];
	}
}

void ProcessKmerHashFile( ifstream * inf, KMerNoHashMap & km ) {
	char buf[1024]; // a line contains two numbers, 1024 should be enough for everybody
	uint64_t pos; double prob;
	while (!inf->eof()) {
		inf->getline(buf, 1024);
		sscanf(buf, "%lu\t%lf", &pos, &prob);
		KMerNo kmerno(pos, prob);
		KMerNoHashMap::iterator it_hash = km.find( kmerno );
		if ( it_hash == km.end() ) {
			KMerCount * kmc = new KMerCount( PositionKMer(kmerno.index), KMerStat(1, KMERSTAT_GOODITER, kmerno.errprob) );
			for (uint32_t j=0; j<K; ++j) {
				kmc->second.qual.set(j, Globals::blobquality[kmerno.index + j] - (char)Globals::qvoffset);
			}
			km.insert( make_pair( kmerno, kmc ) );
		} else {
			it_hash->second->second.count++;
			it_hash->second->second.totalQual *= kmerno.errprob;
			for (uint32_t j=0; j<K; ++j) {
				const int cur_qual = (int)it_hash->second->second.qual[j];
				if (cur_qual + (int)Globals::blobquality[kmerno.index + j] - (int)Globals::qvoffset < MAX_SHORT) {
					it_hash->second->second.qual.set(j, cur_qual + (short)Globals::blobquality[kmerno.index + j] - (short)Globals::qvoffset);
				} else {
					it_hash->second->second.qual.set(j, MAX_SHORT);
				}
			}
		}
	}
}

void PrintProcessedKmerHashFile( ofstream * outf, hint_t & kmer_num, KMerNoHashMap & km ) {
	for (KMerNoHashMap::iterator it = km.begin(); it != km.end(); ++it) {
		(*outf) << it->second->first.start() << "\t"
				<< string(Globals::blob + it->second->first.start(), K) << "\t"
				<< it->second->second.count << "\t"
				<< setw(8) << it->second->second.totalQual << "\t";
		for (size_t i=0; i < K; ++i) (*outf) << it->second->second.qual[i] << " ";
		(*outf) << "\n";
		delete it->second;
		++kmer_num;
	}
	km.clear();
}

void PrintKMerFileWithChangeTo( ofstream * outf, const vector<KMerCount *> & kmers ) {
	for (vector<KMerCount *>::const_iterator it = kmers.begin(); it != kmers.end(); ++it) {
		(*outf) << (*it)->first.start() << "\t"
				<< string(Globals::blob + (*it)->first.start(), K) << "\t"
				<< (*it)->second.count << "\t"
				<< (*it)->second.changeto << "\t"
				<< setw(8) << (*it)->second.totalQual << "\t";
		for (size_t i=0; i < K; ++i) (*outf) << (*it)->second.qual[i] << " ";
		(*outf) << "\n";
	}

}

string getFilename( const string & dirprefix, const string & suffix ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << suffix.data();
	return tmp.str();
}

string getFilename( const string & dirprefix, int iter_count, const string & suffix ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << "." << suffix.data();
	return tmp.str();
}

string getFilename( const string & dirprefix, int iter_count, const string & suffix, int suffix_num ) {
	ostringstream tmp;
	tmp.str(""); tmp << dirprefix.data() << "/" << std::setfill('0') << std::setw(2) << iter_count << "." << suffix.data() << "." << suffix_num;
	return tmp.str();
}

void fillInKmersFromFile( const string & fname, vector<hint_t> *kmernos ) {
	kmernos->clear();
	ifstream ifs(fname);
	char buf[16000];
	hint_t pos;
	while (!ifs.eof()) {
		ifs.getline(buf, 16000);
		sscanf(buf, "%lu", &pos);
		kmernos->push_back(pos);
	}
	ifs.close();
}

void fillInKmersAndNosFromFile( const string & fname, vector<KMerCount*> *kmers, vector<hint_t> *kmernos ) {
	kmernos->clear();
	kmers->clear();
	ifstream ifs(fname);
	char buf[16000];
	while (!ifs.eof()) {
		ifs.getline(buf, 16000);
		hint_t pos; int cnt; double qual;
		sscanf(buf, "%lu\t%*s\t%u\t%lf", &pos, &cnt, &qual);
		//cout << pos << "\t" << cnt << "\t" << qual << endl;
		kmernos->push_back(pos);
		kmers->push_back(new KMerCount(PositionKMer(pos), KMerStat(cnt, KMERSTAT_BAD, qual)));
	}
	ifs.close();
}

void fillInKmersWithChangeToFromFile( const string & fname, vector<KMerCount*> *kmers, vector<hint_t> *kmernos ) {
	kmernos->clear();
	kmers->clear();
	ifstream ifs(fname);
	char buf[16000];
	while (!ifs.eof()) {
		ifs.getline(buf, 16000);
		hint_t pos; int cnt; double qual; hint_t chg;
		sscanf(buf, "%lu\t%*s\t%u\t%lu\t%lf", &pos, &cnt, &chg, &qual);
		//cout << pos << "\t" << cnt << "\t" << qual << endl;
		kmernos->push_back(pos);
		kmers->push_back(new KMerCount(PositionKMer(pos), KMerStat(cnt, chg, qual)));
	}
	ifs.close();
}

void fillInSolidKmersFromFile( const string & fname, vector<KMerCount*> *kmers ) {
	ifstream ifs(fname);
	char buf[16000];
	while (!ifs.eof()) {
		ifs.getline(buf, 16000);
		ifs.getline(buf, 16000);
		hint_t pos; int cnt; double qual;
		sscanf(buf, ">%lu  cnt=%u  tql=%lf", &pos, &cnt, &qual);
		//cout << "kmerno =" << Globals::kmernos << endl;
		//cout << "kmerno size=" << Globals::kmernos->size() << endl;
		vector<hint_t>::const_iterator it_vec = lower_bound(Globals::kmernos->begin(), Globals::kmernos->end(), pos, PositionKMer::compareKMersDirect );
		//cout << "lower_bound=" << (it_vec - Globals::kmernos->begin()) << "  kmerssize=" << kmers->size() << "  nossize=" << Globals::kmernos->size() << endl;
		if ( *it_vec == pos ) {
			//cout << "solid: pos=" << pos << "\tcnt=" << cnt << "\tqual=" << qual << endl;
			kmers->at(it_vec - Globals::kmernos->begin())->second.makeGoodForIterative();
		}
	}
	ifs.close();
}


void fillInBadKmersFromFile( const string & fname, vector<KMerCount*> *kmers ) {
	ifstream ifs(fname);
	char buf[16000];
	while (!ifs.eof()) {
		ifs.getline(buf, 16000);
		ifs.getline(buf, 16000);
		hint_t pos; hint_t poscent;
		if ( sscanf(buf, ">%lu part of cluster %lu", &pos, &poscent) < 2 ) continue;
		vector<hint_t>::const_iterator it_vec = lower_bound(Globals::kmernos->begin(), Globals::kmernos->end(), pos, PositionKMer::compareKMersDirect );
		if ( *it_vec == pos ) {
			//cout << "pos=" << pos << "\tchangeto=" << poscent << endl;
			vector<hint_t>::const_iterator it_vec_tochange = lower_bound(Globals::kmernos->begin(), Globals::kmernos->end(), poscent, PositionKMer::compareKMersDirect);
			if ( *it_vec_tochange != poscent) cout << "   BAD GOOD KMER " << poscent << "\tlower_bound=" << it_vec_tochange - Globals::kmernos->begin() << endl;
			kmers->at(it_vec - Globals::kmernos->begin())->second.changeto = it_vec_tochange - Globals::kmernos->begin();
		}
	}
	ifs.close();
}


void dumpBlob( const string & fname ) {
	ofstream ofs(fname);
	ofs << Globals::blob_max_size << "\t" << Globals::blob_size << endl;
	ofs << Globals::blob << endl;
	ofs << Globals::blobquality << endl;
	ofs.close();
	cout << "deleting... " << endl;
	delete [] Globals::blob;
	delete [] Globals::blobquality;
	Globals::blob = NULL; Globals::blobquality = NULL;
	cout << "ok" << endl;
}

void loadBlob( const string & fname ) {
	cout << "deleting... " << endl;
	if (Globals::blob != NULL) delete [] Globals::blob;
	if (Globals::blobquality != NULL) delete [] Globals::blobquality;
	cout << "ok" << endl;
	ifstream ifs(fname);
	char buf[1024];
	ifs.getline(buf, 1024);
	sscanf(buf, "%lu\t%lu", &Globals::blob_max_size, &Globals::blob_size);
	Globals::blob = new char[Globals::blob_max_size];
	Globals::blobquality = new char[Globals::blob_max_size];
	ifs.getline(Globals::blob, Globals::blob_max_size);
	ifs.getline(Globals::blobquality, Globals::blob_max_size);
	ifs.close();
	cout << "New max_size=" << Globals::blob_max_size << "\tsize=" << Globals::blob_size << "\treal_size=" << strlen(Globals::blob) << endl;
}


void getGlobalConfigParameters( const string & config_file ) {
	TIMEDLN("Loading config from " << config_file.c_str());
	cfg::create_instance(config_file);
	Globals::working_dir = cfg::get().working_dir;
	Globals::qvoffset = cfg::get().quality_offset;
	Globals::error_rate = cfg::get().error_rate;
	Globals::blocksize_quadratic_threshold = cfg::get().blocksize_quadratic_threshold;
	Globals::good_cluster_threshold = cfg::get().good_cluster_threshold;
	Globals::blob_margin = cfg::get().blob_margin;
	Globals::trim_quality = cfg::get().trim_quality;
	Globals::trim_left_right = cfg::get().trim_left_right;
	Globals::use_iterative_reconstruction = cfg::get().use_iterative_reconstruction;
	Globals::iterative_reconstruction_threshold = cfg::get().iterative_reconstruction_threshold;
	Globals::max_reconstruction_iterations = cfg::get().max_reconstruction_iterations;
	Globals::reconstruction_in_full_iterations = cfg::get().reconstruction_in_full_iterations;
	Globals::read_kmers_after_clustering = cfg::get().read_kmers_after_clustering;
	Globals::write_kmers_after_clustering = cfg::get().write_kmers_after_clustering;
	Globals::kmers_after_clustering = cfg::get().kmers_after_clustering;
	Globals::write_each_iteration_kmers = cfg::get().write_each_iteration_kmers;
	Globals::regular_threshold_for_correction = cfg::get().regular_threshold_for_correction;
	Globals::discard_only_singletons = cfg::get().discard_only_singletons;
	Globals::special_nonsingleton_threshold = cfg::get().special_nonsingleton_threshold;
	Globals::use_true_likelihood = cfg::get().use_true_likelihood;
	Globals::num_of_tmp_files = cfg::get().num_of_tmp_files;
	Globals::conserve_memory = cfg::get().conserve_memory;
	Globals::paired_reads = cfg::get().paired_reads;
	Globals::skip_to_clustering = cfg::get().skip_to_clustering;
	Globals::skip_to_subvectors = cfg::get().skip_to_subvectors;
	Globals::skip_sorting_subvectors = cfg::get().skip_sorting_subvectors;
	Globals::skip_cluster_merging = cfg::get().skip_cluster_merging;
	Globals::hamming_class_buffer = cfg::get().hamming_class_buffer;
	Globals::unload_blob_before_merge = cfg::get().unload_blob_before_merge;
	Globals::debug_output_clustering = cfg::get().debug_output_clustering;
	Globals::debug_output_likelihood = cfg::get().debug_output_likelihood;
	Globals::skip_iterative = cfg::get().skip_iterative;
	Globals::likelihood_e_step = cfg::get().likelihood_e_step;
	Globals::subtract_simplex_volume = cfg::get().subtract_simplex_volume;
	Globals::change_n_to_random = cfg::get().change_n_to_random;
}

