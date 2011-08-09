/*
 * position_kmer.cpp
 *
 *  Created on: 01.08.2011
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
#include "position_kmer.hpp"


SubKMerPQ::SubKMerPQ( vector<hint_t> * vec, int nthr, subkmer_comp_type sort_routine ) : boundaries(nthr + 1), v(vec), nthreads(nthr), pq(sort_routine), it(nthr), it_end(nthr), cur_min(0, 0) {
	// find boundaries of the pieces
	size_t sub_size = (size_t)(v->size() / nthreads);
	for (int j=0; j<nthreads; ++j) {
		boundaries[j] = j * sub_size;
	}
	boundaries[nthreads] = v->size();
}

void SubKMerPQ::doSort(int j, const boost::function< bool (const hint_t & kmer1, const hint_t & kmer2)  > & sub_sort) {
	sort(v->begin() + boundaries[j], v->begin() + boundaries[j+1], sub_sort);
}

void SubKMerPQ::initPQ() {
	for (int j=0; j<nthreads; ++j) {
		it[j] = v->begin() + boundaries[j];
		it_end[j] = v->begin() + boundaries[j+1];
		pq.push( SubKMerPQElement(*(it[j]), j) );
	}
	cur_min = pq.top();	
}

hint_t SubKMerPQ::nextPQ() {
	SubKMerPQElement pqel = pq.top(); pq.pop();
	++it[pqel.n];
	if ( it[pqel.n] != it_end[pqel.n] ) pq.push( SubKMerPQElement(*(it[pqel.n]), pqel.n) );	
	return pqel.ind;
}

