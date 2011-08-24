#include <omp.h>
#include "subkmers.hpp"

void SubKMerSorter::runSort() {

	int subkmer_nthreads = max ( (tau_ + 1) * ( (int)(nthreads_ / (tau_ + 1)) ), tau_+1 );
	int effective_subkmer_threads = min(subkmer_nthreads, nthreads_);

	vector< SubKMerPQ > * vskpq = &vskpq_;

	// we divide each of (tau+1) subkmer vectors into nthreads/(tau+1) subvectors
	// as a result, we have subkmer_nthreads threads for sorting
	#pragma omp parallel for shared(vskpq) num_threads( effective_subkmer_threads )
	for (int j=0; j < subkmer_nthreads; ++j) {
		// for each j, we sort subvector (j/(tau_+1)) of the vector of subkmers at offset (j%(tau+1))
		(*vskpq)[ (j % (tau_+1)) ].doSort( j / (tau_+1), sub_less[(j % (tau_+1))] );
	}

	for (int j=0; j < tau_+1; ++j) {
		vskpq_[j].initPQ();
	}
}

bool SubKMerSorter::getNextBlock( int i, vector<hint_t> & block ) {
	block.clear();
	if ( vskpq_[i].emptyPQ() ) return false;
	hint_t last = vskpq_[i].peekPQ();
	while (!vskpq_[i].emptyPQ()) {
		hint_t cur = vskpq_[i].peekPQ();
		if ( sub_equal[i](last, cur) ) { //add to current reads
			block.push_back(cur);
			vskpq_[i].popPQ();
		} else {
			return true;
		}
	}
	return (block.size() > 0);
}

SubKMerSorter::SubKMerSorter( size_t kmers_size, vector<KMerCount*> * k, int nthreads, int tau, SubKMerSorterType type ) :
	nthreads_(nthreads), tau_(tau), kmers_size_(kmers_size), kmers_(NULL) {
	// we set the sorting functions depending on the type
	// here the sorting functions are regular sorting functions predefined in PositionKMer
	switch( type ) {
		case SorterTypeStraight:
			for (int j=0; j < tau+1; ++j) {
				sub_less.push_back(    boost::bind(PositionKMer::compareSubKMers,        _1, _2, k, tau, 
					Globals::subKMerPositions->at(j), Globals::subKMerPositions->at(j+1) ) );
				sub_greater.push_back( boost::bind(PositionKMer::compareSubKMersGreater, _1, _2, k, tau, 
					Globals::subKMerPositions->at(j), Globals::subKMerPositions->at(j+1) ) );
				sub_equal.push_back(   boost::bind(PositionKMer::equalSubKMers,          _1, _2, k, tau, 
					Globals::subKMerPositions->at(j), Globals::subKMerPositions->at(j+1) ) );
			}
			break;
		case SorterTypeChequered:
			for (int j=0; j < tau+1; ++j) {
				sub_less.push_back(    boost::bind(PositionKMer::compareSubKMersCheq,        _1, _2, k, tau, j) );
				sub_greater.push_back( boost::bind(PositionKMer::compareSubKMersGreaterCheq, _1, _2, k, tau, j) );
				sub_equal.push_back(   boost::bind(PositionKMer::equalSubKMersCheq,          _1, _2, k, tau, j) );
			}
			break;
	}

	// initialize the vectors
	initVectors();
}

SubKMerSorter::SubKMerSorter( vector< hint_t > * kmers, vector<KMerCount*> * k, int nthreads, int tau, int jj,
	SubKMerSorterType type, SubKMerSorterType parent_type ) : nthreads_(nthreads), tau_(tau), kmers_size_(kmers->size()), kmers_(kmers) {

	//cout << "    constructor nthreads=" << nthreads << " tau=" << tau << " kmerssize=" << kmers_size_ << " jj=" << jj << endl;
	// we set the sorting functions depending on the type
	// for sorting a specific block, we use sorting functions with specific exemptions depending on jj
	if( type == SorterTypeStraight) {
		assert(parent_type == SorterTypeStraight);

		vector< pair<uint32_t, uint32_t> > my_positions(tau+1);
		uint32_t left_size = Globals::subKMerPositions->at(jj);
		uint32_t right_size = K - Globals::subKMerPositions->at(jj+1);
		uint32_t total_size = left_size + right_size;
		uint32_t left_end = ( (tau + 1) * left_size ) / ( total_size );
		uint32_t increment = total_size / (tau+1);

		for (uint32_t i=0; i < left_end; ++i) my_positions[i] = make_pair( i * increment, (i+1) * increment );
		if (left_end > 0) my_positions[left_end-1].second = left_size;
		for (uint32_t i=left_end; i < (uint32_t)tau+1; ++i) my_positions[i] = make_pair(
			Globals::subKMerPositions->at(jj+1) + (i  -left_end) * increment,
			Globals::subKMerPositions->at(jj+1) + (i+1-left_end) * increment );
		if (jj < tau) my_positions[tau].second = K;

		for (int j=0; j < tau+1; ++j) {
			sub_less.push_back(    boost::bind(PositionKMer::compareSubKMers,        _1, _2, k, tau,
				my_positions[j].first, my_positions[j].second ) );
			sub_greater.push_back( boost::bind(PositionKMer::compareSubKMersGreater, _1, _2, k, tau,
				my_positions[j].first, my_positions[j].second ) );
			sub_equal.push_back(   boost::bind(PositionKMer::equalSubKMers,          _1, _2, k, tau,
				my_positions[j].first, my_positions[j].second ) );
		}
	} else if ( type == SorterTypeChequered ) {
		assert(parent_type == SorterTypeStraight); // yet to implement a chequered second level
		for (int j=0; j < tau+1; ++j) {
			sub_less.push_back(    boost::bind(PositionKMer::compareSubKMersCheq,        _1, _2, k, tau+1, j) );
			sub_greater.push_back( boost::bind(PositionKMer::compareSubKMersGreaterCheq, _1, _2, k, tau+1, j) );
			sub_equal.push_back(   boost::bind(PositionKMer::equalSubKMersCheq,          _1, _2, k, tau+1, j) );
		}
	}

	// initialize the vectors
	initVectors();
}

void SubKMerSorter::initVectors() {
	v_ = new vector< vector<hint_t> >(tau_+1);
	int nthreads_per_subkmer = max( (int)(nthreads_ / (tau_ + 1)), 1);

	for (int j=0; j<tau_+1; ++j) {
		(*v_)[j].resize( kmers_size_ );
		for (size_t m = 0; m < kmers_size_; ++m) (*v_)[j][m] = m;
		SubKMerCompType sort_greater = boost::bind(SubKMerPQElement::functionSubKMerPQElement, _1, _2, sub_greater[j]);
		SubKMerPQ skpq( &((*v_)[j]), nthreads_per_subkmer, sort_greater );
		vskpq_.push_back(skpq);
	}
}

SubKMerPQ::SubKMerPQ( vector<hint_t> * vec, int nthr, SubKMerCompType sort_routine ) : boundaries(nthr + 1), v(vec), nthreads(nthr), pq(sort_routine), it(nthr), it_end(nthr) {
	// find boundaries of the pieces
	size_t sub_size = (size_t)(v->size() / nthreads);
	for (int j=0; j<nthreads; ++j) {
		boundaries[j] = j * sub_size;
	}
	boundaries[nthreads] = v->size();
}

void SubKMerPQ::doSort(int j, const SubKMerFunction & sub_sort) {
	sort(v->begin() + boundaries[j], v->begin() + boundaries[j+1], sub_sort);
}

void SubKMerPQ::initPQ() {
	for (int j=0; j<nthreads; ++j) {
		it[j] = v->begin() + boundaries[j];
		it_end[j] = v->begin() + boundaries[j+1];
		pq.push( SubKMerPQElement(*(it[j]), j) );
	}
}

hint_t SubKMerPQ::nextPQ() {
	hint_t res = peekPQ(); popPQ();
	return res;
}

void SubKMerPQ::popPQ() {
	SubKMerPQElement pqel = pq.top(); pq.pop();
	++it[pqel.n];
	if ( it[pqel.n] != it_end[pqel.n] ) pq.push( SubKMerPQElement(*(it[pqel.n]), pqel.n) );	
}

