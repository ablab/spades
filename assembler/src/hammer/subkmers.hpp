//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef SUBKMERS_HPP_
#define SUBKMERS_HPP_

#include <vector>
#include <queue>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "read/read.hpp"
#include "kmer_stat.hpp"
#include "position_read.hpp"
#include "position_kmer.hpp"
#include "hammer_tools.hpp"


struct SubKMerPQElement; // forward declaration
typedef boost::function< bool (const hint_t & kmer1, const hint_t & kmer2) > SubKMerFunction;
typedef boost::function< bool (const SubKMerPQElement & kmer1, const SubKMerPQElement & kmer2) > SubKMerCompType;

// these are classes for the subkmer priority queue -- a result of parallel sort

struct SubKMerPQElement {
	hint_t ind;
	int n;
	SubKMerPQElement( hint_t index, int no) : ind(index), n(no) { }

	static bool compareSubKMerPQElements( const SubKMerPQElement & kmer1, const SubKMerPQElement & kmer2, const std::vector<KMerCount> * km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return PositionKMer::compareSubKMersGreater( kmer1.ind, kmer2.ind, km, tau, start_offset, end_offset );
	}

	static bool compareSubKMerPQElementsCheq( const SubKMerPQElement & kmer1, const SubKMerPQElement & kmer2, const std::vector<KMerCount> * km, const uint32_t tau, const uint32_t start) {
		return PositionKMer::compareSubKMersGreaterCheq( kmer1.ind, kmer2.ind, km, tau, start );
	}

	static bool compareSubKMerPQElementsDirect( const SubKMerPQElement & kmer1, const SubKMerPQElement & kmer2, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return PositionKMer::compareSubKMersGreaterDirect( kmer1.ind, kmer2.ind, tau, start_offset, end_offset );
	}

	static bool compareSubKMerPQElementsHInt( const SubKMerPQElement & kmer1, const SubKMerPQElement & kmer2, const std::vector<hint_t> * km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return PositionKMer::compareSubKMersGreaterHInt( kmer1.ind, kmer2.ind, km, tau, start_offset, end_offset );
	}

	static bool functionSubKMerPQElement(const SubKMerPQElement & kmer1, const SubKMerPQElement & kmer2, SubKMerFunction kmer_function) {
		return kmer_function( kmer1.ind, kmer2.ind );
	}
};

class SubKMerPQ {
  private:
	vector< size_t > boundaries;
	vector<hint_t> * v;
	int nthreads;

	SubKMerCompType sort_routine;
	std::priority_queue< SubKMerPQElement, vector<SubKMerPQElement>, SubKMerCompType  > pq;
	vector< vector<hint_t>::iterator > it;
	vector< vector<hint_t>::iterator > it_end;

	vector< string > fnames_;
	vector< boost::shared_ptr<FIStream> > ifs_;
	string buf_;

	hint_t getNextElementFromFile(size_t j);

  public:
	/**
	  * constructor
	  */
	SubKMerPQ( vector<hint_t> * vec, int nthr, SubKMerCompType sort_routine );

	/**
	  * constructor
	  */
	SubKMerPQ( const vector<string> & fnames, int nthr, SubKMerCompType sort_routine );

	/**
	  * sort one subvector array j (only one for easy parallelization)
	  */
	void doSort(int j, const SubKMerFunction & sub_sort);

	/**
	  * initialize priority queue
	  */
	void initPQ();

	/**
	  * get next priority queue element and pop the top
	  */
	hint_t nextPQ();

	/**
	  * peek at next priority queue element
	  */
	hint_t peekPQ() const { return pq.top().ind; }

	/**
	  * pop the top
	  */
	void popPQ();

	/**
	  * close all files and clean up
	  */
	void closePQ();

	/**
	  * is priority queue empty
	  */
	bool emptyPQ() const { return ( pq.size() == 0 ); }

	/// get boundaries
	const vector< size_t > & get_boundaries() const { return boundaries; }
};


/**
  * @class SubKMerSorter
  * Abstraction for parallel sorting of sub-kmers in kmers.
  */
class SubKMerSorter {
  public:

	enum SubKMerSorterType {
		SorterTypeStraight,
		SorterTypeChequered,
		SorterTypeChequeredDirect,
		SorterTypeFileBasedStraight
	};
	
	/**
	  * constructor for sorting a new kmers vector
	  */
	SubKMerSorter( size_t kmers_size, vector<KMerCount> * k, int nthreads, int tau, SubKMerSorterType type );

	/**
	  * constructor for sorting a new kmers vector
	  */
	SubKMerSorter( size_t kmers_size, vector<hint_t> * k, int nthreads, int tau, SubKMerSorterType type );

	/**
	  * constructor for sorting a specific block of kmers
	  * @param jj means the index of subkmer we want to *exclude*
	  * (this block results from the sorting of new kmers, so this block is automatically equal for all kmers in this block)
	  * @param parent_type type of the sorter with which this block was produced
	  */
	SubKMerSorter( vector< hint_t > * kmers, vector<KMerCount> * k, int nthreads, int tau, int jj,
		SubKMerSorterType type, SubKMerSorterType parent_type );

	/**
	  * constructor; similar to the one above but with hint_t kmerno vector
	  */
	SubKMerSorter( vector< hint_t > * kmers, vector<hint_t> * k, int nthreads, int tau, int jj,
		SubKMerSorterType type, SubKMerSorterType parent_type );

	/**
	  * destructor
	  */
	~SubKMerSorter();

	/**
	  * run parallel sort
	  */
	void runSort(std::string inputFile = "");

	/**
	  * get next block w.r.t. the sub_equal function from the i-th subkmer (in place)
	  * run this only after runSort();
	  * @return true if a block is returned; false is the queue is empty
	  */
	bool getNextBlock( int i, vector<hint_t> & block );

	/**
	 * is it file-based
	 */
	bool isFileBased() const { return type_ == SorterTypeFileBasedStraight; }

  private:
	SubKMerSorterType type_; // type of this sorter
	vector<SubKMerFunction> sub_less;	// subkmer comparison: less
	vector<SubKMerFunction> sub_greater;	// subkmer comparison: greater
	vector<SubKMerFunction> sub_equal;	// subkmer comparison: equal
	vector< vector<hint_t> > * v_;		// vector of kmer indices
	vector< SubKMerPQ > vskpq_;		// vector of subkmer priority queues for outputting the results
	int nthreads_;				// number of threads
	int tau_;				// tau
	size_t kmers_size_;			// size of kmer vector -- if we are sorting a newly made kmers vector
	vector< hint_t > * kmers_;		// the kmer vector itself -- if we are sorting a specific block
	vector< string > fnames_; // vector of file names for file-based versions
	vector< string > sorted_fnames_; // vector of file names for file-based versions

	void initVectors();
	void runMemoryBasedSort();
	void runFileBasedSort(std::string inputFile);
};

#endif

