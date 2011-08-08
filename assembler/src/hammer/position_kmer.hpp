#ifndef POSITION_KMER_HPP_
#define POSITION_KMER_HPP_

#include <vector>
#include <queue>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include "read/read.hpp"
#include "kmer_stat.hpp"
#include "position_read.hpp"

#define K 55
#define GOOD_SINGLETON_THRESHOLD 1 
#define CONSENSUS_BLOB_MARGIN 0.1


typedef std::pair<std::string, uint32_t> StringCount;

class PositionKMer {
	hint_t start_;

  public:
	static std::vector<PositionRead> * pr;
	static std::vector<Read> * rv;
	static std::vector<bool> * rv_bad;
	static hint_t revNo;

	static char* blob;
	static hint_t blob_max_size;
	static hint_t blob_size;

	static hint_t* blobkmers;

	static std::vector<uint32_t> * subKMerPositions;

	static bool compareSubKMers( const hint_t & kmer1, const hint_t & kmer2, const std::vector<KMerCount> * km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return ( strncmp( blob + km->at(kmer1).first.start_ + start_offset,
			  	  blob + km->at(kmer2).first.start_ + start_offset,
				  end_offset - start_offset ) < 0 );
	}

	static bool compareSubKMersGreater( const hint_t & kmer1, const hint_t & kmer2, const std::vector<KMerCount> * km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return ( strncmp( blob + km->at(kmer1).first.start_ + start_offset,
			  	  blob + km->at(kmer2).first.start_ + start_offset,
				  end_offset - start_offset ) > 0 );
	}

	static bool equalSubKMers( const hint_t & kmer1, const hint_t & kmer2, const std::vector<KMerCount> * km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return ( strncmp( blob + km->at(kmer1).first.start_ + start_offset,
			  	  blob + km->at(kmer2).first.start_ + start_offset,
				  end_offset - start_offset ) == 0 );
	}

  	static hint_t readNoFromBlobPosInternal( hint_t blobpos, hint_t start, hint_t end ) {
		if (start >= end - 1) return start;
		hint_t mid = start + (end - start) / 2;
		if ( blobpos < pr->at(mid).start() ) {
			return readNoFromBlobPosInternal( blobpos, start, mid );
		} else {
			return readNoFromBlobPosInternal( blobpos, mid, end );
		}
	}

	static hint_t readNoFromBlobPos( hint_t blobpos ) {
		return readNoFromBlobPosInternal ( blobpos, 0, pr->size() );
	}



	PositionKMer( hint_t readno, uint32_t startpos ) {
		start_ = pr->at(readno).start() + startpos;
	}

	PositionKMer( hint_t startpos ) {
		start_ = startpos;
	}

	char at(uint32_t pos) const {
		return blob[ start_ + pos ];
	}

	char operator [] (hint_t pos) const {
		return blob[ start_ + pos ];
	}

	bool operator < ( const PositionKMer & kmer ) const {
		return ( strncmp( blob + start_, blob + kmer.start_, K)  < 0 );
	}
	
	bool operator == ( const PositionKMer & kmer ) const {
		return ( strncmp( blob + start_, blob + kmer.start_, K) == 0 );
	}

	hint_t start() const { return start_; }

	string str() const {
		string res = "";
		for (uint32_t i = 0; i < K; ++i) {
			res += at(i);
		}
		return res;
	}

	string strSub(uint32_t tau, uint32_t offset) const {
		string res = "";
		for (uint32_t i = PositionKMer::subKMerPositions->at(offset); i < PositionKMer::subKMerPositions->at(offset+1); ++i) {
			res += at(i);
		}
		return res;
	}

};

inline bool KCgreater ( const KMerCount & l, const KMerCount & r ) {
	return l.first < r.first;
}

struct KMerNo {
	hint_t index;

	KMerNo( hint_t no ) : index(no) { } 

	bool equal(const KMerNo & kmerno) const {
		return ( strncmp( PositionKMer::blob + index, PositionKMer::blob + kmerno.index, K) == 0 );
	}

	bool test_equal(const KMerNo & kmerno) const {
		return ( index == kmerno.index );
	}

	string str() const {
		string res = "";
		for (uint32_t i = 0; i < K; ++i) {
			res += PositionKMer::blob[ index + i ];
		}
		return res;
	}

	static bool less(const KMerNo &l, const KMerNo &r) {
		return ( strncmp( PositionKMer::blob + l.index, PositionKMer::blob + r.index, K) < 0 );

	}

	static bool greater(const KMerNo &l, const KMerNo &r) {
		return ( strncmp( PositionKMer::blob + l.index, PositionKMer::blob + r.index, K) > 0 );

	}

	static bool test_less(const KMerNo &l, const KMerNo &r) {
		return ( l.index < r.index );

	}
	static bool test_greater(const KMerNo &l, const KMerNo &r) {
		return ( l.index > r.index );

	}

};


// these are classes for the subkmer priority queue -- a result of parallel sort

struct SubKMerPQElement {
	hint_t ind;
	int n;
	SubKMerPQElement( hint_t index, int no) : ind(index), n(no) { }

	static bool compareSubKMerPQElements( const SubKMerPQElement & kmer1, const SubKMerPQElement & kmer2, const std::vector<KMerCount> * km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return PositionKMer::compareSubKMersGreater( kmer1.ind, kmer2.ind, km, tau, start_offset, end_offset );
	}
	
};

typedef boost::function< bool (const SubKMerPQElement & kmer1, const SubKMerPQElement & kmer2) > subkmer_comp_type;

class SubKMerPQ {
  private:
	vector< size_t > boundaries;
	vector<hint_t> * v;
	int nthreads;

	subkmer_comp_type sort_routine;
	std::priority_queue< SubKMerPQElement, vector<SubKMerPQElement>, subkmer_comp_type  > pq;
	vector< vector<hint_t>::iterator > it;
	vector< vector<hint_t>::iterator > it_end;
	SubKMerPQElement cur_min;

  public:
	/**
	  * constructor
	  */
	SubKMerPQ( vector<hint_t> * vec, int nthr, subkmer_comp_type sort_routine );

	/**
	  * sort one subvector array j (only one for easy parallelization)
	  */
	void doSort(int j, const boost::function< bool (const hint_t & kmer1, const hint_t & kmer2)  > & sub_sort);

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
	hint_t peekPQ() { return cur_min.ind; }

	/**
	  * is priority queue empty
	  */
	bool emptyPQ() { return ( pq.size() == 0 ); }

	/// get boundaries
	const vector< size_t > & get_boundaries() { return boundaries; }
};


#endif

