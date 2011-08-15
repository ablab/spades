#ifndef POSITION_KMER_HPP_
#define POSITION_KMER_HPP_

#include <math.h>
#include <vector>
#include <queue>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include "read/read.hpp"
#include "kmer_stat.hpp"
#include "position_read.hpp"

#define K 55

/**
  * convert quality value to actual probability
  */
inline double qual2prob(uint8_t qual) {
	if (qual < 3) return 0.25;
	static std::vector<double> prob(255, -1);
	if (prob[qual] < -0.1) {
		prob[qual] = 1 - pow(10.0, - qual / 10.0);
	}
	return prob[qual];
}

typedef std::pair<std::string, uint32_t> StringCount;

class PositionKMer {
	hint_t start_;

  public:
	static std::vector<PositionRead> * pr;
	static std::vector<Read> * rv;
	static std::vector<bool> * rv_bad;
	static std::vector<Read> * rvLeft;
	static std::vector<Read> * rvRight;
	static std::vector<Read> * rvLeft_bad;
	static std::vector<Read> * rvRight_bad;
	static hint_t revNo;
	static hint_t lastLeftNo;

	static char* blob;
	static char* blobquality;
	static hint_t blob_max_size;
	static hint_t blob_size;

	static hint_t* blobkmers;

	static std::vector<uint32_t> * subKMerPositions;

	static void writeBlob( const char * fname );
	static void readBlob( const char * fname );
	static void writeBlobKMers( const char * fname );
	static void readBlobKMers( const char * fname );
	static void writeKMerCounts( const char * fname, const vector<KMerCount> & kmers );
	static void readKMerCounts( const char * fname, vector<KMerCount> * kmers );

	static double getKMerQuality( const hint_t & index, const int qvoffset );

	static bool compareSubKMersCheq( const hint_t & kmer1, const hint_t & kmer2, const std::vector<KMerCount> * km, const uint32_t tauplusone, const uint32_t start) {
		for ( hint_t i = start; i < K; i += tauplusone ) {
			if ( blob[ km->at(kmer1).first.start_ + i ] != blob [ km->at(kmer2).first.start_ + i ] ) {
				return ( blob[ km->at(kmer1).first.start_ + i ] < blob [ km->at(kmer2).first.start_ + i ] );
			}
		}
		return false;
	}

	static bool compareSubKMersGreaterCheq( const hint_t & kmer1, const hint_t & kmer2, const std::vector<KMerCount> * km, const uint32_t tauplusone, const uint32_t start) {
		for ( hint_t i = start; i < K; i += tauplusone ) {
			if ( blob[ km->at(kmer1).first.start_ + i ] != blob [ km->at(kmer2).first.start_ + i ] ) {
				return ( blob[ km->at(kmer1).first.start_ + i ] > blob [ km->at(kmer2).first.start_ + i ] );
			}
		}
		return false;
	}

	static bool equalSubKMersCheq( const hint_t & kmer1, const hint_t & kmer2, const std::vector<KMerCount> * km, const uint32_t tauplusone, const uint32_t start) {
		for ( hint_t i = start; i < K; i += tauplusone ) {
			if ( blob[ km->at(kmer1).first.start_ + i ] != blob [ km->at(kmer2).first.start_ + i ] ) {
				return false;
			}
		}
		return true;
	}


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


#endif

