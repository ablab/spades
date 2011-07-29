#ifndef POSITION_KMER_HPP_
#define POSITION_KMER_HPP_

#include <vector>
#include "read/read.hpp"
#include "kmer_stat.hpp"
#include "position_read.hpp"

#define K 55
#define GOOD_SINGLETON_THRESHOLD 1 
#define CONSENSUS_BLOB_MARGIN 0.1


typedef std::pair<std::string, uint32_t> StringCount;

class PositionKMer {
	uint64_t start_;

  public:
	static std::vector<PositionRead> * pr;
	static std::vector<Read> * rv;
	static uint64_t revNo;

	static char* blob;
	static uint64_t blob_max_size;
	static uint64_t blob_size;

	static std::vector<uint32_t> * subKMerPositions;

	static bool compareSubKMers( const uint64_t kmer1, const uint64_t kmer2, const std::vector<KMerCount> * km, const uint32_t tau, const uint32_t offset) {
		return ( strncmp( blob + km->at(kmer1).first.start_ + subKMerPositions->at(offset),
			  	  blob + km->at(kmer2).first.start_ + subKMerPositions->at(offset),
				  subKMerPositions->at(offset + 1) - subKMerPositions->at(offset)) < 0 );
	}

	static bool equalSubKMers( const uint64_t kmer1, const uint64_t kmer2, const std::vector<KMerCount> * km, const uint32_t tau, const uint32_t offset) {
		for (uint32_t i = offset; i < K; i += tau+1) {
			if (km->at(kmer1).first[i] != km->at(kmer2).first[i]) {
				return false;
			}
			i += tau+1;
		}
		return true;
	}

  	static uint64_t readNoFromBlobPosInternal( uint64_t blobpos, uint64_t start, uint64_t end ) {
		if (start >= end - 1) return start;
		uint64_t mid = start + (end - start) / 2;
		if ( blobpos < pr->at(mid).start() ) {
			return readNoFromBlobPosInternal( blobpos, start, mid );
		} else {
			return readNoFromBlobPosInternal( blobpos, mid, end );
		}
	}

	static uint64_t readNoFromBlobPos( uint64_t blobpos ) {
		return readNoFromBlobPosInternal ( blobpos, 0, pr->size() );
	}



	PositionKMer( uint64_t readno, uint32_t startpos ) {
		start_ = pr->at(readno).start() + startpos;
	}

	PositionKMer( uint64_t startpos ) {
		start_ = startpos;
	}

	char at(uint32_t pos) const {
		return blob[ start_ + pos ];
	}

	char operator [] (uint32_t pos) const {
		return blob[ start_ + pos ];
	}

	bool operator < ( const PositionKMer & kmer ) const {
		return ( strncmp( blob + start_, blob + kmer.start_, K)  < 0 );
	}
	
	bool operator == ( const PositionKMer & kmer ) const {
		return ( strncmp( blob + start_, blob + kmer.start_, K) == 0 );
	}

	string str() const {
		string res = "";
		for (uint32_t i = 0; i < K; ++i) {
			res += at(i);
		}
		return res;
	}

	string strSub(uint32_t tau, uint32_t offset) const {
		string res = "";
		for (uint32_t i = offset; i < K; i+=tau+1) {
			res += at(i);
		}
		return res;
	}

};

inline bool KCgreater ( const KMerCount & l, const KMerCount & r ) {
	return l.first < r.first;
}

struct KMerNo {
	uint64_t index;

	KMerNo( uint64_t no ) : index(no) { } 

	bool equal(const KMerNo & kmerno) {
		return ( strncmp( PositionKMer::blob + index, PositionKMer::blob + kmerno.index, K) == 0 );
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
};


#endif

