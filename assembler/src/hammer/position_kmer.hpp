#ifndef POSITION_KMER_HPP_
#define POSITION_KMER_HPP_

#include <vector>
#include "kmer_stat.hpp"

#define K 55
#define GOOD_SINGLETON_THRESHOLD 1 
#define CONSENSUS_BLOB_MARGIN 0.1


class PositionKMer;

typedef map<PositionKMer, KMerStat> KMerStatMap;
typedef pair<PositionKMer, KMerStat> KMerCount;
typedef vector<KMerCount> KMerStatVector;

typedef pair<std::string, uint32_t> StringCount;


class PositionKMer {
	// uint64_t read_;
	uint32_t start_;

  public:
	static std::vector<ReadStat> * rv;
	static uint64_t revNo;
	static char* blob;

	static bool compareSubKMers( const uint64_t kmer1, const uint64_t kmer2, const std::vector<KMerCount> * km, const uint32_t tau, const uint32_t offset) {
		for (uint32_t i = offset; i < K; i += tau+1) {
			if (km->at(kmer1).first[i] != km->at(kmer2).first[i]) {
				return (  km->at(kmer1).first[i] < km->at(kmer2).first[i] );
			}
		}
		return false;
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

	PositionKMer( uint64_t readno, uint32_t startpos ) {
		start_ = rv->at(readno).blobpos + startpos;
	}

	char at(uint32_t pos) const {
		return blob[ start_ + pos ];
//		return rv->at(read_).read.getSequenceString()[ start_ + pos ];
	}

	char operator [] (uint32_t pos) const {
		return blob[ start_ + pos ];
//		return at(pos);
	}

	// virtual uint32_t size() const { return K; }

	bool operator < ( const PositionKMer & kmer ) const {
		for (uint32_t i = 0; i < K; i++) {
			if ( at(i) != kmer.at(i) ) {
				return ( at(i) < kmer.at(i) );
			}
		}
		return false;
	}
	
	bool operator == ( const PositionKMer & kmer ) const {
		for (uint32_t i = 0; i < K; ++i) {
			if ( at(i) != kmer.at(i) ) {
				return false;
			}
		}
		return true;
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

#endif

