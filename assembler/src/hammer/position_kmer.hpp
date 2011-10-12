#ifndef POSITION_KMER_HPP_
#define POSITION_KMER_HPP_

#include <math.h>
#include <vector>
#include <queue>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/unordered_map.hpp>
#include "sequence/seq.hpp"
#include "read/read.hpp"
#include "kmer_stat.hpp"
#include "position_read.hpp"

class PositionKMer {
	hint_t start_;

  public:

	static bool compareSubKMersCheq( const hint_t & kmer1, const hint_t & kmer2, const std::vector<KMerCount*> * km, const uint32_t tauplusone, const uint32_t start) {
		for ( uint32_t i = start; i < K; i += tauplusone ) {
			if ( Globals::blob[ km->at(kmer1)->first.start_ + i ] != Globals::blob [ km->at(kmer2)->first.start_ + i ] ) {
				return ( Globals::blob[ km->at(kmer1)->first.start_ + i ] < Globals::blob [ km->at(kmer2)->first.start_ + i ] );
			}
		}
		return false;
	}

	static bool compareSubKMersGreaterCheq( const hint_t & kmer1, const hint_t & kmer2, const std::vector<KMerCount*> * km, const uint32_t tauplusone, const uint32_t start) {
		for ( uint32_t i = start; i < K; i += tauplusone ) {
			if ( Globals::blob[ km->at(kmer1)->first.start_ + i ] != Globals::blob [ km->at(kmer2)->first.start_ + i ] ) {
				return ( Globals::blob[ km->at(kmer1)->first.start_ + i ] > Globals::blob [ km->at(kmer2)->first.start_ + i ] );
			}
		}
		return false;
	}

	static bool equalSubKMersCheq( const hint_t & kmer1, const hint_t & kmer2, const std::vector<KMerCount*> * km, const uint32_t tauplusone, const uint32_t start) {
		for ( uint32_t i = start; i < K; i += tauplusone ) {
			if ( Globals::blob[ km->at(kmer1)->first.start_ + i ] != Globals::blob [ km->at(kmer2)->first.start_ + i ] ) {
				return false;
			}
		}
		return true;
	}


	static bool compareSubKMers( const hint_t & kmer1, const hint_t & kmer2, const std::vector<KMerCount*> * km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return ( strncmp( Globals::blob + km->at(kmer1)->first.start_ + start_offset,
			  	  Globals::blob + km->at(kmer2)->first.start_ + start_offset,
				  end_offset - start_offset ) < 0 );
	}

	static bool compareSubKMersGreater( const hint_t & kmer1, const hint_t & kmer2, const std::vector<KMerCount*> * km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return ( strncmp( Globals::blob + km->at(kmer1)->first.start_ + start_offset,
			  	  Globals::blob + km->at(kmer2)->first.start_ + start_offset,
				  end_offset - start_offset ) > 0 );
	}

	static bool equalSubKMers( const hint_t & kmer1, const hint_t & kmer2, const std::vector<KMerCount*> * km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return ( strncmp( Globals::blob + km->at(kmer1)->first.start_ + start_offset,
			  	  Globals::blob + km->at(kmer2)->first.start_ + start_offset,
				  end_offset - start_offset ) == 0 );
	}

  	static hint_t readNoFromBlobPosInternal( hint_t blobpos, hint_t start, hint_t end ) {
		if (start >= end - 1) return start;
		hint_t mid = start + (end - start) / 2;
		if ( blobpos < Globals::pr->at(mid).start() ) {
			return readNoFromBlobPosInternal( blobpos, start, mid );
		} else {
			return readNoFromBlobPosInternal( blobpos, mid, end );
		}
	}

	static hint_t readNoFromBlobPos( hint_t blobpos ) {
		return readNoFromBlobPosInternal ( blobpos, 0, Globals::pr->size() );
	}



	PositionKMer( hint_t readno, uint32_t startpos ) : start_(Globals::pr->at(readno).start() + startpos) { }
	PositionKMer( hint_t startpos ) : start_(startpos) { }
	PositionKMer() : start_(-1) { }


	char at(uint32_t pos) const {
		return Globals::blob[ start_ + pos ];
	}

	char operator [] (hint_t pos) const {
		return Globals::blob[ start_ + pos ];
	}

	int totalQual (hint_t pos) const {
		return Globals::totalquality[ start_ + pos ];
	}

	bool operator < ( const PositionKMer & kmer ) const {
		return ( strncmp( Globals::blob + start_, Globals::blob + kmer.start_, K)  < 0 );
	}
	
	bool operator == ( const PositionKMer & kmer ) const {
		return ( strncmp( Globals::blob + start_, Globals::blob + kmer.start_, K) == 0 );
	}

	hint_t start() const { return start_; }

	std::string str() const {
		std::string res = "";
		for (uint32_t i = 0; i < K; ++i) {
			res += at(i);
		}
		return res;
	}

	std::string strQual() const {
		std::string res = "";
		for (uint32_t i = 0; i < K; ++i) {
			res += Globals::blobquality[ start_ + i ];
		}
		return res;
	}

	std::string strSub(uint32_t tau, uint32_t offset) const {
		std::string res = "";
		for (uint32_t i = Globals::subKMerPositions->at(offset); i < Globals::subKMerPositions->at(offset+1); ++i) {
			res += at(i);
		}
		return res;
	}

};

inline bool KCgreater ( const KMerCount & l, const KMerCount & r ) {
	return l.first < r.first;
}


#endif

