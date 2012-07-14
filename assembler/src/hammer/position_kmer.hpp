//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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
#include "config_struct_hammer.hpp"

class PositionKMer {
	hint_t start_;

  public:

	static bool compareSubKMersCheq(hint_t kmer1, hint_t kmer2, const std::vector<KMerCount> *km, const uint32_t tauplusone, const uint32_t start) {
		for ( uint32_t i = start; i < K; i += tauplusone ) {
			if ( Globals::blob[ (*km)[kmer1].first.start_ + i ] != Globals::blob [ (*km)[kmer2].first.start_ + i ] ) {
				return ( Globals::blob[ (*km)[kmer1].first.start_ + i ] < Globals::blob [ (*km)[kmer2].first.start_ + i ] );
			}
		}
		return false;
	}

	static bool compareSubKMersGreaterCheq( hint_t kmer1, hint_t kmer2, const std::vector<KMerCount> *km, const uint32_t tauplusone, const uint32_t start) {
		for ( uint32_t i = start; i < K; i += tauplusone ) {
			if ( Globals::blob[ (*km)[kmer1].first.start_ + i ] != Globals::blob [ (*km)[kmer2].first.start_ + i ] ) {
				return ( Globals::blob[ (*km)[kmer1].first.start_ + i ] > Globals::blob [ (*km)[kmer2].first.start_ + i ] );
			}
		}
		return false;
	}

	static bool equalSubKMersCheq( hint_t kmer1, hint_t kmer2, const std::vector<KMerCount> *km, const uint32_t tauplusone, const uint32_t start) {
		for ( uint32_t i = start; i < K; i += tauplusone ) {
			if ( Globals::blob[ (*km)[kmer1].first.start_ + i ] != Globals::blob [ (*km)[kmer2].first.start_ + i ] ) {
				return false;
			}
		}
		return true;
	}

	static bool compareSubKMersCheqDirect( hint_t kmer1, hint_t kmer2, const uint32_t tauplusone, const uint32_t start) {
		for ( uint32_t i = start; i < K; i += tauplusone ) {
			if ( Globals::blob[ kmer1 + i ] != Globals::blob [ kmer2 + i ] ) {
				return ( Globals::blob[ kmer1 + i ] < Globals::blob [ kmer2 + i ] );
			}
		}
		return false;
	}

	static bool compareSubKMersGreaterCheqDirect( hint_t kmer1, hint_t kmer2, const uint32_t tauplusone, const uint32_t start) {
		for ( uint32_t i = start; i < K; i += tauplusone ) {
			if ( Globals::blob[ kmer1 + i ] != Globals::blob [ kmer2 + i ] ) {
				return ( Globals::blob[ kmer1 + i ] > Globals::blob [ kmer2 + i ] );
			}
		}
		return false;
	}

	static bool equalSubKMersCheqDirect( hint_t kmer1, hint_t kmer2, const uint32_t tauplusone, const uint32_t start) {
		for ( uint32_t i = start; i < K; i += tauplusone ) {
			if ( Globals::blob[ kmer1 + i ] != Globals::blob [ kmer2 + i ] ) {
				return false;
			}
		}
		return true;
	}

	static bool compareSubKMersCheqHInt( hint_t kmer1, hint_t kmer2, const std::vector<hint_t> *km, const uint32_t tauplusone, const uint32_t start) {
		for ( uint32_t i = start; i < K; i += tauplusone ) {
			if ( Globals::blob[ (*km)[kmer1] + i ] != Globals::blob [ (*km)[kmer2] + i ] ) {
				return ( Globals::blob[ (*km)[kmer1] + i ] < Globals::blob [ (*km)[kmer2] + i ] );
			}
		}
		return false;
	}

	static bool compareSubKMersGreaterCheqHInt( hint_t kmer1, hint_t kmer2, const std::vector<hint_t> *km, const uint32_t tauplusone, const uint32_t start) {
		for ( uint32_t i = start; i < K; i += tauplusone ) {
			if ( Globals::blob[ (*km)[kmer1] + i ] != Globals::blob [ (*km)[kmer2] + i ] ) {
				return ( Globals::blob[ (*km)[kmer1] + i ] > Globals::blob [ (*km)[kmer2] + i ] );
			}
		}
		return false;
	}

	static bool equalSubKMersCheqHInt( hint_t kmer1, hint_t kmer2, const std::vector<hint_t> *km, const uint32_t tauplusone, const uint32_t start) {
		for ( uint32_t i = start; i < K; i += tauplusone ) {
			if ( Globals::blob[(*km)[kmer1] + i ] != Globals::blob [(*km)[kmer2] + i ] ) {
				return false;
			}
		}
		return true;
	}

	static bool compareSubKMers( hint_t kmer1, hint_t kmer2, const std::vector<KMerCount> *km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return ( strncmp( Globals::blob + (*km)[kmer1].first.start_ + start_offset,
                      Globals::blob + (*km)[kmer2].first.start_ + start_offset,
                      end_offset - start_offset ) < 0 );
	}

	static bool compareSubKMersGreater( hint_t kmer1, hint_t kmer2, const std::vector<KMerCount> *km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return ( strncmp( Globals::blob + (*km)[kmer1].first.start_ + start_offset,
                      Globals::blob + (*km)[kmer2].first.start_ + start_offset,
                      end_offset - start_offset ) > 0 );
	}

	static bool compareSubKMersGreaterSimple( const pair<hint_t, pair< double, size_t > > & kmer1, const pair<hint_t, pair< double, size_t > > & kmer2) {
		return ( strncmp( Globals::blob + kmer1.first, Globals::blob + kmer2.first, K ) > 0 );
	}

	static bool compareSubKMersLessSimple( const pair<hint_t, pair< double, size_t > > & kmer1, const pair<hint_t, pair< double, size_t > > & kmer2) {
		return ( strncmp( Globals::blob + kmer1.first, Globals::blob + kmer2.first, K ) < 0 );
	}

	static bool compareSubKMersGFirst( const pair<hint_t, pair< double, size_t > > & kmer1, const pair<hint_t, pair< double, size_t > > & kmer2) {
		for ( uint32_t i = 0; i < K; ++i ) {
			if ( Globals::blob[ kmer1.first + i ] != Globals::blob [ kmer2.first + i ] ) {
				switch ( Globals::blob[ kmer1.first + i ] ) {
				case 'G': return true;
				case 'A': return ( Globals::blob [ kmer2.first + i ] != 'G' );
				case 'T': return ( Globals::blob [ kmer2.first + i ] == 'C' );
				case 'C': return false;
				default: return false;
				}
			}
		}
		return false;
	}

	static bool compareSubKMersCFirst( const pair<hint_t, pair< double, size_t > > & kmer1, const pair<hint_t, pair< double, size_t > > & kmer2) {
		for ( uint32_t i = 0; i < K; ++i ) {
			if ( Globals::blob[ kmer1.first + i ] != Globals::blob [ kmer2.first + i ] ) {
				switch ( Globals::blob[ kmer1.first + i ] ) {
				case 'C': return true;
				case 'T': return ( Globals::blob [ kmer2.first + i ] != 'C' );
				case 'A': return ( Globals::blob [ kmer2.first + i ] == 'G' );
				case 'G': return false;
				default: return false;
				}
			}
		}
		return false;
	}

	static bool equalSubKMers( hint_t kmer1, hint_t kmer2, const std::vector<KMerCount> *km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return ( strncmp( Globals::blob + (*km)[kmer1].first.start_ + start_offset,
			  	  Globals::blob + (*km)[kmer2].first.start_ + start_offset,
				  end_offset - start_offset ) == 0 );
	}

	static bool compareSubKMersHInt( hint_t kmer1, hint_t kmer2, const std::vector<hint_t> *km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return ( strncmp( Globals::blob + (*km)[kmer1] + start_offset, Globals::blob + (*km)[kmer2] + start_offset,
                      end_offset - start_offset ) < 0 );
	}

	static bool compareSubKMersGreaterHInt( hint_t kmer1, hint_t kmer2, const std::vector<hint_t> *km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return ( strncmp( Globals::blob + (*km)[kmer1] + start_offset, Globals::blob + (*km)[kmer2] + start_offset,
                      end_offset - start_offset ) > 0 );
	}

	static bool equalSubKMersHInt( hint_t kmer1, hint_t kmer2, const std::vector<hint_t> *km, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return ( strncmp( Globals::blob + (*km)[kmer1] + start_offset, Globals::blob + (*km)[kmer2] + start_offset,
                      end_offset - start_offset ) == 0 );
	}

	static bool compareSubKMersDirect( hint_t kmer1, hint_t kmer2, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return ( strncmp( Globals::blob + kmer1 + start_offset, Globals::blob + kmer2 + start_offset, end_offset - start_offset ) < 0 );
	}

	static bool equalKMersDirect( hint_t kmer1, hint_t kmer2) {
		return ( strncmp( Globals::blob + kmer1, Globals::blob + kmer2, K ) == 0 );
	}

	static bool compareKMersDirect( hint_t kmer1, hint_t kmer2) {
		return ( strncmp( Globals::blob + kmer1, Globals::blob + kmer2, K ) < 0 );
	}

	static bool compareSubKMersGreaterDirect( hint_t kmer1, hint_t kmer2, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		return ( strncmp( Globals::blob + kmer1 + start_offset, Globals::blob + kmer2 + start_offset, end_offset - start_offset ) > 0 );
	}

	static bool equalSubKMersDirect( hint_t kmer1, hint_t kmer2, const uint32_t tau, const uint32_t start_offset, const uint32_t end_offset) {
		//cout << "      equalSubKMersDirect: kmer1=" << kmer1 << " kmer2=" << kmer2 << " start_offset=" << start_offset << " end_offset=" << end_offset << " max=" << strlen(Globals::blob) << endl;
		return ( strncmp( Globals::blob + kmer1 + start_offset, Globals::blob + kmer2 + start_offset, end_offset - start_offset ) == 0 );
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

	char operator [] (hint_t pos) const {
		return Globals::blob[ start_ + pos ];
	}

	char at(hint_t pos ) const
	{
		if (pos >= Globals::blob_max_size)
			throw std::out_of_range((boost::format("PositionKMer, max: %d, pos: %d") % Globals::blob_max_size % pos).str());

		return Globals::blob[ start_ + pos ];
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
    std::string res(Globals::blobquality + start_, K);
    char qv_offset = cfg::get().input_qvoffset;

    for (uint32_t i = 0; i < K; ++i)
      res[i] += qv_offset;

    return res;
  }

	std::string strSub(uint32_t tau, uint32_t offset) const {
		std::string res = "";
		for (uint32_t i = Globals::subKMerPositions->at(offset); i < Globals::subKMerPositions->at(offset+1); ++i) {
			res += at(i);
		}
		return res;
	}

  friend std::ostream& binary_write(std::ostream &os, const PositionKMer &pos);
  friend void binary_read(std::istream &is, PositionKMer &pos);
};

inline bool KCgreater ( const KMerCount & l, const KMerCount & r ) {
	return l.first < r.first;
}

inline std::ostream& binary_write(std::ostream &os, const PositionKMer &pos) {
  os.write((char*)&pos.start_, sizeof(pos.start_));

  return os;
}

inline void binary_read(std::istream &is, PositionKMer &pos) {
  is.read((char*)&pos.start_, sizeof(pos.start_));
}

inline std::ostream& binary_write(std::ostream &os, const KMerCount &k) {
  return binary_write(binary_write(os, k.first), k.second);
}

inline void binary_read(std::istream &is, KMerCount &k) {
  binary_read(is, k.first);
  binary_read(is, k.second);
}

inline char getQual(const KMerCount & kmc, size_t i) {
  if (Globals::use_common_quality)
    return Globals::common_quality * kmc.second.count;
  if (kmc.second.count == 1)
    return Globals::blobquality[kmc.first.start() + i];
  else
    return kmc.second.qual[i];
}

inline double getProb(const KMerCount &kmc, size_t i, bool log) {
  uint8_t qual = getQual(kmc, i);

  return (log ? Globals::quality_lprobs[qual] : Globals::quality_probs[qual]);
}

#endif

