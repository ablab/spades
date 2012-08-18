//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef POSITION_KMER_HPP_
#define POSITION_KMER_HPP_

#include <math.h>
#include <vector>
#include "sequence/seq.hpp"
#include "read/read.hpp"
#include "kmer_stat.hpp"
#include "position_read.hpp"
#include "config_struct_hammer.hpp"

class PositionKMer {
	hint_t start_;

  public:

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

	static bool equalKMersDirect( hint_t kmer1, hint_t kmer2) {
		return ( strncmp( Globals::blob + kmer1, Globals::blob + kmer2, K ) == 0 );
	}

	static bool compareKMersDirect( hint_t kmer1, hint_t kmer2) {
		return ( strncmp( Globals::blob + kmer1, Globals::blob + kmer2, K ) < 0 );
	}

  PositionKMer( hint_t readno, uint32_t startpos ) : start_(Globals::pr->at(readno).start() + startpos) { }
	explicit PositionKMer(hint_t startpos): start_(startpos) { }
	PositionKMer() : start_(-1ULL) { }

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

  template<class Writer>
  friend Writer& binary_write(Writer &os, const PositionKMer &pos);
  template<class Reader>
  friend void binary_read(Reader &is, PositionKMer &pos);
};

template<class Writer>
inline Writer& binary_write(Writer &os, const PositionKMer &pos) {
  os.write((char*)&pos.start_, sizeof(pos.start_));

  return os;
}

template<class Reader>
inline void binary_read(Reader &is, PositionKMer &pos) {
  is.read((char*)&pos.start_, sizeof(pos.start_));
}

#endif

