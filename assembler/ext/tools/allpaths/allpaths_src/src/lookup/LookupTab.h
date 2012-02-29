/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file LookupTab.h
 * \author tsharpe
 * \date Jun 15, 2009
 *
 * \brief A Feudal lookup table.
 * This is the class you use to read a table of lookups.  To create one you use LookupTabBuilder.
 * The outer-vector index is a kmer, and the inner-vector contains Locations where the kmer is found.
 */
#ifndef LOOKUP_LOOKUPTAB_H_
#define LOOKUP_LOOKUPTAB_H_
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"

class Location
{
public:
    Location()
    : mContigID(~0U), mOffset(~0U)
    {}

    Location( unsigned int contigID, unsigned int offset )
    : mContigID(contigID), mOffset(offset)
    {}

    // compiler-supplied copying and destructor are OK

    unsigned int getContig() const
    { return mContigID; }

    unsigned int getOffset() const
    { return mOffset; }

private:
    unsigned int mContigID;
    unsigned int mOffset;
};

TRIVIALLY_SERIALIZABLE(Location);
typedef SerfVec<Location> LocationVec;
typedef MasterVec<LocationVec> VecLocationVec;

class LookupTab : public VecLocationVec
{
public:
    LookupTab() { init(); }
    LookupTab( int K ) : VecLocationVec(1 << 2*K), mK(K) {}

    /// File must exist.  This constructor reads in a lookup table from the file.
    LookupTab( String const& fileName ) : VecLocationVec(fileName) { init(); }

    /// K, the length of the kmer, for this table.
    unsigned int getK() const
    { return mK; }

    /// the usual definition: 2-bits per base with the more 5' bases occupying the most significant bits
    typedef unsigned int kmer_t;

    /// Find all locations for a specified kmer.
    LocationVec& getLocations( kmer_t kmer )
    { return operator[](kmer); }

    /// Find all locations (const version) for a specified kmer.
    LocationVec const& getLocations( kmer_t kmer ) const
    { return operator[](kmer); }

    /// Find the number of locations for a specified kmer.
    unsigned int getFreq( kmer_t kmer ) const
    { return operator[](kmer).size(); }

    /// Helper function to turn a base sequence into a kmer.
    /// The Itr class is an input iterator that derefs to a code for a base (i.e., a number from 0 to 3).
    template <class Itr> static unsigned int getKmer( Itr start, Itr const& end )
    { unsigned int result = 0;
      while ( start != end )
      { result = (result << 2) | *start; ++start; }
      return result; }

    /// Helper function to turn a base sequence into a kmer.
    /// The Itr class is a random-access iterator that derefs to a code for a base (i.e., a number from 0 to 3).
    template <class Itr> unsigned int getKmer( Itr start ) const
    { return getKmer(start,start+mK); }

private:
    void init()
    { unsigned int nnn = size();
      mK = 0;
      while ( nnn >>= 2 ) mK += 1; }

    unsigned int mK;
};


#endif /* LOOKUP_LOOKUPTAB_H_ */
