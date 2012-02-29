///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This file defines "fastavector", which stores bases in fasta format.

#ifndef FASTA_VECTOR
#define FASTA_VECTOR

#include "Basevector.h"
#include "Bitvector.h"
#include "Charvector.h"
#include "String.h"
#include "Superb.h"
#include "Vec.h"
#include "dna/Bases.h"
#include <cstddef>

class fastavector : public charvector
{
    typedef charvector::alloc_type Alloc;

public:
    fastavector() {}

    explicit fastavector( Alloc const& alloc ) : charvector(alloc) {}

    explicit fastavector( size_type n, char val = GeneralizedBase::X.asChar(),
                          size_type capacity = 0, Alloc const& alloc = Alloc() )
    : charvector(n,val,capacity,alloc) {}

    template <class Itr>
    fastavector( Itr first, Itr const& last, size_type capacity = 0,
              Alloc const& alloc = Alloc() )
    : charvector(first,last,capacity,alloc) {}

    fastavector( const String & s ) { assign(s.begin(),s.end()); }

    explicit fastavector( const bvec& b ) { SetFrom(b); }

    // compiler-supplied copying and destructor are OK

    fastavector& resize( size_type n, char val = GeneralizedBase::X.asChar() )
    { charvector::resize(n,val); return *this; }

    void SetFrom( const bvec& b )
    { clear().reserve(b.size());
    bvec::const_iterator end(b.end());
    for ( bvec::const_iterator itr(b.begin()); itr != end; ++itr )
       push_back(Base::val2Char(*itr)); }

    String ToString() const
    { String s; s.Set(&(*this)[0],size()); return s; }

    /// SetToSubOf(fv, start, len):  Set this to the length len piece of fv,
    /// starting at position start.  The case where this == &fv is allowed.
    /// If len == -1, go all the way to the end of fv.
    fastavector& SetToSubOf( fastavector const& fv, size_type start,
                             size_type len )
    { AssertLe( start, fv.size() );
      if ( len == ~0U ) len = fv.size() - start;
      AssertLe( len, fv.size()-start );
      reserve(len);
      memcpy(data(),fv.data()+start,len);
      setSize(len);
      return *this; }

    bool HasGaps() const
    { const_iterator stop(end());
      for ( const_iterator itr(begin()); itr != stop; ++itr )
          if ( *itr == 'n' ) return true;
      return false; }

    /// Split this into chunks, separated by gaps ('n'), and return each chunk
    /// as a gapless fastavector.
    vec<fastavector> SplitOnGaps() const;

    /// Returns a basevector, and sets a bitvector of ambiguous bases, if supplied
    /// possibly represent.  Fails if this fastavector has any gaps.
    /// Returns an empty set if the number of possibilities is more than max.
    basevector ToBasevector( bitvector* ambiguous = NULL) const;

    /// Returns the set of all basevectors that could this fastavector could
    /// possibly represent.  Fails if this fastavector has any gaps.
    /// Returns an empty set if the number of possibilities is more than max.
    vecbasevector AllBasevectors( size_t max = 1000 ) const;

    /// Returns the set of all kmers contained in AllBasevectors (but if a kmer
    /// location has more than small_max possibilities, ignore it.)
    /// Much more likely to return a reasonably sized set than AllBasevectors.
    vecbasevector AllKmers( size_type J, size_t small_max = 100 ) const;

    /// Prints in a fasta format: "><string_id>\n" followed by the full base
    /// sequence stored in the basevector; breaks
    /// long sequences nicely into 80-character lines
    void Print( ostream& out, const String& id = "" ) const;

    void ReverseComplement();

    void Append( const fastavector& f )
    { append(f.begin(),f.end()); }

    int AmbCount() const
    { int result = 0;
      const_iterator stop(end());
      for ( const_iterator itr(begin()); itr != stop; ++itr )
          if ( GeneralizedBase::fromChar(*itr).isAmbiguous() ) result += 1;
      return result; }

  int MaxWindowAmbCount(const int sz_window) const 
  {
    int max_counts = 0;
    const int nb = size();
    if (nb < sz_window) {
      for (int ib = 0; ib < nb; ib++)
        if (GeneralizedBase::fromChar((*this)[ib]).isAmbiguous())
          max_counts++;
    }
    else {
      int counts = 0;
      for (int ib = 0; ib < sz_window; ib++)
        if (GeneralizedBase::fromChar((*this)[ib]).isAmbiguous())
          counts++;
    
      max_counts = counts;

      for (int ib = sz_window; ib < nb; ib++) {
        if (GeneralizedBase::fromChar((*this)[ib            ]).isAmbiguous()) counts++;
        if (GeneralizedBase::fromChar((*this)[ib - sz_window]).isAmbiguous()) counts--;
        if (counts > max_counts)
          max_counts = counts;
      }
    }
    return max_counts;
  }
    



    fastavector& combine( fastavector const& fv )
    { if ( fv.size() > size() ) resize(fv.size());
      iterator dst(begin());
      const_iterator stop(fv.end());
      for ( const_iterator itr(fv.begin()); itr != stop; ++itr, ++dst )
          *dst = GeneralizedBase::ambiguityCode(*itr,*dst);
      return *this; }

    /// Compute the concatenation of two fastavectors
    friend fastavector Cat( const fastavector& left, const fastavector& right )
    { fastavector result(left.begin(),left.end(),left.size()+right.size());
      result.append(right.begin(),right.end());
      return result; }

    friend fastavector Cat( const fastavector& v1, const fastavector& v2,
         const fastavector& v3 )
    {    fastavector v( v1.begin( ), v1.end( ), v1.size() + v2.size() + v3.size() );
         v.append( v2.begin( ), v2.end( ) );
         v.append( v3.begin( ), v3.end( ) );
         return v;    }

    /// Combine: merge fastavectors starting at their beginnings.
    friend fastavector Combine( const vec<fastavector>& F )
    { fastavector result;
      vec<fastavector>::const_iterator end(F.end());
      for ( vec<fastavector>::const_iterator itr(F.begin()); itr != end; ++itr )
          result.combine(*itr);
      return result; }

     friend fastavector Combine( const fastavector& f1, const fastavector& f2 )
     { fastavector result(f1); result.combine(f2); return result; }

     // Load from fasta file.

     friend void LoadFromFastaFile( const String& f, vec<fastavector>& v );
  
     friend void LoadFromFastaFile( const String& f, vec<fastavector>& v, 
				    vec<String>& names );
     

     friend void BinaryWrite( int fd, const fastavector& b );
     friend void BinaryRead( int fd, fastavector& b );

     // Write a fasta-format file which contains the bases in the
     // input fasta scaffolded together (with gaps, etc.) as defined
     // by the scaffolds.  If supplied, rc defines which contigs are
     // reverse-complement. Gaps < min_gap will be reset at min_gap.
     friend void WriteScaffoldedFasta( const String &out_file,
				       const vec<fastavector> &fasta,
				       const vec<superb> &scaffolds,
				       const vec<Bool> &rc = vec<Bool>( ),
				       const int min_gap = 1, 
				       const char gap_char = 'N',
				       const Bool ncbi_format = False);

     // Remove contigs that don't appear in a scaffold.  

     friend void RenumberAndMinimize( vec<fastavector>& f, vec<superb>& s );
  
};

SELF_SERIALIZABLE(fastavector);

typedef MasterVec<fastavector> vecfastavector;
typedef vecfastavector vecfvec;

typedef fastavector fvec;
typedef fastavector FastaVec;
typedef vecfastavector FastaVecVec;


#endif
