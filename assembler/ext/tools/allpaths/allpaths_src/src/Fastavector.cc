///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Fastavector.h"
#include "Bitvector.h"
#include "feudal/OuterVecDefs.h"
#include "graph/DigraphTemplate.h"
#include "system/System.h"
#include <cstring>
#include <istream>
#include <string>

// Split this into chunks, separated by gaps ('n'), and return each chunk as a
// gapless fastavector.
// TODO: generalize this into a templatized STL algorithm.
vec<fastavector>
fastavector::SplitOnGaps() const
{
    vec<fastavector> chunks;
    const_iterator stop(end());
    for ( const_iterator itr(begin()); itr != stop; ++itr )
    {
        if ( *itr != 'n' ) // found a chunk's start
        {
            const_iterator itr2(itr);
            while ( ++itr2 != stop && *itr2 != 'n' )
                ; // just keep going until we find the chunk's end

            chunks.push_back(fastavector(itr,itr2));

            if ( itr2 == stop )
                break;
            itr = itr2;
        }
    }
    return chunks;
}

// Returns a basevector, and sets a bitvector of ambiguous bases, if supplied.
basevector 
fastavector::ToBasevector( bitvector* ambiguous) const
{
  basevector bv(size());
  if (ambiguous) {
    ambiguous->resize(size());
    ambiguous->Zero();
  }

  // loop over bases
  for ( size_type idx = 0; idx < size(); ++idx ) {
    GeneralizedBase const& gb = GeneralizedBase::fromChar((*this)[idx]);
    char const* endBase = gb.end();
    char const* curBase = gb.bases();

    // set basevec base to first possible base code
    bv.set(idx, Base::char2Val(*curBase));

    // mark if this base is ambiguous
    if (endBase - curBase > 1 && ambiguous)
      ambiguous->set(idx, 1);
  }
  return bv;
}



// Returns the set of all basevectors that this fastavector could
// possibly represent.  Fails if this fastavector has any gaps.
// Returns an empty set if the number of possibilities is more than max.
vecbasevector
fastavector::AllBasevectors( size_t maxVecs ) const
{
  ForceAssert( !HasGaps() );
  vecbasevector all;
  
  size_t nVecs = 1;
  const_iterator stop(end());
  for ( const_iterator itr(begin()); itr != stop; ++itr ) {
    GeneralizedBase const& gb = GeneralizedBase::fromChar(*itr);
    nVecs *= gb.end() - gb.bases();
    if ( nVecs > maxVecs )
      return all; // too many options - return an empty vecbasevector
  }
  
  all.reserve(nVecs);
  
  // Create an initial copy -- ambiguity codes will be squashed randomly,
  // which is OK.  We'll nail 'em down later.
  all.push_back(bvec(begin(),end(),GenCharToRandomBaseMapper()));
  
  // If there is no ambiguity anywhere, we can return the single vecbasevector.
  if ( nVecs == 1 ) return all;
  
  // Find each ambiguous position, with more than one base option
  for ( size_type idx = 0; idx < size(); ++idx ) {
    
    GeneralizedBase const& gb = GeneralizedBase::fromChar((*this)[idx]);
    char const* endBase = gb.end();
    char const* curBase = gb.bases();
    if ( endBase - curBase == 1 ) continue; // only one base option here
    
    // take the current set of bvecs and set the code in the
    // ambiguous position to the first legal code
    bvec::value_type code = Base::char2Val(*curBase);
    vecbvec::iterator end(all.end());
    for ( vecbvec::iterator itr(all.begin()); itr != end; ++itr )
      itr->set(idx,code);
    
    // for the remaining legal codes in the ambiguity set at this idx
    for ( ++curBase; curBase < endBase; ++curBase ) {
      // append a copy of the current set, as it was initially,
      // for each additional legal base, and set the code in the
      // ambiguous position to the next legal value
      code = Base::char2Val(*curBase);
      
      for ( vecbvec::iterator itr(all.begin()); itr != end; ++itr ) {
	// this could theoretically invalidate the iterators,
	// but it doesn't because we've reserved enough space
	all.push_back(*itr);
	all.back().set(idx,code);
      }
    }
  }
  
  return all;
}

// Returns the set of all kmers contained in AllBasevectors (but if a kmer
// location has more than small_max possibilities, ignore it.)
// Much more likely to return a reasonably sized set than AllBasevectors.
vecbasevector
fastavector::AllKmers( fastavector::size_type K, size_t small_max ) const
{
  vecbasevector result;
  if ( size() < K ) return result; // empty result
  
  size_type nKmers = size() - K + 1;
  result.reserve(nKmers);
  
  const_iterator stop(begin(nKmers));
  const_iterator itr2(begin(K));
  
  bool unambig = false;
  
  for ( const_iterator itr(begin()); itr != stop; ++itr, ++itr2 ) {
    
    // As a time-saving step, we use the unambig flag to keep track of when
    // the fastavector contains no ambiguous bases.  In these areas we can
    // create the basevector directly from the fastavector, instead of calling
    // AllBasevectors.
    int n_ambig;
    if ( unambig ) {
      result.push_back( basevector( itr, itr2, GenCharToRandomBaseMapper() ) );
      n_ambig = 1;
    }
    else {
      vecbasevector vbv(0);
      vbv = fastavector(itr,itr2).AllBasevectors(small_max);
      result.append( vbv.begin(), vbv.end() );
      n_ambig = vbv.size();
    }
    
    // Update the unambig flag.
    unambig = false;
    if ( n_ambig == 1 )
      if ( itr2 != end() ) {
	GeneralizedBase const& gb = GeneralizedBase::fromChar(*itr2);
	if ( gb.end() - gb.bases() == 1 )
	  unambig = true;
      }
  }
  
  return result;
}

void fastavector::Print( ostream& out, const String& id ) const
{
    out << '>' << id;
    for ( size_type i = 0; i < size(); ++i )
    {
        if ( !(i % 80) ) out << '\n';
        out << (*this)[i];
    }
    out << '\n';
}

void fastavector::ReverseComplement()
{
    GeneralizedBase::reverseComplement(begin(), end());
}

void BinaryWrite( int fd, const fastavector& b )
{
    fastavector::size_type n = b.size();
    WriteBytes(fd, &n, sizeof(n));
    if ( n )
        WriteBytes(fd, &b[0], b.size());
}

void BinaryRead( int fd, fastavector& b )
{
    fastavector::size_type n;
    ReadBytes(fd, &n, sizeof(n));
    b.resize(n);
    if ( n > 0 )
        ReadBytes(fd, &b[0], b.size());
}

void LoadFromFastaFile( const String& f, vec<fastavector>& vfv )
{
    ifstream in(f.c_str());
    std::string buf;
    buf.reserve( 8192 );
    fastavector fv;
    fv.reserve( 100000 );

    vfv.clear();

    bool dataPresent = false;
    while ( getline(in,buf) )
    {
        size_t sz = buf.size();
        if ( !sz )
            continue; // NON-STRUCTURED!  ignore blank lines

        if ( buf[0] != '>' ) // if not a FASTA comment
        {
            if ( fv.size() + sz > fv.capacity() )
                // grow rapidly to avoid a realloc per line on whole-genome
                // files (which typically have a gigantic chr1 coming first)
                fv.reserve( 4*fv.capacity() );

            // buffer the FASTA
            fv.append(buf.begin(),buf.end());
        }
        else if ( dataPresent )
        {
            vfv.push_back(fv);
            fv.clear();
        }
        dataPresent = true;
    }

    if ( !dataPresent )
    {
        if ( !in.eof() )
            FatalErr("Unable to read FASTA file " << f);
    }
    else
    {
        vfv.push_back(fv);
        if ( !in.eof() )
            FatalErr("Unable to read FASTA file " << f << " to completion");
    }
}

void LoadFromFastaFile( const String& f, vec<fastavector>& vfv, vec<String>& vnames )
{
    ifstream in(f.c_str());
    std::string buf;
    buf.reserve( 8192 );
    fastavector fv;
    fv.reserve( 100000 );

    vfv.clear();

    bool dataPresent = false;
    while ( getline(in,buf) )
    {
        size_t sz = buf.size();
        if ( !sz )
            continue; // NON-STRUCTURED! ignore blank lines

        if ( buf[0] != '>' ) // if not a FASTA comment
        {
            if ( fv.size() + sz > fv.capacity() )
                // grow rapidly to avoid a realloc per line on whole-genome
                // files (which typically have a gigantic chr1 coming first)
                fv.reserve(4 * fv.capacity());

            // buffer the FASTA
            fv.append(buf.begin(), buf.end());

            // if we're jumping right into FASTA without having seen a comment
            if ( !dataPresent )
                vnames.push_back("");
        }
        else
        {
            String name(buf.begin()+1,buf.end());
            DeleteTrailingWhiteSpace(name);
            vnames.push_back(name);
            if ( dataPresent )
            {
                vfv.push_back(fv);
                fv.clear();
            }
        }
        dataPresent = true;
    }

    if ( !dataPresent )
    {
        if ( !in.eof() )
            FatalErr("Unable to read FASTA file " << f);
    }
    else
    {
        vfv.push_back(fv);
        if ( !in.eof() )
            FatalErr("Unable to read FASTA file " << f << " to completion");
    }

    ForceAssertEq( vfv.size(), vnames.size() );
}

void WriteScaffoldedFasta( const String &out_file,
			   const vec<fastavector> &fasta,
			   const vec<superb> &scaffolds,
			   const vec<Bool> &rc,
			   const int min_gap, 
			   const char gap_char,
			   const Bool ncbi_format)
{
  Ofstream( out, out_file );
  
  // Loop over all superbs (scaffolds).
  for ( size_t i = 0; i < scaffolds.size(); i++ ) {
    const superb& S = scaffolds[i];
    if (ncbi_format) {	  // one-based, zero padding for lexical order
      out << ">scaffold";
      out.width(5);
      out.fill('0');
      out << i+1 << "\n";
    } else {
      out << ">scaffold_" << i << "\n";
    }
    
    // Find the base sequence of this scaffold, including 'n's for
    // gaps. WARNING: gaps < min_gap will be reset at min_gap.
    vec<char> s;
    for ( int j = 0; j < S.Ntigs( ); j++ ) {
      bool need_rcing = rc.nonempty() && rc[ S.Tig(j) ];
      fastavector b_rc;
      if ( need_rcing ) {
	b_rc = fasta[S.Tig(j)];
	b_rc.ReverseComplement( );
      }
      const fastavector &b = need_rcing ? b_rc : fasta[S.Tig(j)];
      for ( unsigned int l = 0; l < b.size( ); l++ )
	s.push_back( b[l] );
      if ( j < S.Ntigs( ) - 1 )
	s.push_back_copies( gap_char, Max( min_gap, S.Gap(j) ) );
    }
    
    // Print the bases, adding line breaks where necessary.
    int printed = 1;
    for ( unsigned int j = 0; j < s.size( ); j++, printed++ ) {
      out << s[j];
      if ( printed % 80 == 0 ) out << "\n";
    }
    if ( printed == 1 || printed % 80 != 1 ) out << "\n";
  }
  
  out.close( );
}  

template int digraphE<fastavector>::EdgeObjectIndexByIndexTo(int, int) const;
template void digraphE<fastavector>::DeleteEdges(vec<int> const&);
template void digraphE<fastavector>::AddVertices(int);
template vec<fastavector> digraphE<fastavector>::EdgeObjectsBetween(int, int) const;
template vec<int> digraphE<fastavector>::EdgesBetween(int, int) const;
template const fastavector& digraphE<fastavector>::EdgeObject(int) const;

void RenumberAndMinimize( vec<fastavector>& f, vec<superb>& s )
{    vec<Bool> used( f.size( ), False );
     for ( size_t i = 0; i < s.size( ); i++ )
     {    for ( int j = 0; j < s[i].Ntigs( ); j++ )
          {    int m = s[i].Tig(j);
               ForceAssertLt( m, f.isize( ) );
               used[m] = True;    }    }
     vec<int> to_new( f.size( ), -1 );
     size_t count = 0;
     for ( size_t i = 0; i < f.size( ); i++ )
     {    if ( used[i] ) 
          {    f[count] = f[i];
               to_new[i] = count++;    }    }
     f.resize(count);
     for ( size_t i = 0; i < s.size( ); i++ )
     {    for ( int j = 0; j < s[i].Ntigs( ); j++ )
          {    int m = s[i].Tig(j);
               s[i].SetTig( j, to_new[m] );    }    }    }

template class OuterVec<fastavector>;
