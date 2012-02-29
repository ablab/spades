#ifndef __INCLUDE_SortKmersImpl_H
#define __INCLUDE_SortKmersImpl_H

// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "math/Functions.h"
#include "dna/Bases.h"
#include "system/StaticAssert.h"
#include "kmers/SortKmers.h"
#include "kmers/KmerRecord.h"

/**
   File: SortKmersImpl

   Extraction of kmers from reads, and sorting of the result; see SortKmers() .

   Upon entry, reads and read_ids are given (as is pass).  There are (at
   present) passes 0 (only), passes 0 through 9, or 0 through 99,
   depending on whether Passes is 1 or 10 or 100.  In each pass, we
   extract some of the k-mers b from each read.

   Compute the reverse complement c of b, and pick the "minimum" of b, c.
   The exact definition of minimum doesn't matter, except that it is symmetric:
   min(b,c) = min(c,b).

   We then form the triple (minimum (b or c as above), read_id, k-mer position),
   as a kmer_record.  These are placed in R (in total S of them), and sorted,
   using only the k-mer part.  Thus the output values are R and S.

   We store the position as pos + 1 (if we don't use the reverse complement),
   and as -(pos + 1) if we do.

   When a kmer is palindromic two entries are created in R, one for the
   forward direction (pos+1) and one for the reverse complement (-(pos+1)).

   In the 10 pass case, the passes are defined by the following table:

   0: A..A or T..T      5: C..C or G..G
   1: C..A or T..G      6. G..C or G..C
   2. G..A or T..C      7. A..T or A..T
   3. T..A or T..A      8. A..G or C..T
   4: A..C or G..T      9. C..G or C..G

   For example, pass 0 includes all k-mers of the form A..A or T..T.

   In the 100 pass case, pass 10 would include all k-mers of the form
   AC..AA or TC..AT or AT..GA or TT..GT.  Etc.  In the 100 pass case, the variable
   'pass' may be thought of as a collection of two independent digits, each of
   which restrains the chosen kmers to a particular choice at a particular
   location.  For example:
   
   -- All pass numbers with the ones digit 2 (2, 12, 22, 32, etc.) force a kmer to
   match the form G..A or T..C.
   -- All pass numbers with the tens digit 4 (40, 41, etc.) force a kmer to match
   the form *A..C* or *G..T*, where * can be any base.
   -- Hence, pass number 42 will accept only kmers with one of the following forms:
   GA...CA, GG...TA, TA...CC, TG...CT.
   
   Memory utilization has been reduced for the case Passes=100, for huge
   "reads"; the same improvement has not been implemented for the other cases.

   \ingroup grp_kmerGathering
*/

#define SORT_CORE \
{  KSHAPE::extractKmer(b, *reads[l], q); \
   basevector b_rc; \
   b_rc.ReverseComplement( b );  \
   if ( b_rc < b  ||  b == b_rc ) R[S++].Set( b_rc, read_ids[l], -(q+1) ); \
     if ( b < b_rc  ||  b == b_rc ) R[S++].Set( b, read_ids[l], q+1 );    }


#define BASE_EQ(BYTE,SHIFT,BASE) ((BYTE & (3<<SHIFT)) == (BASE)<<SHIFT)

// 1 PASS VERSION

/**
   \copydoc SortKmers(dummy<100>,const vecbasevector&,const vec<int>&,int,vec<RECORD>&,unsigned int&,bool)

   This routine may not scale well; for a better-scaling version, see the 100-pass version:
   SortKmers(dummy<100>,const vecbasevector&,const vec<int>&,int,vec<RECORD>&,unsigned int&) .

   \callergraph
*/
template< class KSHAPE, class RECORD >
void SortKmers( dummy<1>,
                vec<bvec const*> const& reads,
		const vec<int>& read_ids,
		int,
		vec<RECORD>& R,
		unsigned int& S,
		Bool use_stable_sort,
		Bool palind_both_dirs )
{
  const int K = KSHAPE::KSIZE;
  const unsigned int KSPAN = (unsigned int)KSHAPE::KSPAN;
  basevector b(K);
  S = 0;
  for ( size_t l = 0; l < reads.size( ); l++ )
    {    if ( reads[l]->size( ) < KSPAN ) continue;
    unsigned int N = reads[l]->size( ) - KSPAN + 1;
    if ( S + 2*N >= R.size( ) ) {
      unsigned nn = Max( (long unsigned) ( 1.2 * R.size( ) ), (R.size( ) + 2*N) );
      if ( nn < R.size( ) ) FatalErr( "SortKmers<1>: Unsigned-int overflow (R=" << R.size( ) << ")" );
      R.resize(nn);
    }
    
    Bool use_b, use_c;
    for ( unsigned int q = 0; q < N; q++ )
      {    SORT_CORE    }
    }
  if (use_stable_sort)
    stable_sort( R.begin( ), R.begin( ) + S );
  else
    sort( R.begin( ), R.begin( ) + S );
}

// *********************** SIMPLE 10 PASS VERSION ************************
// Intended only for quality control, but not maintained.

/*
  template<int K, int I> void SortKmers( dummy<10>, const vecbasevector& reads,
  const vec<int>& read_ids, int pass, vec< kmer_record<K, I> >& R,
  unsigned int& S )
  {    basevector b(K), c(K);
  S = 0;
  for ( unsigned int l = 0; l < reads.size( ); l++ )
  {    if ( reads[l].size( ) < K ) continue;
  unsigned int N = reads[l].size( ) - K + 1;
  if ( S + 2*N >= R.size( ) )
  R.resize( MAX( 6*R.size( )/5, R.size( ) + 2*N ) );
  Bool use_b, use_c;

  if ( pass == 7 ) pass = 12;
  unsigned char bases[2];
  bases[0] = pass % 4;
  bases[1] = (pass >> 2);

  for ( unsigned int q = 0; q < N; q++ )
  {    if ( (reads[l](q) == bases[0] && reads[l](q+K-1) == bases[1]) ||
  (reads[l](q) == 3-bases[1] && reads[l](q+K-1) == 3-bases[0]) )
  {    SORT_CORE    }    }    }
  sort( R.begin( ), R.begin( ) + S );    }
*/

// *********************** 10 PASS VERSION ************************

/**
   \copydoc SortKmers(dummy<100>,const vecbasevector&,const vec<int>&,int,vec<RECORD>&,unsigned int&,bool)

   \note This simple 10-pass version is intended only for quality control, but not maintained.
*/
template<class KSHAPE, class RECORD>
void SortKmers( dummy<10>,
                vec<bvec const*> const& reads,
                const vec<int>& read_ids,
		int pass,
		vec<RECORD>& R,
		unsigned int& S,
		Bool use_stable_sort,
		Bool palind_both_dirs )
{
  const int K = KSHAPE::KSIZE;
  const unsigned int KSPAN = (unsigned int)KSHAPE::KSPAN;
  basevector b(K);
  S = 0;
  for ( size_t l = 0; l < reads.size( ); l++ ) {
    if ( reads[l]->size( ) < KSPAN ) continue;
    unsigned int N = reads[l]->size( ) - KSPAN + 1;

    if ( pass == 7 ) pass = 12;
    unsigned char ba = pass % 4, bb = (pass >> 2) % 4;

    Bool use_c, use_b;

    unsigned int q = 0;
    const basevector& r = *reads[l];
    while(1) {
      while( q < N ) {
        if ( r[q] == ba && r[q+K-1] == bb ) break;
        if ( r[q] == 3 - bb && r[q+K-1] == 3 - ba ) break;
        ++q;
      }
      if ( q == N ) break;
      if ( S + 2*N >= R.size( ) ) {
        unsigned nn = Max( (long unsigned) ( 1.2 * R.size( ) ), (R.size( ) + 2*N) );
        if ( nn < R.size( ) ) FatalErr( "SortKmers<10>: Unsigned-int overflow (R=" << R.size( ) << ")" );
        R.resize(nn);
      }

      SORT_CORE;
      ++q;
    }
  }
  if ( use_stable_sort)
    stable_sort( R.begin( ), R.begin( ) + S );
  else
    sort( R.begin( ), R.begin( ) + S );
}

// *********************** 100 PASS VERSION ************************

template<class KSHAPE, class RECORD>
void one_pass_of_sort( vec<bvec const*> const& reads,
		       const vector< read_id_t >& read_ids,
		       int pass1, int pass2, Bool palind_both_dirs,
		       Bool use_stable_sort,
		       vector<RECORD>& R, unsigned int& S )
{
    const nbases_t K = KSHAPE::KSIZE;
    const nbases_t KSPAN = KSHAPE::KSPAN;

    basevector b(K);

    S = 0;

    // Bases for this pass: All test kmers must match these bases or their rc
    const base_t bA = pass1 % 4;
    const base_t bB = pass2 % 4;
    const base_t bY = (pass2 >> 2) % 4;
    const base_t bZ = (pass1 >> 2) % 4;

    // Reverse complements of each base
    // Note that the RC first base is the complement of the last base, etc.
    const base_t bA_rc = GetComplementaryBase( bZ );
    const base_t bB_rc = GetComplementaryBase( bY );
    const base_t bY_rc = GetComplementaryBase( bB );
    const base_t bZ_rc = GetComplementaryBase( bA );

    for ( size_t l = 0; l < reads.size( ); l++ ) {
      const basevector& r = *reads[l];
      if ( r.isize( ) < KSPAN ) continue;

      read_id_t thisReadId = read_ids.empty() ? l : read_ids[l];
      
      const read_pos_t N = r.size( ) - KSPAN + 1;

      typedef typename basevector::const_iterator cit_t;

      read_pos_t q = 0;     /* start of kmer in the read */

      // kA = first base in kmer, kB = second base, ..., kZ = last base
      base_t kA, kB, kY, kZ;
      cit_t iter_beg = r.Begin();
      cit_t iter_end = r.Begin( KSPAN - 2 );
      kA = *( iter_beg++ );
      kB = *iter_beg;
      kY = *( iter_end++ );
      kZ = *iter_end;

      int top = 0;
      while( top < N )
      {
	  AssertEq( r[ q ], kA );
	  AssertEq( r[ q+1 ], kB );
	  AssertEq( r[ q+K-2 ], kY );
	  AssertEq( r[ q+K-1 ], kZ );

	  top += Min( N - top, 10000 );
	  if ( S + 2*(top - q) >= R.size( ) ) {
	    unsigned nn = Max( (unsigned) ( 1.2 * R.size( ) ), (S + 2*(top - q)) );
	    if ( nn < R.size( ) ) FatalErr( "SortKmers<100>: Unsigned-int overflow (R=" << R.size( ) << ")" );
	    R.resize(nn);
	  }

	  while( True ) {
	    if ( q == top ) break;

	    AssertEq( r[ q ]    , kA );
	    AssertEq( r[ q+1 ]  , kB );
	    AssertEq( r[ q+K-2 ], kY );
	    AssertEq( r[ q+K-1 ], kZ );

	    // Determine if this read belongs in this pass
	    if ( ( kA == bA     && kZ == bZ  ||
		   kA == bA_rc  && kZ == bZ_rc ) &&
		 ( kB == bB    && kY == bY  ||
		   kB == bB_rc && kY == bY_rc ) )
	      {
		KSHAPE::extractKmer(b, r, q);

		basevector::CanonicalForm canon_form = b.Canonicalize();
		if ( canon_form == basevector::rc_form ||
		     canon_form == basevector::palindrome && palind_both_dirs ) R[S++].Set( b, thisReadId, -(q+1) );
		if ( canon_form == basevector::fw_form ||
		     canon_form == basevector::palindrome ) R[S++].Set( b, thisReadId, q+1 );
	      }

	    ++q;

	    // Advance the kmer by incrementing its bases
	    kA = kB;
	    kB = *( ++iter_beg );
	    kY = kZ;
	    if ( q < N )
	      kZ = *( ++iter_end );

	    if ( q == top ) break;
          }
      }
    }

    // Sort the kmers in R
    if (use_stable_sort)
      stable_sort( R.begin( ), R.begin( ) + S );
    else
      sort( R.begin( ), R.begin( ) + S );
}



/**
   Extracts and sorts all kmers from a read set.

   \copydoc SortKmersImpl.h

   Template arguments:
      \li \c K - the size of the kmers
      \li \c RECORD - the type of kmer occurrence records to generate, to put into the output array \p R .  The type must have
              a Set method to record the kmer sequence, read id and position; see kmer_record::Set.  Two example types are
             \link kmer\endlink and kmer_record.
      \li \c KSHAPE - the \link KmerShape.h shape\endlink of the kmers to extract

   \param[in] reads the reads from which to extract the kmers
   \param[in] read_ids an array parallel to \p reads, giving an integer id to each corresponding read in \p reads.
   \param[in] pass the pass; this function must be called once for each pass, where the number of passes is indicated
       by the first argument to this function (and chooses the right function version for the given number of passes).
   \param[out] R the kmer records, one for each occurrence of a kmer in a read, referencing
       back to the reads.
   \param[out] S the number of records added to \p R by this call
   \param[in]  use_stable_sort if \c true, records for a given kmer that \link kmer_record::GetId refer\endlink to the same read will be
               kept adjacent in the output array \p R .

   \callergraph

   Implementation notes:
      \li on any invocation of this routine, exactly one call to CALL_ONE_PASS_OF_SORT() happens.  Remember that this routine must be called multiple times,
          once for each \p pass value.
*/
template<class KSHAPE, class RECORD>
void SortKmers( dummy<100>, vec<bvec const*> const& reads, const vec<int>& read_ids,
		int pass, vec<RECORD>& R, unsigned int& S, Bool use_stable_sort,
		Bool palind_both_dirs )
{
  // We separate the pass number into its constituent digits and put each one
  // in the range of [0-6],8,9,12.  These numbers encode (in base 4) the ten
  // pairs of bases that are unique modulo RC: 0 = AA, 1 = AC, 4 = CA, etc.
  int pass1 = pass % 10;
  int pass2 = pass / 10;
  if ( pass1 == 7 ) pass1 = 12;
  if ( pass2 == 7 ) pass2 = 12;

  // Fill R with a sorted list of all kmers in moreReads that match pass1 and pass2
  one_pass_of_sort<KSHAPE,RECORD>( reads, read_ids, pass1, pass2, palind_both_dirs, use_stable_sort, R, S );
  
}  // SortKmers()


#endif
// #ifndef __INCLUDE_SortKmersImpl_H
