// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


#ifndef SORT_KMERS_H
#define SORT_KMERS_H

/**
   File: SortKmers.h

   Extracting and sorting of kmers.  See <SortKmers()>.

   @file
*/

#include <iosfwd>

#include "Basevector.h"
#include "Vec.h"
#include "kmers/KmerRecord.h"
#include "kmers/KmerShape.h"

/**
   Class: dummy
   
   Used to select the right version of <SortKmers()>: 1-pass, 10-pass or 100-pass.
 */
template<int Passes> class dummy { };

/**
   Macro: SortKmersNumPasses
   
   A slightly better name for the class (<dummy>) used to statically
   select the version of <SortKmers()> that is called.
 */
#define SortKmersNumPasses dummy

/**
   Type Concept: SortKmersOutputRecord

   A record representing one occurrence of a kmer in a read.  Classes
   that model this concept can be used as the _RECORD_ template
   argument to <SortKmers()>, which extracts all kmer occurrences from
   the reads and records each occurrence in a record of this class.
   To model this type concept, a type needs to define the members
   described below.  Not all these members are used by SortKmers(),
   but we still include them in this type concept as they are used
   by related code and it's useful to have a single type concept
   that captures all the related members.
   <kmer_record> and <kmer> model this type concept.

   Constant: BASES_SIZE

   The size of the kmer stored in this record.
   
   >static const int BASES_SIZE;
   
   Method: Set

   Records the occurrence of a kmer in a read.

   >  void Set( const kmer_t& b, int read_id, int read_pos );
   
   Parameters:

      b - the kmer (its sequence)
      read_id - the <read id> of the read in which this kmer occurs
      read_pos - the position in the read of this kmer's occurrence;
        it is 1-offset, and is positive if the <kmer's canonical form>
	occurs in the read, and negative if the reverse complement
	of the kmer's canonical form occurs in the read.

   Method: Bytes

   Return the kmer data as an array of bytes.

   > const unsigned char* Bytes( ) const;

   Method: Ints

   Return the kmer data as an array of unsigned ints.

   > const unsigned int* Ints( ) const;

   Method: GetBasevector

   Sets a <basevector> to the content of this kmer.
   
   > void GetBasevector( basevector& kmer ) const;
   
*/

// End: Section

/**
   Declare two versions of SortKmers() for the given number of passes:
   the version where all template arguments need to be specified,
   and the version where all but the first two are set to defaults.
*/
#define DECLARE_SORT_KMERS(numPasses)                                          \
 template<class KSHAPE, class RECORD>                                          \
 void SortKmers ( SortKmersNumPasses<numPasses> d1,                            \
                  vec<bvec const*> const& reads,                               \
                  const vec<int>& read_ids, int pass, vec<RECORD>& R,          \
                  unsigned int& S, Bool use_stable_sort = False,               \
                  Bool palind_both_dirs = True );                              \
                                                                               \
 template<int K, int I>                                                        \
 inline void SortKmers ( SortKmersNumPasses<numPasses> d1,                     \
                  vec<bvec const*> const& reads,                               \
                  const vec<int>& read_ids, int pass,                          \
                  vec< kmer_record<K,I> >& R,                                  \
                  unsigned int& S, Bool use_stable_sort = False,               \
                  Bool palind_both_dirs = True )                               \
 { SortKmers< KmerShapeDefaultClass(K), kmer_record<K,I> >(d1,reads,read_ids,  \
   pass,R,S,use_stable_sort,palind_both_dirs); }                               \
                                                                               \
 template<class KSHAPE, class RECORD>                                          \
 inline void SortKmers ( SortKmersNumPasses<numPasses> d1,                     \
                    const vecbasevector& reads,                                \
                    const vec<int>& read_ids, int pass, vec<RECORD>& R,        \
                    unsigned int& S, Bool use_stable_sort = False,             \
                    Bool palind_both_dirs = True )                             \
 { vec<bvec const*> vvv; vvv.reserve(reads.size());                            \
   vecbvec::const_iterator end(reads.end());                                   \
   for ( vecbvec::const_iterator itr(reads.begin()); itr != end; ++itr )       \
       vvv.push_back(&*itr);                                                   \
   SortKmers<KSHAPE,RECORD>(d1,vvv,read_ids,pass,R,S,use_stable_sort,          \
                            palind_both_dirs); }                               \
                                                                               \
 template<int K, int I>                                                        \
 inline void SortKmers ( SortKmersNumPasses<numPasses> d1,                     \
                  const vecbasevector& reads,                                  \
                  const vec<int>& read_ids, int pass,                          \
                  vec< kmer_record<K,I> >& R,                                  \
                  unsigned int& S, Bool use_stable_sort = False,               \
                  Bool palind_both_dirs = True )                               \
  { SortKmers< KmerShapeDefaultClass(K), kmer_record<K,I> >(d1,reads,read_ids, \
           pass,R,S,use_stable_sort,palind_both_dirs); }                       \
                                                                               \
 typedef int eatSemicolon ## numPasses

DECLARE_SORT_KMERS(1);
DECLARE_SORT_KMERS(10);
DECLARE_SORT_KMERS(100);

const Bool USE_STABLE_SORT = True;
const Bool DONT_USE_STABLE_SORT = False;
const Bool PALIND_BOTH_DIRS = True;
const Bool PALIND_FW_ONLY = False;

/**
   Macro for simplifying explicit template instantiation of SortKmers().
   See also: FOR_ALL_K(), FOR_ALL_GAPS(), FOR_ALL_K_GAP(), 
*/
#define INSTANTIATE_SORTKMERS(_KSHAPE, RECORD, _passes) \
  template void SortKmers<_KSHAPE,RECORD> ( \
     dummy<_passes>, vec<bvec const*> const& reads, const vec<int>& read_ids, \
     int pass, vec<RECORD>& R, unsigned int& S, Bool use_stable_sort, Bool palind_both_dirs )

#define INSTANTIATE_SORTKMERS_FOR_KSHAPE(KSHAPE, dummy) \
    INSTANTIATE_SORTKMERS(KSHAPE, KmerRecordType(KSHAPE, 1), 1); \
    INSTANTIATE_SORTKMERS(KSHAPE, KmerRecordType(KSHAPE, 2), 1); \
    INSTANTIATE_SORTKMERS(KSHAPE, KmerRecordType(KSHAPE, 1), 10); \
    INSTANTIATE_SORTKMERS(KSHAPE, KmerRecordType(KSHAPE, 2), 10); \
    INSTANTIATE_SORTKMERS(KSHAPE, KmerRecordType(KSHAPE, 1), 100); \
    INSTANTIATE_SORTKMERS(KSHAPE, KmerRecordType(KSHAPE, 2), 100); \
    INSTANTIATE_SORTKMERS(KSHAPE, kmer<KSHAPE::KSIZE>, 100)

#define TYPE_NAME_I_K(prefix1, prefix2, I, K) prefix1 ## _ ## prefix2 ## _ ## I ## _ ## K ## _t

/**
   Macro: INSTANTIATE_SORTKMERS_FOR_I_K

   Instantiate the SortKmers routine for the <default kmer form> and the specified
   kmer size.  Useful if a particular program needs to use kmer sizes outside the
   <supported kmers>.  See <FindMatches.cc> for an example.

   *NOTE*: before invoking this macro in your .cc file, make sure to include the following
   header files:
   
  >#include "SortKmersImpl.h"
  >#include "VecTemplate.h"
*/
#define INSTANTIATE_SORTKMERS_FOR_I_K(I, K, prefix)                            \
    typedef KmerShapeDefaultClass(K) TYPE_NAME_I_K(prefix, kmer_shape, I, K);  \
    typedef kmer_record<K,I> TYPE_NAME_I_K(prefix, kmer_record, I, K);         \
    INSTANTIATE_SORTKMERS(TYPE_NAME_I_K(prefix, kmer_shape, I, K),             \
			  TYPE_NAME_I_K(prefix, kmer_record, I, K), 1);        \
    INSTANTIATE_SORTKMERS(TYPE_NAME_I_K(prefix, kmer_shape, I, K),             \
			  TYPE_NAME_I_K(prefix, kmer_record, I, K), 10);       \
    INSTANTIATE_SORTKMERS(TYPE_NAME_I_K(prefix, kmer_shape, I, K),             \
			  TYPE_NAME_I_K(prefix, kmer_record, I, K), 100)

extern void (*SORT_KMERS_END_PASS)( int pass );

/**
    Macros: FOR_KMER_OCCS_BEG / FOR_KMER_OCCS_END

    Macros that for each kmer, let you do something for all occurrences of that
    kmer.  The occurrences are presented as a pair of iterators specifying
    a range of kmer records.

    Synopsis:

    (begin example)
    
    typedef kmer_record< K, 2 > mykrec_t;

    FOR_KMER_OCCS_BEG( KSHAPE, mykrec_t, reads, false, from, to ) {
    
       // here, from and to are of type vec< mykrec_t >::const_iterator.
       // the range [from, to) specifies the range of kmer records,
       // of type mykrec_t, that give the occurrences of one kmer in the reads.

       // Standard C++ looping constructs (break and continue) may be used.
       
    }  // FOR_KMER_OCCS_BEG()
    FOR_KMER_OCCS_END( );
    
    (end example)

    For an example of use see MarkTrustedA1.cc .
 */
#define FOR_KMER_OCCS_BEG2(KSHAPE, KRECORD, reads, use_stable_sort, palind_both_dirs, BEG_PASS, FROM, TO)  do { \
    const vecbasevector& _seqs = reads; \
    \
    int _N = _seqs.size( ); \
    const int _K = KSHAPE::KSIZE; \
    vec<read_id_t> _rid; \
    unsigned int _S = 0; \
    for ( size_t l = 0; l < _seqs.size( ); l++ ) \
      _S += _seqs[l].size( ) - _K + 1; \
    _S += _S/4; \
    _S /= 33; \
    vec< KRECORD > _R(_S); \
    do { BEG_PASS; } while( 0 );		  \
    for ( int _pass = 0; _pass < 100; _pass++ ) { \
      SortKmersNumPasses<100> _d100; \
      SortKmers< KSHAPE, KRECORD >( _d100, _seqs, _rid, _pass, _R, _S, use_stable_sort, palind_both_dirs ); \
      int _j = 0; \
      for ( int _i = 0; _i < (int) _S; _i = _j ) { \
	for ( _j = _i+1; _j < (int) _S; _j++ ) { \
	  if ( _R[_i] < _R[_j] ) break; \
	} \
	typename vec< KRECORD >::const_iterator FROM = _R.begin() + _i; \
	typename vec< KRECORD >::const_iterator TO = _R.begin() + _j;

#define FOR_KMER_OCCS_END2(END_PASS) \
      } \
  do { END_PASS; } while( 0 );			\
    } \
  } while( 0 )


#define FOR_KMER_OCCS_BEG(KSHAPE, KRECORD, reads, use_stable_sort, palind_both_dirs, FROM, TO) \
  FOR_KMER_OCCS_BEG2(KSHAPE, KRECORD, reads, use_stable_sort, palind_both_dirs, cout << "pass " << flush, FROM, TO )

#define FOR_KMER_OCCS_END(DUMMY) \
  FOR_KMER_OCCS_END2( if ( SORT_KMERS_END_PASS ) (*SORT_KMERS_END_PASS)(_pass ); )







#endif
//#ifndef SORT_KMERS_H
