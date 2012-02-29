///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LOOKUP_KMER_INDEX_H
#define LOOKUP_KMER_INDEX_H

#include "CoreTools.h"
#include "Basevector.h"

/// \file KmerIndex.h: Utilities for working with kmers using
/// straightforward index encoding.  Currently no .cc file.

inline unsigned int Index(const vec<char>& bases, unsigned int pos, unsigned int K)
{   unsigned int index = 0;
    while ( K-- )
        index = (index << 2) | Base::char2Val(bases[pos++]);
     return index;    }


/// Converts Kmer index back to the base sequence, to which this index corresponds.
/// @param bases[out] the bases will be written into this basevector
/// @param index Kmer index (numerical representation) to be converted to base sequence
/// @param K size of the Kmer
/// @return returns the reference to the same \c bases vector passed as the first argument
inline basevector &  KmerIndexToBases(basevector &bases, unsigned int index, unsigned int K) {
  bases.Setsize(K);

  for ( int l = K-1 ; l >=0 ; l-- ) {
    //    char base;
    //    switch (Kmer & 0x3) {
    //    case BASE_A: base = 'A';
    bases.Set(l,(unsigned char)index & 0x3);
    index >>=2;
  }
  return bases;
}


/// Returns numeric representation (index) of a Kmer of length \c K that starts at position
/// \c pos in the sequence \c bases. This method is \em unchecked: if pos+K runs out of
/// basevector boundary, the result is undefined!
inline unsigned int Index(const basevector& bases, const unsigned int pos, const unsigned int K)
{    unsigned int index = bases[pos];
     for ( unsigned int l = pos + 1; l < pos + K; l++ )
     {
       index <<= 2;
       index ^= bases[ l ];
     }
     return index;
}

/// Returns numeric representation (index) of a Kmer of length \c K that starts at specified
/// basevector iterator. The method is \em unchecked: if the iterator can not be legitimately
/// dereferenced and then advanced and dereferenced K-1 more times, the result is undefined
template <class ITER>
inline unsigned int Index(const ITER & b_iter, const unsigned int K) {
  ITER _iter(b_iter);
  unsigned int index = 0;
  for ( unsigned int l = 0 ; l < K; l++ ) {
       index <<= 2;
       index ^= _iter.Next();
  }
  return index;
}


/// Calculates new index value from the one for the previous kmer in the basevector
inline void NextIndex(unsigned int& index, const basevector& bases,
		      const unsigned int pos, const unsigned int K)
{    index <<= 2;
     index ^= bases[pos + K - 1];
     index &= (1 << K * 2 ) - 1;    }

/// Somewhat faster way to calculate new index value from the one for
/// the previous kmer in the basevector.  If index is initialized with
/// Index(bases, 0, K) then the first value of nextpos for NextIndex()
/// should be K.  Kmask should be KmerBitmask(K).
inline void NextIndex2(unsigned int& index, const basevector& bases,
		       const unsigned int nextpos, const unsigned int Kmask)
{
  index <<= 2;
  index ^= bases[nextpos];
  index &= Kmask;
}

/// The unsigned int that has 1s in the positions used for Kmer
/// numbers, and 0s elsewhere.
inline unsigned int KmerBitmask(const unsigned int K)
{
  return (1 << K*2 ) - 1;
}

/// Wrap NextIndex2 to be more convenient to use
struct KmerIndexSeq {
  unsigned int Kmask_, K_;
  unsigned int index_, nextpos_, tmp_;

  KmerIndexSeq(unsigned int K) : Kmask_(KmerBitmask(K)), K_(K) { }
  /// Initialize for new basevector starting at position <pos>
  void Reset(const basevector &bases, unsigned int pos = 0)  {
    index_ = Index(bases, pos, K_);
    nextpos_=pos+K_;
  }

  /// Returns the current kmer index (and moves on to next)
  unsigned int operator()(const basevector &bases)
  {
    tmp_=index_;
    if (nextpos_<bases.size()) {
      NextIndex2(index_, bases, nextpos_, Kmask_);
      ++nextpos_;
    }
    return tmp_;
  }

};

/// This is an adaptor-like class that can be wrapped around
/// any basevector iterator; advancing/dereferencing this iterator
/// will give numeric representation (index) of the subsequent Kmers
/// read from the underlying basevector. Note: current implementation
/// gives a Java-like rather than an STL-like iterator: HasNext() and Next()
/// methods are used for traversal, and there is no such thing as "end()" iterator
/// to compare a running iterator to.
template <class ITER>
class KmerIterator {
 public:
   /// Constructor: initializes the iterator
   KmerIterator(unsigned int K, const ITER & iter, const ITER & iter_end) :
     iter_(iter), iter_end_(iter_end), Kmask_(KmerBitmask(K)) {

         index_ = 0;
         // note: we retrieve K-1 bases in the loop below; when
	 // next() is invoked, it will always retrieve one additional
	 // base
         for ( unsigned int l = 1 ; l < K ; l++ ) {
	     if ( iter_ == iter_end_ ) {
	         FatalErr("Basevector iterator out of bound");
	     }
	     index_ <<= 2;
	     index_ ^= iter_.Next();
	 }
	 // the locally stored iterator already points to the next base - the one
	 // right after the last base used to initialize the first Kmer index
   };


   /// Constructor: initializes the iterator from a basevector; this is
   /// a shortcut - generic, iterator-based constructor is much more flexible
   /// since it can work with different traversal directions and adapters
   KmerIterator(unsigned int K, const basevector & b) :
     iter_(b.Begin()), iter_end_(b.End()), Kmask_(KmerBitmask(K)) {
         index_ = 0;
         // note: we retrieve K-1 bases in the loop below; when
	 // next() is invoked, it will always retrieve one additional
	 // base
	 for ( unsigned int l = 1 ; l < K; l++ ) {
	     if ( iter_ == iter_end_ ) {
	         FatalErr("Basevector iterator out of bound");
	     }
	     index_ <<= 2;
	     index_ ^= *iter_;
	     ++iter_;
	 }
	 // the locally stored iterator already points to the next base - the one
	 // right after the last base used to initialize the first Kmer index
   };

   /// returns true if more Kmers are available from the basevector we are
   /// currently scanning
   Bool HasNext() { return iter_ != iter_end_ ; }

   /// Returns next Kmer index from the basevector we are scanning; the very
   /// first invocation of this method retirns Kmer at position 0 in the basevector
   /// this iterator was initialized with (if basevector was used in constructor) or
   /// at the position pointed to by the iterator (if iterator was passed to the
   /// constructor); next invocation returns Kmer starting at the next position, etc.
   unsigned int Next() {
     index_ <<= 2;
     index_ ^= *iter_;
     ++iter_;
     index_ &= Kmask_;
     return index_;
   }

 private:
  ITER iter_;
  ITER iter_end_;
  unsigned int Kmask_;
  unsigned int index_;
};


/// Helper function (similar to e.g. back_inserter in STL - makes compiler deduce
/// the return type and may save typing in some occasions).
template <class ITER>
KmerIterator<ITER> kmer_iterator(unsigned int K, const ITER & start, const ITER & end) {
  return KmerIterator<ITER>(K,start,end);
}


#endif
