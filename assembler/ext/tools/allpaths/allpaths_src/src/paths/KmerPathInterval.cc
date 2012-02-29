// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "paths/KmerPathInterval.h"
#include "VecTemplate.h"

#include <algorithm> // upper_bound

const unsigned int tagged_rpint::LOOKBACK_MAX;
const unsigned int tagged_rpint::POSITION_MAX;
const unsigned int tagged_rpint::LENGTH_MAX;
const unsigned int big_tagged_rpint::LOOKBACK_MAX;
const unsigned int big_tagged_rpint::POSITION_MAX;
const unsigned int big_tagged_rpint::LENGTH_MAX;
const unsigned int new_tagged_rpint::LOOKBACK_MAX;
const unsigned int new_tagged_rpint::POSITION_MAX;
const unsigned int new_tagged_rpint::LENGTH_MAX;
const int KmerPathInterval::maxDiff;

// Methods and friends of KmerPathInterval:

template<class TAG> void KmerPathInterval::AppendToDatabase( 
     vec<TAG>& segs, int i, int j ) const
{    if ( !isGap( ) )
     {    ForceAssertLe( (unsigned int) Length( ), TAG::LENGTH_MAX );// class-appropriate limit
          ForceAssertLe( uint(j), TAG::POSITION_MAX );// class-appropriate limit
          TAG t( Start( ), Length( ), i, j );
          segs.push_back(t);    }    }

template void KmerPathInterval::AppendToDatabase( vec<tagged_rpint>& segs, 
						  int i, int j ) const;
template void KmerPathInterval::AppendToDatabase( vec<big_tagged_rpint>& segs, 
						  int i, int j ) const;
template void KmerPathInterval::AppendToDatabase( vec<new_tagged_rpint>& segs, 
						  int i, int j ) const;

// Friend functions of tagged_rpint:

// Function: Prepare
// Given a vec<tagged_rpint>, this sorts it and sets the lookback
// entries for all of its elements.  segs[j-lookback] is the index
// of the first segment which intersects segs[j].
template<class TAG> void Prepare( vec<TAG>& segs ) {

  if ( segs.empty() ) return;

  // Note that the ONLY sort criterion is the starting kmer.
  // Hence this sort is not stable.
  Sort(segs);

  // Find lookback values.  Cap the lookback at LOOKBACK_MAX.
  typename vec<TAG>::iterator i = segs.begin(), j = segs.begin();
  while( j != segs.end() ) {
    if( i->Stop() < j->Start() )
      i++;
    else {
      ForceAssert( j >= i );
      uint lookback = j - i;
      lookback = Min( lookback, TAG::LOOKBACK_MAX ); // class-appropriate limit
      j->SetLookback( lookback );
      j++;
    }
  }
}

template void Prepare( vec<tagged_rpint>& segs );
template void Prepare( vec<big_tagged_rpint>& segs );
template void Prepare( vec<new_tagged_rpint>& segs );

/**
   Function: Contains

   Given a kmer, find all occurrences of that kmer in the reads.
   More specifically, given a kmer and a vector of <tagged read path intervals> that
   has been <Prepare()>'ed, return the list of indices in that vector
   of the tagged read path intervals containing the kmer.
*/
template<class TAG> void Contains( const vec<TAG>& segs, kmer_id_t index, 
               vec<longlong>& answer, bool append, int cap )
{    
  if ( ! append )
    answer.clear( );

  // Find index of first interval that starts with a kmer greater than the given kmer.
  TAG target;
  target.data1_ = index << 24;
  longlong to = upper_bound( segs.begin( ), segs.end( ), target ) - segs.begin( );

  // If it's the first interval, there are no intervals containing the given kmer.
  if ( to == 0 ) return;

  // Otherwise, decrement "to" so it points to the last interval that
  // starts with a kmer less than or equal to the given kmer.
  --to;

  // Use lookback to find the first interval that contains the "to"
  // interval's last kmer (which may be further back than we really
  // need to go).
  uint lookback = segs[to].Lookback();
  longlong from = to - lookback;
  
  // NOTE: The TAG data structure sets a cap on lookback size: LOOKBACK_MAX.
  // If lookback == LOOKBACK_MAX, then we need to recurse to make sure we are
  // going sufficiently far back.
  while ( lookback == TAG::LOOKBACK_MAX ) {
    lookback = segs[from].Lookback();
    from -= lookback;
  }
  

  // Copy all the indices between "from" and "to" that contain the given kmer.
  if ( cap < 0 )
  {    for ( longlong i = from; i <= to; i++ )
       {
         if ( index < segs[i].Start( ) ) continue;
         if ( segs[i].Stop( ) < index ) continue;
         answer.push_back(i);
       }
  }
  else
  {    int count = 0;
       for ( longlong i = from; i <= to; i++ )
       {
         if ( index < segs[i].Start( ) ) continue;
         if ( segs[i].Stop( ) < index ) continue;
         answer.push_back(i);
         if ( ++count == cap ) break;
       }
  }
}

template void Contains( const vec<tagged_rpint>& segs, kmer_id_t index, 
               vec<longlong>& answer, bool append, int cap );
template void Contains( const vec<big_tagged_rpint>& segs, kmer_id_t index, 
               vec<longlong>& answer, bool append, int cap );
template void Contains( const vec<new_tagged_rpint>& segs, kmer_id_t index, 
               vec<longlong>& answer, bool append, int cap );

// Overload Contains() with a second version that looks for all intervals overlapping
// a given KmerPathInterval, instead of a single kmer.

template<class TAG> void Contains( const vec<TAG>& segs, KmerPathInterval rpi, 
     vec<longlong>& answer, bool append, int cap ) 
{
  if ( ! append )
    answer.clear( );

  // Find index of first interval that starts with a kmer greater than
  // the interval's Stop.
  TAG target;
  target.data1_ = rpi.Stop() << 24;
  longlong to = upper_bound( segs.begin(), segs.end(), target ) - segs.begin();

  // If it's the first interval, there are no intervals that overlap
  // the given interval.
  if ( to == 0 ) return;

  // Otherwise, decrement "to" so it points to the last interval that
  // starts with a kmer less than or equal to the given interval's
  // Stop.
  --to;

  // Find the index of the first interval that starts with a kmer
  // greater than the given interval's Start.
  target.data1_ = rpi.Start() << 24;
  longlong from = upper_bound( segs.begin(), segs.begin()+to, target ) - segs.begin();

  // If "from" is not the first interval, decrement "from" so it
  // points to the last interval that starts with a kmer less than or
  // equal to the given interval's Start.  Use lookback to find the
  // first interval that contains containing that kmer.
  if ( from > 0 ) {
    --from;
    
    uint lookback = segs[from].Lookback();
    from -= lookback;
    
    // NOTE: The TAG data structure sets a cap on lookback size: LOOKBACK_MAX.
    // If lookback == LOOKBACK_MAX, then we need to recurse to make sure we are
    // going sufficiently far back.
    while ( lookback == TAG::LOOKBACK_MAX ) {
      lookback = segs[from].Lookback();
      from -= lookback;
    }
  }
  
  ForceAssertGe( from, 0 );

  // Copy all the indices between "from" and "to" that overlap the given KPI.
  if( cap < 0 ) {
    for( longlong i = from; i <= to; i++ ) {
      if( segs[i].Overlaps(rpi) )
	answer.push_back(i);
    }
  }
  else {
    int count=0;
    for( longlong i = from; i <= to; i++ ) {
      if( segs[i].Overlaps(rpi) ) {
	answer.push_back(i);
	if( ++count == cap ) break;
      }
    }
  }
}

template void Contains( const vec<tagged_rpint>& segs, KmerPathInterval rpi, 
     vec<longlong>& answer, bool append, int cap );
template void Contains( const vec<big_tagged_rpint>& segs, KmerPathInterval rpi, 
     vec<longlong>& answer, bool append, int cap );
template void Contains( const vec<new_tagged_rpint>& segs, KmerPathInterval rpi, 
     vec<longlong>& answer, bool append, int cap );

// This will efficiently find a single instance of the requested kmer.  
// Intended for base lookup, where you don't need to find all of them.

template<class TAG> longlong Instance( const vec<TAG>& segs, kmer_id_t k ) {

  TAG target;
  target.data1_ = k << 24;

  // Find the first interval that starts with a kmer greater than k
  typename vec<TAG>::const_iterator i = 
    upper_bound( segs.begin(), segs.end(), target );

  // If it's the first interval, there are no intervals containing k
  if ( i == segs.begin() ) return -1;
  
  // Otherwise, decrement to point to the last interval that starts
  // with a kmer <=k
  i--;

  // Use lookback to find the first interval that reaches past
  // the current i's last kmer.
  uint lookback = i->Lookback();
  i -= lookback;
  
  // NOTE: The TAG data structure sets a cap on lookback size: LOOKBACK_MAX.
  // If lookback == LOOKBACK_MAX, then we need to recurse to make sure we are
  // going sufficiently far back.
  while ( lookback == TAG::LOOKBACK_MAX ) {
    lookback = i->Lookback();
    i -= lookback;
  }
  

  // Now search forwards to find an interval containing k.
  // Usually i itself will work, but not always.
  for( ; i->Start() <= k && i < segs.end( ); i++ )
    if( i->Stop() >= k )
      return (i-segs.begin());

  // If we exit this loop without returning, then k wasn't found.
  return -1;
}

template longlong Instance( const vec<tagged_rpint>& segs, kmer_id_t k );
template longlong Instance( const vec<big_tagged_rpint>& segs, kmer_id_t k );
template longlong Instance( const vec<new_tagged_rpint>& segs, kmer_id_t k );


BINARY2_DEF(tagged_rpint);
BINARY2_DEF(big_tagged_rpint);
BINARY2_DEF(new_tagged_rpint);

BINARY3_DEF(tagged_rpint);
BINARY3_DEF(big_tagged_rpint);
BINARY3_DEF(new_tagged_rpint);

