// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#include "paths/KmerPathDatabase.h"

#include <map>

// Methods of class KmerPathDatabaseTemplate.

template <class TAG>
KmerPathDatabaseTemplate<TAG>::KmerPathDatabaseTemplate()
  : mp_taggedRpints( new vec<TAG>() ),
    m_ownRpints( true ),
    m_to( 0 ),
    m_hits( 0 ),
    m_misses( 0 ),
    m_missDistance( 0 )
{
}

template <class TAG>
void KmerPathDatabaseTemplate<TAG>::CopyFrom( const KmerPathDatabaseTemplate<TAG>& other )
{
  if ( other.m_ownRpints )
  {
    mp_taggedRpints = new vec<TAG>( *(other.mp_taggedRpints) );
    m_ownRpints = true;
  }
  else
  {
    mp_taggedRpints = other.mp_taggedRpints;
    m_ownRpints = false;
  }
  
  m_to = 0;
  m_hits = 0;
  m_misses = 0;
  m_missDistance = 0;
}

template <class TAG>
KmerPathDatabaseTemplate<TAG>::KmerPathDatabaseTemplate( const String& filename )
  : mp_taggedRpints( 0 ),
    m_ownRpints( false ),
    m_to( 0 ),
    m_hits( 0 )
{
  this->Read( filename );
}

template <class TAG>
KmerPathDatabaseTemplate<TAG>::KmerPathDatabaseTemplate( const KmerPathDatabaseTemplate<TAG>& other )
  : mp_taggedRpints( 0 ),
    m_ownRpints( false ),
    m_to( 0 ),
    m_hits( 0 ),
    m_misses( 0 ),
    m_missDistance( 0 )
{
  this->CopyFrom( other );
}

template <class TAG>
KmerPathDatabaseTemplate<TAG>::KmerPathDatabaseTemplate( const vec<TAG>* p_vecOfTaggedRpints )
  : mp_taggedRpints( p_vecOfTaggedRpints ),
    m_ownRpints( false ),
    m_to( 0 ),
    m_hits( 0 ),
    m_misses( 0 ),
    m_missDistance( 0 )
{
}

template <class TAG>
template <class KmerPathVector>
void KmerPathDatabaseTemplate<TAG>::ConstructFromKmerPaths( const KmerPathVector& fwdPaths, 
							    const KmerPathVector& revPaths )
{
  vec<TAG>* p_taggedRpints = new vec<TAG>();
  unsigned int numSegs = 0;
  for ( int pass = 1; pass <= 2; ++pass ) {
    for ( int rc = 0; rc < 2; ++rc ) {
      const KmerPathVector& paths = ( rc==0 ? fwdPaths : revPaths );
      for ( int pathIdx = 0; pathIdx < (int) paths.size(); ++pathIdx )
        for ( int segIdx = 0; segIdx < paths[pathIdx].NSegments(); ++segIdx ) {
          const KmerPathInterval& I = paths[pathIdx].Segment(segIdx);
          if ( I.isSeq() ) {
            // In the first pass, just count the number of segments for later reserve().
            if ( pass == 1 ) {
              ++numSegs;
            }
            // In the second pass, validate the data and store it in the vector.
            else {
              ForceAssertLe( (unsigned int)I.Length( ), TAG::LENGTH_MAX );
              ForceAssertLe( (unsigned int)segIdx,      TAG::POSITION_MAX );
              int idx = ( rc==0 ? pathIdx : -pathIdx-1 );
              p_taggedRpints->push_back( TAG( I.Start(),
                                                       I.Length(),
                                                       idx, segIdx ) );
            }
          }
        }
    }
    if ( pass == 1 ) {
      p_taggedRpints->reserve( numSegs );
    }
  }
  
  Prepare( *p_taggedRpints );
  mp_taggedRpints = p_taggedRpints;
}

template <class TAG>
KmerPathDatabaseTemplate<TAG>::KmerPathDatabaseTemplate( const vecKmerPath& fwdPaths,
					 const vecKmerPath& revPaths )
  : mp_taggedRpints( 0 ),
    m_ownRpints( true ),
    m_to( 0 ),
    m_hits( 0 ),
    m_misses( 0 ),
    m_missDistance( 0 )
{
  this->ConstructFromKmerPaths( fwdPaths, revPaths );
}

template <class TAG>
KmerPathDatabaseTemplate<TAG>::KmerPathDatabaseTemplate( const vec<KmerPath>& fwdPaths,
					 const vec<KmerPath>& revPaths )
  : mp_taggedRpints( 0 ),
    m_ownRpints( true ),
    m_to( 0 ),
    m_hits( 0 ),
    m_misses( 0 ),
    m_missDistance( 0 )
{
  this->ConstructFromKmerPaths( fwdPaths, revPaths );
}

template <class TAG>
KmerPathDatabaseTemplate<TAG>& KmerPathDatabaseTemplate<TAG>::operator= ( const KmerPathDatabaseTemplate<TAG>& other )
{
  this->CopyFrom( other );
  return *this;
}

template <class TAG>
KmerPathDatabaseTemplate<TAG>::~KmerPathDatabaseTemplate() 
{
  if ( m_ownRpints && mp_taggedRpints )
    delete mp_taggedRpints;
}

template <class TAG>
void KmerPathDatabaseTemplate<TAG>::Contains( kmer_id_t kmerId, 
				      vec<path_interval_id_t>& answer, 
				      bool append ) const
{
  if ( ! append )
    answer.clear( );

  // Find index of first interval that starts with a kmer greater than the given kmer.
  TAG target;
  target.Set( kmerId, 0, 0, 0 );
  
  // Convenience:
  const vec<TAG>& segs = *mp_taggedRpints;

#define CACHE_TO 0

#if CACHE_TO
  // Here we use the prior value of "m_to" to restrict the area of our
  // search.  After this set of statements, "m_to" will hold the index
  // of the first interval starting with a kmer greater than "kmerId".
  
  if ( segs[m_to].Start() > kmerId )
  {
    const unsigned int lowCheckRange = 32;
    unsigned int lowCheck = m_to;
    if ( lowCheck >= lowCheckRange )
      lowCheck -= lowCheckRange;
    else
      lowCheck = 0;

    if ( segs[lowCheck].Start() <= kmerId )
    {
      ++m_hits;
      m_to = distance( segs.begin(),
                       upper_bound( segs.begin() + lowCheck, 
                                    segs.begin() + m_to + 1, 
                                    target ) );
    }
    else
    {
      ++m_misses;
      longlong to = distance( segs.begin(),
                              upper_bound( segs.begin(),
                                           segs.end(), 
                                           target ) );
      m_missDistance += abs( to - m_to );
      m_to = to;
    }
  }
  
  else // segs[m_to].Start() <= kmerId 
  {
    const unsigned int highCheckRange = 32;
    unsigned int highCheck = m_to + highCheckRange;
    if ( highCheck >= segs.size() )
      highCheck = segs.size()-1;

    if ( segs[highCheck].Start() > kmerId )
    {
      ++m_hits;
      m_to = distance( segs.begin(),
                       upper_bound( segs.begin() + m_to, 
                                    segs.begin() + highCheck + 1, 
                                    target ) );
    }
    else
    {
      ++m_misses;
      longlong to = distance( segs.begin(),
                              upper_bound( segs.begin(),
                                           segs.end(),
                                           target ) );
      m_missDistance += abs( to - m_to );
      m_to = to;
    }
  }
    
#else
  // Simple way with no caching.
  m_to = distance( segs.begin(),
                   upper_bound( segs.begin(),
                                segs.end(),
                                target ) );
#endif

  // If it's the first interval, there are no intervals containing the given kmer.
  if ( m_to == 0 ) return;
  
  // Otherwise, decrement "to" so it points to the last interval that
  // starts with a kmer less than or equal to the given kmer.
  --m_to;
  
  // Use lookback to find the first interval that contains the "to"
  // interval's last kmer (which may be further back than we really
  // need to go).
  unsigned int from = m_to - segs[m_to].Lookback( );

  // Copy all the indices between "from" and "to" that contain the given kmer.
  for ( unsigned int i = from; i <= m_to; i++ )
    if ( kmerId >= segs[i].Start( ) &&
         kmerId <= segs[i].Stop( ) )
      answer.push_back(i);
}




template <class TAG>
void KmerPathDatabaseTemplate<TAG>::GetHighFrequencyKmers( vec<kmer_id_t>& hiFreqKmers,
						   const copy_num_t minFreq ) const
{
  map<kmer_id_t,copy_num_t> kmerCount;

  const size_t passSize = 1000000;
  cout << "Finding high frequency kmers in " 
       << mp_taggedRpints->size()/passSize
       << " passes." << endl;

  size_t dot = 0;
  for ( typename vec<TAG>::const_iterator newIter = mp_taggedRpints->begin();
        newIter != mp_taggedRpints->end(); ++newIter )
  {
    if ( dot % passSize == 0 )
      Dot( cout, dot / passSize );

    ++dot;

    while ( ! kmerCount.empty() )
    {
      if ( kmerCount.begin()->first < newIter->Start() )
      {
        if ( kmerCount.begin()->second >= minFreq )
          hiFreqKmers.push_back( kmerCount.begin()->first );
        kmerCount.erase( kmerCount.begin() );
      }
      else
        break;
    }

    kmer_id_t kmer = newIter->Start();
    map<kmer_id_t,copy_num_t>::iterator countIter = kmerCount.find( kmer );
    while ( countIter != kmerCount.end() &&
            countIter->first <= newIter->Stop() )
    {
      ForceAssertEq( kmer, countIter->first );
      ++(countIter->second);
      ++kmer;
      ++countIter;
    }

    for ( ; kmer <= newIter->Stop(); ++kmer )
      kmerCount.insert( kmerCount.end(), make_pair( kmer, 0 ) );
  }

  map<kmer_id_t,copy_num_t>::iterator countIter = kmerCount.begin();
  for ( ; countIter != kmerCount.end(); ++countIter )
    if ( countIter->second >= minFreq )
      hiFreqKmers.push_back( countIter->first );
}

template <class TAG>
void KmerPathDatabaseTemplate<TAG>::Write( const String& filename ) const
{
  BinaryWrite2( filename, *mp_taggedRpints );
}

template <class TAG>
void KmerPathDatabaseTemplate<TAG>::Read( const String& filename )
{
  if ( m_ownRpints )
    delete mp_taggedRpints;

  vec<TAG>* p_taggedRpints = new vec<TAG>();

  BinaryRead2( filename, *p_taggedRpints );

  mp_taggedRpints = p_taggedRpints;
  m_ownRpints = true;
  m_to = 0;
}


template class KmerPathDatabaseTemplate<tagged_rpint>;
template class KmerPathDatabaseTemplate<big_tagged_rpint>;
