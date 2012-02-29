/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PAIRWISE_ALIGNERS_MATCHLIST_H
#define PAIRWISE_ALIGNERS_MATCHLIST_H

#include "Vec.h"
#include "pairwise_aligners/MaxMutmerFromMer.h"

// Class that tracks info relating to a match between some implicit
// sequence and the given sequence (id2).

class Match {
 public:
  Match()
  {}

  Match( const int id2, const bool rc, const int pos1, const int pos2, const int len )
    : m_id2_rc( rc ? -id2-1 : id2 ), m_begin1( pos1 ), m_end1( pos1+len), m_offset( pos1-pos2 )
  {}
  
  int GetId2() const    { return ( m_id2_rc < 0 ? -m_id2_rc-1 : m_id2_rc ); }
  int GetRc() const     { return ( m_id2_rc < 0 ); }
  int GetPos1() const   { return m_begin1; }
  int GetPos2() const   { return m_begin1 - m_offset; }
  int GetOffset() const { return m_offset; }
  int GetLen() const    { return m_end1 - m_begin1; }

  // Presuming both Matches have the same id1, does this contain the other?
  bool Contains( const Match& other ) const;

  // Presuming both Matches have the same id1, id2, rc, and offset, does this contain the other?
  bool QuickContains( const Match& other ) const;

  bool operator< ( const Match& other ) const
  {
    return ( m_id2_rc < other.m_id2_rc ||
             m_id2_rc == other.m_id2_rc && m_offset < other.m_offset );
  }

  bool operator> ( const Match& other ) const
  {
    return ( m_id2_rc > other.m_id2_rc ||
             m_id2_rc == other.m_id2_rc && m_offset > other.m_offset );
  }

 private:
  // TODO: Potentially dangerous truncation of IDs
  int m_id2_rc;
  int m_begin1;
  int m_end1;
  int m_offset;
};


inline
bool Match::Contains( const Match& other ) const
{
  return ( other.m_id2_rc == m_id2_rc &&
           other.m_offset == m_offset &&
           other.m_begin1 >= m_begin1 &&
           other.m_end1 <= m_end1 );
}


inline
bool Match::QuickContains( const Match& other ) const
{
  return ( other.m_begin1 >= m_begin1 &&
           other.m_end1 <= m_end1 );
}


// Class to track matches among a set of sequences.

// The basic idea is to keep matches in two vectors: one sorted (quick to
// search) and one unsorted (quick to insert).  
//
// When ProcessMatchingKmers() is called, we check for a pre-existing match that
// subsumes those matching kmers.  We check the sorted list first [O(log n) in
// the number of sorted matches].  If that fails, we check the unsorted list
// [O(n) in the number of unsorted matches].  If that fails, then it's a new
// match, and it gets appended to the unsorted list.
//
// If the unsorted list gets larger than sortedBatchSize, then the unordered
// list gets sorted, and the two lists are merged into a new, larger ordered
// list, and the unordered list is emptied.
//
// Matches between sequences id1 and id2 are stored under either id1 OR id2.
// For any given sequences id1 and id2, the choice is deterministic (so we only
// have to check one set of lists) but arbitrary.  This means that getting the
// list of matches for some id will not necessarily (and in fact will rarely)
// return ALL the matches involving that id.  MatchList is intended for
// processing all matches among a set of sequences, not finding all the matches
// of a set of sequences to a given sequence.

class MatchList {
 public:
  MatchList( const int numReads, const unsigned int sortedBatchSize = 256 )
    : m_sortedMatches( numReads ),
      m_unsortedMatches( numReads ),
      m_sortedBatchSize( sortedBatchSize )
  {}
  
  void ProcessMatchingKmers( int pos1, int pos2, int len,
                             int id1, int id2,
                             const basevector* pSeq1, const basevector* pSeq2 );

  void GetMatches( const int id1, vec<Match>& matches ) const;

  void SetSortedBatchSize( const unsigned int size ) { m_sortedBatchSize = size; }

  void GetMatchCounts( vec<int>& matchCounts ) const;

 private:
  bool FindMatchInSorted( const int id1, const Match& newMatch ) const;

  bool FindMatchInUnsorted( const int id1, const Match& newMatch ) const;

  void AddMatch( const int id1, const Match& newMatch );

  void UnsortedToSorted( const int id1 );

  vec< vec<Match> > m_sortedMatches;
  vec< vec<Match> > m_unsortedMatches;
  unsigned int m_sortedBatchSize;
};

#endif
